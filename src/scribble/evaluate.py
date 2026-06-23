#!/usr/bin/env python

from pathlib import Path

def classify_cluster(row, thresholds):
    n = row["n_cells"]
    stability = row["mean_stability"]
    entropy = row["sample_entropy"]

    if n < thresholds["min_cells"]:
        return "keep", "small_cluster", "low"

    if stability < thresholds["low_stability"]:
        return "trajectory", "low_stability_continuum", "medium"

    if n > thresholds["large_cells"] and stability < thresholds["high_stability"]:
        return "subset", "heterogeneous_large_cluster", "high"

    if entropy < thresholds["low_entropy"]:
        return "flag_bias", "sample_specific_cluster", "medium"

    return "keep", "well_defined_cluster", "low"


def find_merge_candidates(decision_df, thresholds):
    merge_pairs = []

    subset_df = decision_df[decision_df["action"] == "subset"]

    for i, row_a in subset_df.iterrows():
        for j, row_b in subset_df.iterrows():
            if row_a["cluster"] >= row_b["cluster"]:
                continue

            if row_a["action"] != "subset":
                continue
            if row_b["action"] != "subset":
                continue

            # Parse values from detail string
            def parse(detail, key):
                return float(detail.split(key + "=")[1].split(";")[0])

            n_a = parse(row_a["detail"], "n")
            n_b = parse(row_b["detail"], "n")
            s_a = parse(row_a["detail"], "stability")
            s_b = parse(row_b["detail"], "stability")
            e_a = parse(row_a["detail"], "entropy")
            e_b = parse(row_b["detail"], "entropy")

            size_ratio = max(n_a, n_b) / min(n_a, n_b)

            low_entropy_cutoff = 0.5

            if (
                size_ratio < thresholds["merge_size_ratio"]
                and abs(s_a - s_b) < thresholds["merge_stability_tol"]
                and (
                    abs(e_a - e_b) < thresholds["merge_entropy_tol"]
                    or (
                        (e_a < low_entropy_cutoff and e_b > low_entropy_cutoff)
                        or (e_b < low_entropy_cutoff and e_a > low_entropy_cutoff)
                    )
                )
            ):
                merge_pairs.append((row_a["cluster"], row_b["cluster"]))

    return merge_pairs


def run_evaluate(args):
    import scanpy as sc
    import numpy as np
    import pandas as pd
    import random
    from scribble.import_data import setup_environment

    PROJECT_DIR = Path(args.project_dir)
    PLOT_DIR = PROJECT_DIR / "scribble/plots"
    TABLE_DIR = PROJECT_DIR / "scribble/tables"

    setup_environment(sc, np, random, PLOT_DIR)

    # ----------------------------
    # Handle input(s)
    # ----------------------------
    input_files = args.input if isinstance(args.input, list) else [args.input]
    input_files = [Path(f) for f in input_files]

    comparison_rows = []

    # ----------------------------
    # Evaluate multiple runs
    # ----------------------------
    if len(input_files) > 1:

        print("Evaluating multiple cluster summaries...\n")

        for f in input_files:
            df = pd.read_csv(f, sep="\t")

            mean_stability = df["mean_stability"].mean()
            mean_entropy = df["sample_entropy"].mean()
            low_stability_frac = (df["mean_stability"] < 0.7).mean()
            n_clusters = df.shape[0]

            # ---- cluster count penalty (soft constraint) ----
            if n_clusters < 10:
                cluster_penalty = 10 - n_clusters
            elif n_clusters > 30:
                cluster_penalty = n_clusters - 30
            else:
                cluster_penalty = 0

            score = (
                mean_stability
                + 0.5 * mean_entropy
                - 2.0 * low_stability_frac
                - 0.1 * cluster_penalty
            )

            comparison_rows.append({
                "file": f.name,
                "path": str(f),
                "score": score,
                "mean_stability": mean_stability,
                "mean_entropy": mean_entropy,
                "low_stability_fraction": low_stability_frac,
                "n_clusters": n_clusters
            })

        comparison_df = pd.DataFrame(comparison_rows)
        comparison_df = comparison_df.sort_values("score", ascending=False)

        # Save comparison table
        comparison_file = TABLE_DIR / "clustering_comparison.tsv"
        comparison_df.to_csv(comparison_file, sep="\t", index=False)

        print("Comparison summary:")
        print(comparison_df[[
            "file",
            "score",
            "mean_stability",
            "mean_entropy",
            "low_stability_fraction",
            "n_clusters"
        ]])

        best_row = comparison_df.iloc[0]
        best_file = Path(best_row["path"])

        print(f"\nSelected best clustering → {best_file.name}")
        print(f"Score: {best_row['score']:.4f}\n")

        input_file = best_file

    else:
        input_file = input_files[0]

    # ----------------------------
    # Continue existing behaviour
    # ----------------------------
    output_file = input_file.with_name(f"{input_file.stem}_decisions.tsv")

    print(f"Loading cluster summary: {input_file}")

    df = pd.read_csv(input_file, sep="\t")

    # ----------------------------
    # Thresholds
    # ----------------------------
    thresholds = {
        "min_cells": args.min_cells,
        "large_cells": args.large_cells,
        "low_stability": args.low_stability,
        "high_stability": args.high_stability,
        "low_entropy": args.low_entropy,
        "merge_size_ratio": args.merge_size_ratio,
        "merge_stability_tol": args.merge_stability_tol,
        "merge_entropy_tol": args.merge_entropy_tol,
    }

    # ----------------------------
    # Classification
    # ----------------------------
    decisions = []

    for _, row in df.iterrows():
        action, reason, priority = classify_cluster(row, thresholds)

        detail = (
            f"n={row['n_cells']}; "
            f"stability={row['mean_stability']:.2f}; "
            f"entropy={row['sample_entropy']:.2f}"
        )

        decisions.append({
            "cluster": row["cluster"],
            "action": action,
            "reason": reason,
            "detail": detail,
            "priority": priority,
        })

    out_df = pd.DataFrame(decisions)

    # ----------------------------
    # Merge candidates
    # ----------------------------
    merge_pairs = find_merge_candidates(out_df, thresholds)

    merge_groups = []
    visited = set()

    for a, b in merge_pairs:
        if a not in visited and b not in visited:
            merge_groups.append({a, b})
            visited.update([a, b])
        else:
            for group in merge_groups:
                if a in group or b in group:
                    group.update([a, b])
                    visited.update([a, b])

    out_df["merge_group"] = ""

    for idx, group in enumerate(merge_groups):
        group_label = f"group_{idx+1}"
        for cl in group:
            out_df.loc[out_df["cluster"] == cl, "merge_group"] = group_label

    print(out_df["action"].value_counts())

    print(f"Saving decisions → {output_file}")
    out_df.to_csv(output_file, sep="\t", index=False)
