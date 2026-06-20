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

            if (
                size_ratio < thresholds["merge_size_ratio"] and
                abs(s_a - s_b) < thresholds["merge_stability_tol"] and
                abs(e_a - e_b) < thresholds["merge_entropy_tol"]
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

    input_file = Path(args.input)
    output_file = input_file.with_name(f"{input_file.stem}_decisions.tsv")

    print(f"Loading cluster summary: {input_file}")
    df = pd.read_csv(input_file, sep="\t")

    # ----------------------------
    # Thresholds (configurable)
    # ----------------------------
    thresholds = {
        "min_cells": args.min_cells,          # 200
        "large_cells": args.large_cells,      # 800
        "low_stability": args.low_stability,  # 0.75
        "high_stability": args.high_stability,# 0.95
        "low_entropy": args.low_entropy       # 0.5
        "merge_size_ratio": args.merge_size_ratio,
        "merge_stability_tol": args.merge_stability_tol,
        "merge_entropy_tol": args.merge_entropy_tol
    }


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
            "priority": priority
        })

    out_df = pd.DataFrame(decisions)

    merge_pairs = find_merge_candidates(out_df, thresholds)

    # Convert to groups
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
