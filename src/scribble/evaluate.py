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


def find_merge_candidates(
    adata,
    decision_df,
    cluster_key="leiden",
    min_connectivity=0.10
):
    """
    Find candidate merge pairs using normalized graph connectivity.

    Connectivity score is the mean of:

        edges(A->B) / total_edges(A)
        edges(B->A) / total_edges(B)

    Values range approximately from 0 to 1.

    Only clusters marked for 'subset' are evaluated.
    """

    import numpy as np

    conn = adata.obsp["connectivities"]

    subset_clusters = (
        decision_df.loc[
            decision_df["action"] == "subset",
            "cluster"
        ]
        .astype(str)
        .tolist()
    )

    merge_pairs = []

    for i, clust_a in enumerate(subset_clusters):

        cells_a = np.where(
            adata.obs[cluster_key].astype(str) == clust_a
        )[0]

        if len(cells_a) == 0:
            continue

        total_edges_a = conn[cells_a].sum()

        for clust_b in subset_clusters[i + 1:]:

            cells_b = np.where(
                adata.obs[cluster_key].astype(str) == clust_b
            )[0]

            if len(cells_b) == 0:
                continue

            total_edges_b = conn[cells_b].sum()

            if total_edges_a == 0 or total_edges_b == 0:
                continue

            edges_ab = conn[cells_a][:, cells_b].sum()
            edges_ba = conn[cells_b][:, cells_a].sum()

            frac_ab = edges_ab / total_edges_a
            frac_ba = edges_ba / total_edges_b

            connectivity_score = (
                frac_ab + frac_ba
            ) / 2

            if connectivity_score >= min_connectivity:

                merge_pairs.append(
                    (clust_a, clust_b)
                )

                print(
                    f"MERGE "
                    f"{clust_a} <-> {clust_b} "
                    f"(graph connectivity="
                    f"{connectivity_score:.3f})"
                )

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
    adata_file = (PROJECT_DIR / "scribble/adata" / input_file.name.replace("_cluster_summary.tsv", "_clustered.h5ad"))
    adata = sc.read(adata_file)

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
    cluster_col = "leiden"

    merge_pairs = find_merge_candidates(
        adata=adata,
        decision_df=out_df,
        cluster_key=cluster_col,
        min_connectivity=args.merge_connectivity
    )

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
