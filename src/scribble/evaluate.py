#!/usr/bin/env python

from pathlib import Path

def classify_cluster(row, thresholds):

    n = row["n_cells"]
    stability = row["mean_stability"]
    entropy = row["sample_entropy"]
    fraction = row.get("fraction_of_dataset", 0)

    reassignment_fraction = row.get(
        "dominant_reassignment_fraction",
        0.0,
    )

    if n < thresholds["min_cells"]:
        return "keep", "small_cluster", "low"

    # Strong evidence that the cluster is not robustly separated
    if reassignment_fraction >= thresholds["high_reassignment"]:
        return (
            "merge",
            "high_reassignment_instability",
            "high",
        )

    if stability < thresholds["low_stability"]:
        return (
            "trajectory",
            "low_stability_continuum",
            "medium",
        )

    if fraction >= thresholds["large_cluster_fraction"]:
        return (
            "subset",
            "large_cluster_requires_inspection",
            "high",
        )

    if entropy < thresholds["low_entropy"]:
        return (
            "flag_bias",
            "sample_specific_cluster",
            "medium",
        )

    return "keep", "well_defined_cluster", "low"


def find_merge_candidates(
    adata,
    decision_df,
    cluster_key="leiden",
    connectivity_percentile=90
):

    import numpy as np

    conn = adata.obsp["connectivities"]

    candidate_clusters = (
        decision_df.loc[
            decision_df["action"] == "merge",
            "cluster"
        ]
        .astype(str)
        .tolist()
    )


    pair_scores = []

    for i, clust_a in enumerate(candidate_clusters):

        cells_a = np.where(
            adata.obs[cluster_key].astype(str) == clust_a
        )[0]

        if len(cells_a) == 0:
            continue

        for clust_b in candidate_clusters[i + 1:]:

            cells_b = np.where(
                adata.obs[cluster_key].astype(str) == clust_b
            )[0]

            if len(cells_b) == 0:
                continue

            edges_ab = conn[cells_a][:, cells_b].sum()

            score = (
                edges_ab /
                np.sqrt(
                    len(cells_a) * len(cells_b)
                )
            )

            pair_scores.append(
                (clust_a, clust_b, float(score))
            )

    if len(pair_scores) == 0:
        return []

    scores = [x[2] for x in pair_scores]

    cutoff = np.percentile(
        scores,
        connectivity_percentile
    )

    print(
        f"Merge connectivity cutoff "
        f"({connectivity_percentile}th percentile): "
        f"{cutoff:.4f}"
    )

    merge_pairs = []

    for clust_a, clust_b, score in pair_scores:

        print(
            f"{clust_a} <-> {clust_b}: "
            f"{score:.4f}"
        )

        if score >= cutoff:

            merge_pairs.append(
                (clust_a, clust_b)
            )

            print(
                f"MERGE "
                f"{clust_a} <-> {clust_b} "
                f"(score={score:.4f})"
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
            instability_burden = (1.0 - df["mean_stability"]).mean()
            weighted_mean_reassignment = np.average(
                df["dominant_reassignment_fraction"]
                .fillna(0),
                weights=df["n_cells"]
            )
            n_merge_candidates = (
                df["dominant_reassignment_fraction"]
                .fillna(0)
                .gt(args.high_reassignment)
                .sum()
            )
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
                - instability_burden
                - weighted_mean_reassignment
                - 0.20 * n_merge_candidates
                - 0.1 * cluster_penalty
            )

            comparison_rows.append({
                "file": f.name,
                "path": str(f),
                "score": score,
                "mean_stability": mean_stability,
                "mean_entropy": mean_entropy,
                "instability_burden": instability_burden,
                "weighted_mean_reassignment": weighted_mean_reassignment,
                "n_merge_candidates": n_merge_candidates,
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
            "instability_burden",
            "weighted_mean_reassignment",
            "n_merge_candidates",
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
    total_cells = df["n_cells"].sum()

    # ----------------------------
    # Thresholds
    # ----------------------------
    thresholds = {
        "min_cells": args.min_cells,
        "large_cluster_fraction": args.large_cluster_fraction,
        "low_stability": args.low_stability,
        "high_stability": args.high_stability,
        "low_entropy": args.low_entropy,
        "high_reassignment": args.high_reassignment,
    }

    # ----------------------------
    # Classification
    # ----------------------------
    decisions = []

    for _, row in df.iterrows():

        row = row.copy()
        row["fraction_of_dataset"] = (row["n_cells"] / total_cells)
        action, reason, priority = classify_cluster(row,thresholds)


        dominant = row.get(
            "dominant_reassignment",
            None,
        )

        reassign = row.get(
            "dominant_reassignment_fraction",
            None,
        )

        detail = (
            f"n={row['n_cells']}; "
            f"stability={row['mean_stability']:.2f}; "
            f"entropy={row['sample_entropy']:.2f}"
        )

        if not pd.isna(dominant):
            detail += (
                f"; dominant={dominant}; "
                f"reassign={reassign:.1%}"
            )

        suggested_merge_target = ""

        if action == "merge":
            suggested_merge_target = row.get(
                "dominant_reassignment",
                ""
            )

        decisions.append({
            "cluster": row["cluster"],
            "n_cells": row["n_cells"],
            "fraction_of_dataset": row["fraction_of_dataset"],
            "action": action,
            "reason": reason,
            "detail": detail,
            "priority": priority,
            "suggested_merge_target": suggested_merge_target,
        })

    out_df = pd.DataFrame(decisions)

    # ----------------------------
    # Merge candidates
    # ----------------------------
    cluster_col = "leiden"

    merge_pairs = []

    for _, row in out_df.iterrows():

        if (
            row["action"] == "merge"
            and pd.notna(row["suggested_merge_target"])
            and str(row["suggested_merge_target"]) != ""
        ):

            target = row["suggested_merge_target"]

            if pd.notna(target):

                try:
                    target = str(int(float(target)))
                except Exception:
                    target = str(target)

                merge_pairs.append(
                    (
                        str(row["cluster"]),
                        target,
                    )
                )


    graph_merge_pairs = find_merge_candidates(
        adata=adata,
        decision_df=out_df,
        cluster_key=cluster_col,
        connectivity_percentile=args.merge_percentile
    )

    merge_pairs.extend(graph_merge_pairs)

    merge_pairs = list(
        {
            tuple(sorted(pair))
            for pair in merge_pairs
        }
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

    out_df["cluster"] = (
        out_df["cluster"]
        .astype(str)
    )

    for idx, group in enumerate(merge_groups):

        group_label = f"group_{idx+1}"

        group = {
            str(int(float(x)))
            if str(x).replace(".", "", 1).isdigit()
            else str(x)
            for x in group
        }

        print(
            f"[ASSIGN MERGE GROUP] "
            f"{group_label} -> {sorted(group)}"
        )

        out_df.loc[
            out_df["cluster"].astype(str).isin(group),
            "merge_group"
        ] = group_label

    print(out_df["action"].value_counts())

    print("\n[FINAL MERGE GROUPS]")

    print(
        out_df.loc[
            out_df["merge_group"] != "",
            [
                "cluster",
                "action",
                "suggested_merge_target",
                "merge_group",
            ]
        ]
    )

    print(f"Saving decisions → {output_file}")
    out_df.to_csv(output_file, sep="\t", index=False)
