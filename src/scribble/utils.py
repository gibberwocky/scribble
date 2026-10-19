#!/usr/bin/env python

import numpy as np
import pandas as pd
import scanpy as sc
from scipy.stats import entropy


def compute_cluster_stability(
    adata,
    resolution,
    n_repeats,
    cluster_key="leiden",
    random_state_base=0,
):
    """
    Compute per-cell cluster stability using repeated Leiden clustering.

    Returns
    -------
    dict
        {
            "cell_stability": ndarray (n_cells,),
            "aligned_assignments": ndarray (n_repeats, n_cells),
            "reference_labels": ndarray (n_cells,)
        }
    """

    assignments = []

    for i in range(n_repeats):

        temp_key = f"__stability_{i}"

        sc.tl.leiden(
            adata,
            resolution=resolution,
            key_added=temp_key,
            flavor="igraph",
            directed=False,
            n_iterations=2,
            random_state=random_state_base + i,
        )

        assignments.append(
            adata.obs[temp_key].astype(str).values.copy()
        )

    assignments = np.array(assignments)

    ref_labels = assignments[0]
    aligned_assignments = [ref_labels]

    for i in range(1, assignments.shape[0]):

        current = assignments[i]

        contingency = pd.crosstab(
            pd.Series(current, name="current"),
            pd.Series(ref_labels, name="reference"),
        )

        mapping = contingency.idxmax(axis=1).to_dict()

        aligned = np.array(
            [mapping.get(label, label) for label in current]
        )

        aligned_assignments.append(aligned)

    aligned_assignments = np.array(aligned_assignments)

    stability = np.array([
        pd.Series(aligned_assignments[:, j]).value_counts().max()
        / n_repeats
        for j in range(aligned_assignments.shape[1])
    ])

    for i in range(n_repeats):
        del adata.obs[f"__stability_{i}"]

    return {
        "cell_stability": stability,
        "aligned_assignments": aligned_assignments,
        "reference_labels": ref_labels,
    }

def summarise_cluster_reassignments(
    final_labels,
    aligned_assignments,
):
    """
    Summarise how cells move between aligned clusters
    across repeated clustering runs.

    Returns
    -------
    pd.DataFrame

    Columns:
        cluster
        destination
        n_votes
        proportion
    """

    final_labels = np.asarray(final_labels)

    rows = []

    unique_clusters = sorted(pd.unique(final_labels))

    for cluster in unique_clusters:

        mask = final_labels == cluster

        cluster_assignments = aligned_assignments[:, mask]

        if cluster_assignments.size == 0:
            continue

        values = pd.Series(
            cluster_assignments.flatten()
        ).value_counts()

        total = values.sum()

        for destination, count in values.items():

            rows.append({
                "cluster": cluster,
                "destination": destination,
                "n_votes": int(count),
                "proportion": count / total,
            })

    return pd.DataFrame(rows)


def build_reassignment_strings(
    reassignment_df,
    top_n=5,
):
    """
    Convert reassignment table into a compact summary string.

    Example:
        R3:72%;R1:18%;6:10%
    """

    summaries = {}

    for cluster, group in reassignment_df.groupby("cluster"):

        group = group.sort_values(
            "proportion",
            ascending=False,
        ).head(top_n)

        summaries[cluster] = ";".join(
            f"{row.destination}:{row.proportion:.1%}"
            for _, row in group.iterrows()
        )

    return summaries


def compute_cluster_summary(
    adata,
    cluster_key,
    stability_key="cluster_stability",
    sample_key="sample",
    reassignment_df=None,
    dominant_reassignment_threshold=0.05,
    reassignment_display_threshold=0.01,
    merge_candidate_threshold=0.30,
):
    """
    Generate cluster summary table.

    Parameters
    ----------
    adata : AnnData
    cluster_key : str
        Column containing cluster labels.
    stability_key : str
        Column containing per-cell stability values.
    sample_key : str
        Column containing sample IDs.
    reassignment_df : pd.DataFrame, optional
        Output from summarise_cluster_reassignments().
    dominant_reassignment_threshold: float
        Minimal proportion for dominant_reassignment
    reassignment_display_threshold: float
        Minimal reassignment threshold to display

    Returns
    -------
    pd.DataFrame
    """

    summary = (
        adata.obs
        .groupby(cluster_key, observed=True)
        .agg(
            n_cells=(cluster_key, "size"),
            mean_stability=(stability_key, "mean"),
            median_stability=(stability_key, "median"),
        )
        .reset_index()
        .rename(columns={cluster_key: "cluster"})
    )

    # --------------------------------------------------
    # Sample entropy
    # --------------------------------------------------
    if sample_key in adata.obs.columns:

        sample_counts = (
            adata.obs
            .groupby(
                [cluster_key, sample_key],
                observed=True,
            )
            .size()
            .unstack(fill_value=0)
        )

        summary["sample_entropy"] = (
            sample_counts.apply(
                lambda row: entropy(row / row.sum()),
                axis=1,
            ).values
        )

    else:

        summary["sample_entropy"] = np.nan

    # --------------------------------------------------
    # Reassignment diagnostics
    # --------------------------------------------------
    summary["dominant_reassignment"] = np.nan
    summary["dominant_reassignment_fraction"] = np.nan

    if reassignment_df is not None:

        summary_strings = {}

        dominant_targets = {}
        dominant_fractions = {}

        for cluster, group in reassignment_df.groupby("cluster"):

            group = group.sort_values(
                "proportion",
                ascending=False
            )

            # Compact summary string for reports
            group = group[
                group["proportion"] >= 0.01
            ]

            display_group = group[
                group["proportion"]
                >= reassignment_display_threshold
            ]

            summary_strings[str(cluster)] = ";".join(
                f"{row.destination}:{row.proportion:.1%}"
                for _, row in display_group.iterrows()
            )

            # Highest-probability alternative destination
            alternatives = group[
                group["destination"].astype(str)
                != str(cluster)
            ]

            if len(alternatives):

                best = alternatives.iloc[0]

                if (
                    best["proportion"]
                    >= dominant_reassignment_threshold
                ):

                    dominant_targets[str(cluster)] = str(
                        best["destination"]
                    )

                    dominant_fractions[str(cluster)] = float(
                        best["proportion"]
                    )


        summary["reassignment_summary"] = (
            summary["cluster"]
            .astype(str)
            .map(summary_strings)
            .fillna("")
        )

        summary["dominant_reassignment"] = (
            summary["cluster"]
            .astype(str)
            .map(dominant_targets)
        )

        summary["dominant_reassignment_fraction"] = (
            summary["cluster"]
            .astype(str)
            .map(dominant_fractions)
        )

    summary["merge_candidate"] = (
        summary["dominant_reassignment_fraction"]
        .fillna(0)
        .ge(merge_candidate_threshold)
    )

    summary["reassignment_class"] = np.select(
        [
            summary[
                "dominant_reassignment_fraction"
            ].fillna(0) >= merge_candidate_threshold,

            summary[
                "dominant_reassignment_fraction"
            ].fillna(0) >= 0.10,
        ],
        [
            "strong",
            "moderate",
        ],
        default="weak",
    )



    return summary
