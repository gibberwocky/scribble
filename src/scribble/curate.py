#!/usr/bin/env python

from pathlib import Path


def parse_merge_groups(merge_string):

    groups = []

    if not merge_string:
        return groups

    for group in merge_string.split(","):

        clusters = [
            x.strip()
            for x in group.split("+")
            if x.strip()
        ]

        if len(clusters) >= 2:
            groups.append(clusters)

    return groups


def parse_refine_requests(refine_string):

    if not refine_string:
        return []

    return [
        x.strip()
        for x in refine_string.split(",")
        if x.strip()
    ]

def run_curate(args):

    import scanpy as sc
    import pandas as pd
    import numpy as np
    import random

    from scribble.import_data import setup_environment
    from scribble.refine import (
        run_clustering,
        restore_counts
    )

    PROJECT_DIR = Path(args.project_dir)

    PLOT_DIR = PROJECT_DIR / "scribble/plots"
    TABLE_DIR = PROJECT_DIR / "scribble/tables"
    ADATA_DIR = PROJECT_DIR / "scribble/adata"

    setup_environment(sc, np, random, PLOT_DIR)
    input_file = Path(args.input)
    curation_file = Path(args.curation)

    print(f"Loading AnnData → {input_file}")

    adata = sc.read(input_file)

    merge_groups = parse_merge_groups(
        args.merge
    )

    refine_requests = parse_refine_requests(
        args.refine
    )

    print("\n[MERGE GROUPS]")
    print(merge_groups)

    print("\n[REFINE REQUESTS]")
    print(refine_requests)

    adata.obs["curated_cluster"] = (
        adata.obs["refine_label"]
        .astype(str)
    )

    for group in merge_groups:

        target = group[0]

        for cluster in group[1:]:

            print(
                f"[MERGE] "
                f"{cluster} -> {target}"
            )

            adata.obs.loc[
                adata.obs["curated_cluster"]
                == cluster,
                "curated_cluster"
            ] = target

    for cluster in refine_requests:

        mask = (
            adata.obs["curated_cluster"]
            .astype(str)
            == cluster
        )

        adata_sub = adata[mask].copy()

        print(
            f"[REFINE] "
            f"{cluster} "
            f"({adata_sub.n_obs} cells)"
        )

        adata_sub = run_refinement_clustering(
            adata_sub,
            args,
        )

    cluster_map = {}

    for i, cl in enumerate(
        sorted(
            adata_sub.obs[
                "leiden_refined"
            ].astype(str).unique()
        ),
        start=1
    ):

        cluster_map[cl] = (
            f"{cluster}_S{i}"
        )

    refined = (
        adata_sub.obs["leiden_refined"]
        .astype(str)
        .map(cluster_map)
    )

    adata.obs.loc[
        adata_sub.obs_names,
        "curated_cluster"
    ] = refined.values

    if "curation_tree" not in adata.uns:
        adata.uns["curation_tree"] = {}

    adata.uns["curation_tree"][child] = {
        "parent": cluster,
        "source": "manual_refine",
    }

    adata.uns["curation_tree"][target] = {
        "merged_from": group,
        "source": "manual_merge",
    }

    compute_cluster_summary(
        adata=adata_global,
        cluster_key="curated_cluster",
        stability_key="cluster_stability",
        reassignment_df=curated_reassignments,
        dominant_reassignment_threshold=(
            args.dominant_reassignment_threshold
        ),
        reassignment_display_threshold=(
            args.reassignment_display_threshold
        ),
        merge_candidate_threshold=(
            args.high_reassignment
        ),
    )
