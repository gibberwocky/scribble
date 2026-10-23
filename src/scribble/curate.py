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
    import harmony as hm
    import random
    from datetime import datetime
    import matplotlib.pyplot as plt

    from scribble.import_data import setup_environment
    from scribble.refine import (
        run_clustering,
        restore_counts,
        compute_markers,
    )
    from scribble.utils import (
        serializable_args,
        write_run_manifest,
        compute_cluster_summary,
    )

    PROJECT_DIR = Path(args.project_dir)
    PLOT_DIR = PROJECT_DIR / "scribble/plots"
    MANIFEST_DIR = PROJECT_DIR / "scribble/manifests"
    TABLE_DIR = PROJECT_DIR / "scribble/tables"
    ADATA_DIR = PROJECT_DIR / "scribble/adata"

    setup_environment(sc, np, random, PLOT_DIR)
    input_file = Path(args.input)
    tag = datetime.now().strftime("%Y%m%d_%H%M%S")
    manifest_file = (MANIFEST_DIR / f"{input_file.stem}_curate_{tag}_manifest.json")
    umap_png_file = PLOT_DIR / "UMAP_curate_cluster.png"
    markers_file = (TABLE_DIR / "curate_markers.xlsx")
    curate_report_file = (TABLE_DIR / "curation_report.tsv")
    all_markers_file = (TABLE_DIR / "all_curate_markers.tsv")
    reassignments_file = (TABLE_DIR / "curate_cluster_reassignments.tsv")

    print(f"Loading AnnData → {input_file}")
    adata = sc.read(input_file)

    # Use prior params if available and not over-ridden at CLI
    if (
        not args.ignore_saved_params
        and "scribble_refine_params" in adata.uns
    ):

        saved = adata.uns["scribble_refine_params"]

        cli_overrides = set(
            getattr(
                args,
                "cli_overrides",
                [],
            )
        )

        for param in [
            "batch",
            "hvgs",
            "npcs",
            "neighbors",
            "theta",
            "resolution",
            "auto_resolution",
            "res_min",
            "res_max",
            "res_steps",
            "fine_width",
            "n_repeats",
            "min_cells_per_gene",
            "min_cells_per_cluster",
            "nmarkers",
            "skip_integration",
        ]:

            if param not in saved:
                continue

            if param in cli_overrides:
                continue

            setattr(
                args,
                param,
                saved[param]
            )


    master_reassignments = []
    master_marker_rows = []

    #merge_groups = parse_merge_groups(
    #    args.merge
    #)
    #print("\n[MERGE GROUPS]")
    #print(merge_groups)

    refine_requests = parse_refine_requests(args.refine)
    print("\n[REFINE REQUESTS]")
    print(refine_requests)


    adata.obs["curated_cluster"] = (adata.obs["refine_label"].astype(str))

    if "curation_tree" not in adata.uns:
        adata.uns["curation_tree"] = {}

    for cluster in refine_requests:

        mask = (adata.obs["curated_cluster"].astype(str) == cluster)
        adata_sub = adata[mask].copy()

        print(f"[REFINE] | cluster = {cluster} | cells = {adata_sub.n_obs}")

        adata_sub = restore_counts(adata_sub)
        adata_sub = run_clustering(hm, np, sc, pd, adata_sub, args)

        if adata_sub is None:

            adata.uns["curation_tree"][cluster] = {
                "parent": cluster,
                "source": "refinement_failed",
            }

            continue

        cluster_map = {}

        for i, cl in enumerate(
            sorted(
                adata_sub.obs[
                    "leiden_refined"
                ].astype(str).unique()
            ),
            start=1
        ):
            cluster_map[cl] = (f"{cluster}_S{i}")


        refined = (
            adata_sub.obs["leiden_refined"]
            .astype(str)
            .map(cluster_map)
        )

        adata.obs.loc[
            adata_sub.obs_names,
            "curated_cluster"
        ] = refined.values

        marker_tables = compute_markers(
            sc,
            pd,
            np,
            adata_sub,
            args,
            groupby="leiden_refined",
            min_cluster_size=2,
        )

        for cl, df in marker_tables.items():

            node_label = cluster_map.get(str(cl))

            if node_label is None:
                continue

            tmp = df.copy()

            tmp["curated_cluster"] = node_label
            tmp["parent_cluster"] = cluster
            tmp["level"] = 1
            tmp["marker_type"] = "differential"

            master_marker_rows.append(tmp)


        if "curation_tree" not in adata.uns:
            adata.uns["curation_tree"] = {}

        for child in refined.unique():

            adata.uns["curation_tree"][child] = {
                "parent": cluster,
                "source": "manual_refine",
                "tag": tag,
            }

        if "cluster_stability" in adata_sub.obs:

            adata.obs.loc[adata_sub.obs_names, "cluster_stability"] = (
                adata_sub.obs["cluster_stability"]
            )

        if "cluster_reassignments" in adata_sub.uns:

            tmp = (adata_sub.uns["cluster_reassignments"].copy())

            global_map = {
                str(k): v
                for k, v in cluster_map.items()
            }

            tmp["cluster"] = (
                tmp["cluster"]
                .astype(str)
                .map(global_map)
            )

            tmp["destination"] = (
                tmp["destination"]
                .astype(str)
                .map(global_map)
            )

            tmp["source_cluster"] = cluster

            master_reassignments.append(tmp)




    # --------------------------------------------------
    # Save
    # --------------------------------------------------
    curated_reassignments = None

    if len(master_reassignments) > 0:

        curated_reassignments = pd.concat(
            master_reassignments,
            ignore_index=True,
        )

        curated_reassignments.to_csv(
            reassignments_file,
            sep="\t",
            index=False,
        )

    adata_global = restore_counts(adata)
    adata_global.raw = adata_global.copy()

    adata_global.obs["curated_cluster"] = (
        adata.obs["curated_cluster"]
        .astype(str)
    )

    if "cluster_stability" not in adata_global.obs:
        adata_global.obs["cluster_stability"] = (
            adata.obs["cluster_stability"]
        )

    # --------------------------------------------------
    # Markers
    # --------------------------------------------------
    global_markers = compute_markers(sc, pd, np, adata_global, args, groupby="curated_cluster", min_cluster_size=2)
    with pd.ExcelWriter(markers_file) as writer:

        index_df = (
            adata.obs["curated_cluster"]
            .value_counts()
            .rename_axis("cluster")
            .reset_index(name="n_cells")
        )
        index_df.to_excel(writer, sheet_name="INDEX", index=False)

        for cluster, df in sorted(global_markers.items()):

            df.to_excel(writer, sheet_name=str(cluster)[:31], index=False)



    if len(master_marker_rows) > 0:

        all_markers = pd.concat(
            master_marker_rows,
            ignore_index=True
        )

        all_markers.to_csv(
            all_markers_file,
            sep="\t",
            index=False
        )

    else:
        print(
            "No refined marker tables generated."
        )


    # --------------------------------------------------
    # UMAP
    # --------------------------------------------------
    sc.pl.umap(
        adata_global,
        color="curated_cluster",
        legend_loc="right margin",
        legend_fontsize=7,
        show=False
    )
    plt.savefig(umap_png_file, dpi=300, bbox_inches="tight")
    plt.close()

    # --------------------------------------------------
    # Curate summary
    # --------------------------------------------------
    curated_summary = compute_cluster_summary(
        adata=adata_global,
        cluster_key="curated_cluster",
        stability_key="cluster_stability",
        reassignment_df=curated_reassignments,
        dominant_reassignment_threshold=args.dominant_reassignment_threshold,
        reassignment_display_threshold=args.reassignment_display_threshold,
        merge_candidate_threshold=args.high_reassignment,
    )

    curation_rows = []

    if "curation_tree" in adata.uns:

        for cluster, info in (
            adata.uns["curation_tree"]
            .items()
        ):

            curation_rows.append({
                "cluster": cluster,
                "parent": info.get("parent"),
                "curation_source": info.get("source"),
                "tag": info.get("tag"),
            })

    all_clusters = set(
        adata.obs["curated_cluster"]
        .astype(str)
        .unique()
    )

    reported_clusters = {
        row["cluster"]
        for row in curation_rows
    }

    missing = all_clusters - reported_clusters

    for cluster in sorted(missing):

        curation_rows.append({
            "cluster": cluster,
            "parent": cluster,
            "curation_source": "unmodified",
        })

    curation_df = pd.DataFrame(curation_rows)
    curation_report = (
        curation_df
        .merge(
            curated_summary,
            on="cluster",
            how="left",
        )
    )

    for col in [
        "dominant_reassignment",
        "dominant_reassignment_fraction",
        "reassignment_summary",
        "reassignment_class",
        "merge_candidate",
    ]:
        if col not in curation_report.columns:

            if col == "merge_candidate":
                curation_report[col] = False

            else:
                curation_report[col] = np.nan

    curation_report = curation_report[
        [
            "cluster",
            "n_cells",

            "parent",
            "curation_source",
            "tag",

            "mean_stability",
            "median_stability",
            "sample_entropy",

            "dominant_reassignment",
            "dominant_reassignment_fraction",
            "reassignment_summary",
            "reassignment_class",
            "merge_candidate",
        ]
    ]

    curation_report.to_csv(
        curate_report_file,
        sep="\t",
        index=False,
    )



    if "scribble_curation_history" not in adata.uns:
        adata.uns["scribble_curation_history"] = {}

    adata.uns["scribble_curation_history"][tag] = serializable_args(args)
    adata.uns["scribble_manifest"] = str(manifest_file)


    output = input_file.with_name(f"{input_file.stem}_curated.h5ad")
    print(f"Saving → {output}")
    adata.write(output)

    print(f"Writing run manifest to {manifest_file}")
    outputs = {
        "adata": str(output),
        "UMAP png": str(umap_png_file),
        "markers xlsx": str(markers_file),
        "curation report tsv": str(curate_report_file),
    }

    if len(master_reassignments) > 0:
        outputs["reassignments tsv"] = str(
            reassignments_file
        )

    if len(master_marker_rows) > 0:
        outputs["markers tsv"] = str(
            all_markers_file
        )

    write_run_manifest(
        manifest_file=manifest_file,
        module="curate",
        args=args,
        inputs={
            "input_file": str(input_file),
        },
        outputs=outputs
    )
