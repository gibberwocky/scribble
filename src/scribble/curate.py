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
    curation_lineage_file = (TABLE_DIR / "curation_lineage.tsv")
    all_markers_file = (TABLE_DIR / "all_curate_markers.tsv")
    reassignments_file = (TABLE_DIR / "curate_cluster_reassignments.tsv")

    print(f"Loading AnnData → {input_file}")
    adata = sc.read(input_file)

    # Use prior params if available and not over-ridden at CLI
    if (
        not args.ignore_saved_params and
        "scribble_refine_params" in adata.uns
    ):

        saved = adata.uns["scribble_refine_params"]

        print("\nUsing clustering parameters stored in AnnData:")
        print(saved.keys())

        cli_overrides = set(
            getattr(
                args,
                "cli_overrides",
                [],
            )
        )

        for param, value in saved.items():

            # Explicit CLI value wins
            if param in {"manifest", "command"}:
                continue

            if param in cli_overrides:
                continue

            if hasattr(args, param):
                setattr(args, param, value)
                print(f" {param} = {value}")


    lineage_df = None

    if "scribble_lineage" in adata.uns:

        lineage_df = pd.DataFrame.from_dict(
            adata.uns["scribble_lineage"],
            orient="index"
        ).reset_index(drop=True)

        print(
            f"Loaded lineage records "
            f"({len(lineage_df)} rows)"
        )

    master_reassignments = []
    master_marker_rows = []

    refine_requests = parse_refine_requests(args.refine)
    print("\n[REFINE REQUESTS]")
    print(refine_requests)


    adata.obs["curated_cluster"] = (adata.obs["refine_label"].astype(str))

    for cluster in refine_requests:

        mask = (adata.obs["curated_cluster"].astype(str) == cluster)
        adata_sub = adata[mask].copy()

        print(f"[REFINE] | cluster = {cluster} | cells = {adata_sub.n_obs}")

        adata_sub = restore_counts(adata_sub)
        adata_sub = run_clustering(hm, np, sc, pd, adata_sub, args)

        if adata_sub is None:

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


        for child in refined.unique():

            adata.uns["scribble_lineage"][child] = {
                "cluster": child,
                "parent": cluster,
                "level": (
                    lineage_df.loc[
                        lineage_df["cluster"]
                        == cluster,
                        "level"
                    ].iloc[0] + 1
                    if lineage_df is not None
                    and cluster in lineage_df["cluster"].values
                    else 1
                ),
                "lineage_path": (
                    lineage_df.loc[
                        lineage_df["cluster"]
                        == cluster,
                        "lineage_path"
                    ].iloc[0]
                    + f" -> {child}"
                    if lineage_df is not None
                    and cluster in lineage_df["cluster"].values
                    else f"{cluster} -> {child}"
                ),
                "refinement_source": "manual_refine",
                "parent_clusters": cluster,
                "n_parent_clusters": 1,
                "source": "curate",
                "curation_tag": tag,
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
    if "X_umap" in adata.obsm:
        adata_global.obsm["X_umap"] = adata.obsm["X_umap"].copy()
    sc.pp.normalize_total(adata_global)
    sc.pp.log1p(adata_global)
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

    lineage_df = pd.DataFrame.from_dict(
        adata.uns["scribble_lineage"],
        orient="index"
    ).reset_index(
        drop=True
    )

    all_clusters = set(
        adata.obs["curated_cluster"]
        .astype(str)
        .unique()
    )

    recorded_clusters = set(
        lineage_df["cluster"]
    )

    missing = all_clusters - recorded_clusters

    rows = []

    for cluster in sorted(missing):

        rows.append({
            "cluster": cluster,
            "parent": None,
            "level": 0,
            "lineage_path": cluster,
            "refinement_source": "unmodified",
            "parent_clusters": cluster,
            "n_parent_clusters": 1,
            "source": "unmodified",
            "curation_tag": tag,
        })

    if rows:

        lineage_df = pd.concat(
            [
                lineage_df,
                pd.DataFrame(rows),
            ],
            ignore_index=True,
        )

    lineage_df.to_csv(
        curation_lineage_file,
        sep="\t",
        index=False,
    )

    curation_report = (
        lineage_df
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
            "level",
            "lineage_path",

            "refinement_source",
            "parent_clusters",
            "n_parent_clusters",

            "curation_source",
            "curation_tag",

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

    print(
        f"Curated clusters: "
        f"{curation_report['cluster'].nunique()}"
    )

    curation_report.to_csv(
        curate_report_file,
        sep="\t",
        index=False,
    )

    if "curation_reports" not in adata.uns:
        adata.uns["curation_reports"] = {}

    adata.uns["curation_reports"][tag] = (
        str(curate_report_file)
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
        "curation lineage tsv": str(curation_lineage_file),
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
