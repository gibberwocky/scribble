#!/usr/bin/env python

from pathlib import Path

def restore_counts(adata):
    import scanpy as sc
    import pandas as pd

    if "counts_full" not in adata.uns:
        raise ValueError("Full counts not found in adata.uns")

    if "obs_full_names" not in adata.uns:
        raise ValueError("Missing obs_full_names for restoring counts")

    # Full data
    counts_full = adata.uns["counts_full"]
    var_full = adata.uns["var_full"]

    # Map cells to indices
    full_obs_names = pd.Index(adata.uns["obs_full_names"])
    cell_idx = full_obs_names.get_indexer(adata.obs_names)

    if (cell_idx < 0).any():
        raise ValueError("Some cells not found in full dataset")

    # Subset counts to correct cells
    X = counts_full[cell_idx, :]

    return sc.AnnData(
        X=X,
        obs=adata.obs.copy(),
        var=var_full.copy()
    )


def run_refine(args):
    import scanpy as sc
    import pandas as pd
    import numpy as np
    import harmonypy as hm
    from pathlib import Path

    from scribble.import_data import setup_environment
    from scribble.cluster import optimise_resolution

    PROJECT_DIR = Path(args.project_dir)
    PLOT_DIR = PROJECT_DIR / "scribble/plots"
    TABLE_DIR = PROJECT_DIR / "scribble/tables"

    setup_environment(sc, np, None, PLOT_DIR)

    input_file = Path(args.input)
    decisions = pd.read_csv(args.decisions, sep="\t")

    print(f"Loading AnnData → {input_file}")
    adata = sc.read(input_file)

    # --------------------------------------------------
    # Init labels + lineage
    # --------------------------------------------------
    adata.obs["leiden_L2"] = adata.obs["leiden"].astype(str) + "-0"

    if "lineage_tree" not in adata.uns:
        adata.uns["lineage_tree"] = {}

    # --------------------------------------------------
    # Build initial tasks
    # --------------------------------------------------
    tasks = []

    subset_clusters = decisions.loc[
        decisions["action"] == "subset", "cluster"
    ].astype(str)

    for cl in subset_clusters:
        tasks.append({
            "clusters": [cl],
            "level": 2,
        })

    if len(tasks) == 0:
        print("No refinement tasks.")
        return

    # --------------------------------------------------
    # INTERNAL FUNCTIONS
    # --------------------------------------------------

    def _run_clustering(adata_sub):

        sc.pp.filter_genes(adata_sub, min_cells=args.min_cells_per_gene)

        sc.pp.highly_variable_genes(
            adata_sub,
            n_top_genes=args.hvgs,
            batch_key=args.batch,
            flavor="seurat_v3"
        )

        sc.pp.normalize_total(adata_sub)
        sc.pp.log1p(adata_sub)

        adata_sub.raw = adata_sub.copy()
        adata_sub = adata_sub[:, adata_sub.var.highly_variable].copy()

        if not args.no_scale:
            sc.pp.scale(adata_sub, max_value=10)

        sc.tl.pca(adata_sub, n_comps=args.npcs)

        ho = hm.run_harmony(
            adata_sub.obsm["X_pca"],
            adata_sub.obs,
            args.batch,
            theta=args.theta
        )

        adata_sub.obsm["X_pca_harmony"] = ho.Z_corr

        sc.pp.neighbors(
            adata_sub,
            use_rep="X_pca_harmony",
            n_neighbors=args.neighbors
        )

        if args.auto_resolution:
            res, _, _ = optimise_resolution(
                np, pd, sc, adata_sub,
                "X_pca_harmony",
                args.neighbors,
                (args.res_min, args.res_max),
                args.fine_width,
                args.res_steps
            )
        else:
            res = args.resolution

        sc.pp.neighbors(
            adata_sub,
            use_rep="X_pca_harmony",
            n_neighbors=args.neighbors
        )

        sc.tl.leiden(
            adata_sub,
            resolution=res,
            key_added="leiden_refined",
            random_state=0
        )

        return adata_sub


    def _compute_markers(adata_de, groupby):

        sc.tl.rank_genes_groups(
            adata_de,
            groupby,
            method="wilcoxon",
            use_raw=True
        )

        result = adata_de.uns["rank_genes_groups"]

        categories = adata_de.obs[groupby].astype("category").cat.categories

        marker_tables = {}

        for i, label in enumerate(categories):

            genes = np.array(result["names"][str(i)])
            logfc = np.array(result["logfoldchanges"][str(i)])

            df = pd.DataFrame({
                "gene": genes,
                "logFC": logfc
            }).head(args.nmarkers)

            marker_tables[label] = df

        return marker_tables


    def _refine_task(task):

        clusters = task["clusters"]
        level = task["level"]

        mask = adata.obs["leiden"].isin(clusters)
        adata_sub = adata[mask].copy()

        if adata_sub.n_obs < args.min_cells_per_group:
            return []

        adata_sub = restore_counts(adata_sub)

        adata_sub = _run_clustering(adata_sub)

        # -----------------------
        # Map refined labels
        # -----------------------
        refined = (
            adata.obs.loc[adata_sub.obs_names, "leiden"]
            + "-" +
            adata_sub.obs["leiden_refined"].astype(str)
        )

        adata.obs.loc[adata_sub.obs_names, "leiden_L2"] = refined

        # -----------------------
        # Lineage tracking
        # -----------------------
        for parent in clusters:
            children = refined.unique().tolist()

            adata.uns["lineage_tree"].setdefault(parent, [])
            adata.uns["lineage_tree"][parent].extend(children)

        # -----------------------
        # Markers
        # -----------------------
        cluster_sizes = adata_sub.obs["leiden_refined"].value_counts()
        valid = cluster_sizes[cluster_sizes >= 2].index

        if len(valid) >= 2:

            adata_de = adata_sub[
                adata_sub.obs["leiden_refined"].isin(valid)
            ].copy()

            marker_tables = _compute_markers(
                adata_de,
                "leiden_refined"
            )

            out_file = TABLE_DIR / f"{task['clusters']}_markers.xlsx"

            with pd.ExcelWriter(out_file) as writer:
                for cl, df in marker_tables.items():
                    df.to_excel(writer, sheet_name=str(cl), index=False)

        # -----------------------
        # Recursive refinement
        # -----------------------
        new_tasks = []

        if level < args.max_refine_depth:

            counts = adata_sub.obs["leiden_refined"].value_counts()

            for cl in counts.index:
                if counts[cl] >= args.min_cells_per_group:
                    new_tasks.append({
                        "clusters": [f"{clusters[0]}-{cl}"],
                        "level": level + 1
                    })

        return new_tasks


    # --------------------------------------------------
    # MAIN LOOP (safe recursion)
    # --------------------------------------------------

    i = 0
    while i < len(tasks):

        new_tasks = _refine_task(tasks[i])
        tasks.extend(new_tasks)

        i += 1

    # --------------------------------------------------
    # Final global markers
    # --------------------------------------------------
    if len(adata.obs["leiden_L2"].unique()) > 1:

        markers = _compute_markers(adata, "leiden_L2")

        with pd.ExcelWriter(TABLE_DIR / "L2_markers.xlsx") as writer:
            for cl, df in markers.items():
                df.to_excel(writer, sheet_name=str(cl), index=False)

    # --------------------------------------------------
    # Save
    # --------------------------------------------------
    output = input_file.with_name(f"{input_file.stem}_refined.h5ad")
    print(f"Saving → {output}")
    adata.write(output)
