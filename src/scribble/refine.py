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
    import random
    import harmonypy as hm
    from pathlib import Path

    from scribble.import_data import setup_environment
    from scribble.cluster import optimise_resolution

    PROJECT_DIR = Path(args.project_dir)
    PLOT_DIR = PROJECT_DIR / "scribble/plots"
    TABLE_DIR = PROJECT_DIR / "scribble/tables"

    setup_environment(sc, np, random, PLOT_DIR)

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

    def _robust_hvg(adata_sub):
        """
        Adaptive HVG selection with fallback for Seurat v3 loess instability.
        Reduces number of HVGs only if numerical failure occurs.
        """

        current_n = args.hvgs
        min_genes = max(1000, int(args.hvgs * 0.4))  # don't go too low
        step = max(250, args.hvgs // 10)             # adaptive step size

        attempt = 1

        while current_n >= min_genes:
            try:
                sc.pp.highly_variable_genes(
                    adata_sub,
                    n_top_genes=current_n,
                    batch_key=args.batch,
                    flavor="seurat_v3"
                )

                if current_n != args.hvgs:
                    print(f"[HVG] recovered: using {current_n} genes (after fallback)")

                return current_n  # ✅ success

            except Exception as e:
                if "reciprocal condition number" in str(e):
                    print(f"[HVG WARNING] instability at {current_n} genes → retrying")

                    current_n -= step
                    attempt += 1

                    if current_n < min_genes:
                        break
                else:
                    raise e  # unrelated error

        # --------------------------------------------------
        # Last resort: fallback to 'seurat' (loess-free)
        # --------------------------------------------------
        print(f"[HVG WARNING] fallback failed → switching to flavor='seurat'")

        sc.pp.highly_variable_genes(
            adata_sub,
            n_top_genes=args.hvgs,
            batch_key=None,
            flavor="seurat"
        )

        return args.hvgs


    def _run_clustering(adata_sub):

        sc.pp.filter_genes(adata_sub, min_cells=args.min_cells_per_gene)

        _robust_hvg(adata_sub)

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
            flavor="igraph", directed=False, n_iterations=2,
            random_state=0
        )

        return adata_sub


    def _compute_markers(adata_de, groupby):

        import numpy as np
        import pandas as pd

        adata_de.obs[groupby] = adata_de.obs[groupby].astype("category")

        adata_de.obs[groupby] = adata_de.obs[groupby].cat.set_categories(
            sorted(adata_de.obs[groupby].unique()),
            ordered=True
        )

        sc.tl.rank_genes_groups(
            adata_de,
            groupby,
            method="wilcoxon",
            use_raw=True
        )

        result = adata_de.uns["rank_genes_groups"]

        rg_groups = result["names"].dtype.names
        categories = adata_de.obs[groupby].astype("category").cat.categories

        X = adata_de.raw.X
        var_names = adata_de.raw.var_names
        var_index = pd.Series(range(len(var_names)), index=var_names)

        marker_tables = {}

        for i, rg_key in enumerate(rg_groups):

            label = categories[i] if i < len(categories) else rg_key

            genes = np.array(result["names"][rg_key]).ravel()
            logfc = np.array(result["logfoldchanges"][rg_key]).ravel()
            scores = np.array(result["scores"][rg_key]).ravel()
            pvals = np.array(result["pvals"][rg_key]).ravel()
            pvals_adj = np.array(result["pvals_adj"][rg_key]).ravel()

            # filter to genes in var
            valid_mask = np.isin(genes, var_names)

            genes = genes[valid_mask]
            logfc = logfc[valid_mask]
            scores = scores[valid_mask]
            pvals = pvals[valid_mask]
            pvals_adj = pvals_adj[valid_mask]

            gene_idx = var_index.loc[genes].values

            # define groups
            cluster_cells = adata_de.obs[groupby] == label
            other_cells = ~cluster_cells

            expr = X[:, gene_idx].tocsr()

            pct_in = (expr[cluster_cells.values] > 0).mean(axis=0)
            pct_out = (expr[other_cells.values] > 0).mean(axis=0)

            pct_in = np.asarray(pct_in).ravel()
            pct_out = np.asarray(pct_out).ravel()

            mean_in = expr[cluster_cells.values].mean(axis=0)
            mean_out = expr[other_cells.values].mean(axis=0)

            mean_in = np.asarray(mean_in).ravel()
            mean_out = np.asarray(mean_out).ravel()

            df = pd.DataFrame({
                "gene": genes,
                "logfoldchange": logfc,
                "score": scores,
                "pvals": pvals,
                "pvals_adj": pvals_adj,
                "pct_in": pct_in,
                "pct_out": pct_out,
                "pct_diff": pct_in - pct_out,
                "mean_in": mean_in,
                "mean_out": mean_out,
                "cluster_size": cluster_cells.sum()
            })

            df = df.sort_values(
                ["pvals_adj", "pct_diff", "logfoldchange"],
                ascending=[True, False, False]
            ).head(args.nmarkers)

            marker_tables[str(label)] = df

        return marker_tables

    def _is_marker_strong(marker_df, threshold):
        """
        Quick proxy: are top markers clearly separable?
        Uses logfoldchange as proxy for biological signal.
        """
        if marker_df is None or marker_df.empty:
            return False

        return marker_df["logfoldchange"].head(5).mean() >= threshold


    def _refine_task(task):

        clusters = task["clusters"]
        level = task["level"]

        mask = adata.obs["leiden"].isin(clusters)
        adata_sub = adata[mask].copy()

        print(f"[Refine] Level {level} | clusters = {clusters} | cells = {adata_sub.n_obs}")

        if adata_sub.n_obs < args.min_cells_per_group:
            return []

        adata_sub = restore_counts(adata_sub)

        adata_sub = _run_clustering(adata_sub)

        # -----------------------
        # Map refined labels
        # -----------------------
        parent_labels = adata.obs.loc[adata_sub.obs_names, "leiden"].astype(str)

        refined = (
            parent_labels
            + "-"
            + adata_sub.obs["leiden_refined"].astype(str)
        )

        # keep as pandas Series (CRITICAL)
        refined = refined.astype(str)
        refined = "L" + str(level) + "_" + refined

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

            clusters_str = "+".join(task["clusters"])
            out_file = TABLE_DIR / f"L{task['level']}_{clusters_str}_markers.xlsx"
            with pd.ExcelWriter(out_file) as writer:
                for cl, df in marker_tables.items():
                    df.to_excel(writer, sheet_name=str(cl), index=False)

        # -----------------------
        # Recursive refinement
        # -----------------------
        new_tasks = []

        # --------------------------------------------------
        # Adaptive recursive refinement
        # --------------------------------------------------
        new_tasks = []

        if level < args.max_refine_depth:

            counts = adata_sub.obs["leiden_refined"].value_counts()

            # compute mean stability per cluster (if available)
            if "cluster_stability" in adata_sub.obs:
                stability = (
                    adata_sub
                    .obs
                    .groupby("leiden_refined")["cluster_stability"]
                    .mean()
                )
            else:
                stability = {}

            for cl in counts.index:

                cluster_size = counts[cl]
                cluster_label = str(cl)

                marker_df = marker_tables.get(cluster_label, None)

                # ------------------------------
                # DEBUG / DECISION LOGGING
                # ------------------------------
                print(f"  Cluster {cluster_label}: size={cluster_size}", end="")

                if cluster_label in stability:
                    print(f", stability={stability[cluster_label]:.3f}", end="")

                if marker_df is not None:
                    top_fc = marker_df["logfoldchange"].head(5).mean()
                    print(f", marker_strength={top_fc:.2f}", end="")

                print()

                # ------------------------------
                # RULE 1: size filter
                # ------------------------------
                if cluster_size < args.min_cells_per_cluster:
                    print("    → STOP (too small)")
                    continue

                # ------------------------------
                # RULE 2: stability filter
                # ------------------------------
                if cluster_label in stability:
                    if stability[cluster_label] >= args.stability_threshold:
                        print("    → STOP (stable)")
                        continue

                # ------------------------------
                # RULE 3: marker strength filter
                # ------------------------------
                if marker_df is not None:
                    if _is_marker_strong(marker_df, args.marker_strength_threshold):
                        print("    → STOP (already well separated)")
                        continue

                # ------------------------------
                # OTHERWISE → REFINE
                # ------------------------------
                print("    → REFINE further")

                new_tasks.append({
                    "clusters": [f"{clusters[0]}-{cluster_label}"],
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
