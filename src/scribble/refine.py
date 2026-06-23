#!/usr/bin/env python

from pathlib import Path


def run_refine(args):
    import scanpy as sc
    import pandas as pd
    import numpy as np
    import matplotlib.pyplot as plt
    import harmonypy as hm
    import random

    from scribble.import_data import setup_environment
    from scribble.cluster import optimise_resolution
    from scipy.optimize import linear_sum_assignment

    PROJECT_DIR = Path(args.project_dir)
    PLOT_DIR = PROJECT_DIR / "scribble/plots"
    TABLE_DIR = PROJECT_DIR / "scribble/tables"

    setup_environment(sc, np, random, PLOT_DIR)

    input_file = Path(args.input)
    decision_file = Path(args.decisions)

    print(f"Loading AnnData → {input_file}")
    adata = sc.read(input_file)

    print(f"Loading decisions → {decision_file}")
    decisions = pd.read_csv(decision_file, sep="\t")

    # --------------------------------------------------
    # Initialise hierarchical labels
    # --------------------------------------------------
    print("Initialising L2 labels...")
    adata.obs["leiden_L2"] = adata.obs["leiden"].astype(str) + "-0"

    # --------------------------------------------------
    # Identify merge groups
    # --------------------------------------------------
    groups = [
        g for g in decisions["merge_group"].dropna().unique()
        if g != ""
    ]

    if len(groups) == 0:
        print("No merge groups found. Nothing to refine.")
        return

    print(f"Found {len(groups)} merge group(s)")

    # --------------------------------------------------
    # Process each group
    # --------------------------------------------------
    for group in groups:

        print(f"\nProcessing {group}")

        subset_clusters = decisions.loc[
            decisions["merge_group"] == group, "cluster"
        ].astype(str).tolist()

        print(f"Clusters: {subset_clusters}")

        adata_sub = adata[adata.obs["leiden"].isin(subset_clusters)].copy()

        print(f"Cells: {adata_sub.n_obs}")

        if adata_sub.n_obs < args.min_cells_per_group:
            print("Skipping (too few cells)")
            continue

        # --------------------------------------------------
        # Preprocessing
        # --------------------------------------------------
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

        # --------------------------------------------------
        # Harmony
        # --------------------------------------------------
        print(f"Running Harmony (theta={args.theta})")

        ho = hm.run_harmony(
            adata_sub.obsm["X_pca"],
            adata_sub.obs,
            args.batch,
            theta=args.theta
        )

        adata_sub.obsm["X_pca_harmony"] = ho.Z_corr

        # --------------------------------------------------
        # Neighbors
        # --------------------------------------------------
        sc.pp.neighbors(
            adata_sub,
            use_rep="X_pca_harmony",
            n_pcs=args.npcs,
            n_neighbors=args.neighbors
        )

        # --------------------------------------------------
        # Resolution optimisation
        # --------------------------------------------------
        if args.auto_resolution:

            best_res, _, _ = optimise_resolution(
                np, pd, sc,
                adata_sub,
                "X_pca_harmony",
                args.neighbors,
                (args.res_min, args.res_max),
                args.fine_width,
                args.res_steps
            )

            res = best_res
        else:
            res = args.resolution

        print(f"Using resolution: {res}")

        # --------------------------------------------------
        # Clustering
        # --------------------------------------------------
        sc.tl.leiden(
            adata_sub,
            resolution=res,
            key_added="leiden_refined",
            flavor="igraph",
            directed=False,
            n_iterations=2
        )

        # --------------------------------------------------
        # Stability
        # --------------------------------------------------
        if args.n_repeats > 1:

            print(f"Computing stability ({args.n_repeats} repeats)")

            ref = adata_sub.obs["leiden_refined"].astype(str).values
            all_assignments = []

            def align(reference, current):
                contingency = pd.crosstab(current, reference)
                cost = -contingency.values
                r, c = linear_sum_assignment(cost)

                mapping = {
                    contingency.index[r[i]]: contingency.columns[c[i]]
                    for i in range(len(r))
                }

                return (
                    pd.Series(current)
                    .map(mapping)
                    .fillna(current)
                    .values
                )

            for i in range(args.n_repeats):
                sc.tl.leiden(
                    adata_sub,
                    resolution=res,
                    key_added=f"tmp_{i}",
                    random_state=i
                )

                raw = adata_sub.obs[f"tmp_{i}"].astype(str).values
                aligned = align(ref, raw)
                all_assignments.append(aligned)

            all_assignments = np.array(all_assignments)

            stability = np.array([
                np.unique(all_assignments[:, i], return_counts=True)[1].max()
                / args.n_repeats
                for i in range(all_assignments.shape[1])
            ])

            adata_sub.obs["cluster_stability"] = stability

            for i in range(args.n_repeats):
                del adata_sub.obs[f"tmp_{i}"]

        # --------------------------------------------------
        # Within-lineage markers
        # --------------------------------------------------
        print(f"Computing within-lineage markers for {group}")

        sc.tl.rank_genes_groups(
            adata_sub,
            "leiden_refined",
            method="wilcoxon",
            use_raw=True
        )

        result = adata_sub.uns["rank_genes_groups"]
        marker_clusters = result["names"].dtype.names

        markers_file = TABLE_DIR / f"{input_file.stem}_{group}_within_lineage_markers.xlsx"

        with pd.ExcelWriter(markers_file, engine="openpyxl") as writer:

            X = adata_sub.raw.X if adata_sub.raw is not None else adata_sub.X
            var_names = adata_sub.raw.var_names if adata_sub.raw is not None else adata_sub.var_names

            parent_series = adata.obs.loc[adata_sub.obs_names, "leiden"].astype(str)
            parent_label = "+".join(sorted(parent_series.unique(), key=int))

            for cl in marker_clusters:

                genes = result["names"][cl]

                valid = [g for g in genes if g in var_names]
                gene_idx = [var_names.get_loc(g) for g in valid]

                expr = X[:, gene_idx]
                expr = expr.A if hasattr(expr, "A") else expr

                cluster_cells = adata_sub.obs["leiden_refined"] == cl
                other_cells = ~cluster_cells

                df = pd.DataFrame({
                    "gene": valid,
                    "logfoldchange": result["logfoldchanges"][cl][:len(valid)],
                    "score": result["scores"][cl][:len(valid)],
                    "pvals": result["pvals"][cl][:len(valid)],
                    "pvals_adj": result["pvals_adj"][cl][:len(valid)],
                })

                pct_in = (expr[cluster_cells.values] > 0).mean(axis=0)
                pct_out = (expr[other_cells.values] > 0).mean(axis=0)

                df["pct_in"] = pct_in
                df["pct_out"] = pct_out
                df["pct_diff"] = pct_in - pct_out

                df["mean_in"] = expr[cluster_cells.values].mean(axis=0)
                df["mean_out"] = expr[other_cells.values].mean(axis=0)

                df["cluster_size"] = cluster_cells.sum()

                df = df.sort_values(
                    ["pvals_adj", "pct_diff", "logfoldchange"],
                    ascending=[True, False, False]
                ).head(args.nmarkers)

                sheet_name = f"cluster_{parent_label}-{cl}"[:31]
                df.to_excel(writer, sheet_name=sheet_name, index=False)

        # --------------------------------------------------
        # Map refined labels back
        # --------------------------------------------------
        refined_labels = (
            adata.obs.loc[adata_sub.obs_names, "leiden"].astype(str)
            + "-"
            + adata_sub.obs["leiden_refined"].astype(str)
        )

        adata.obs.loc[adata_sub.obs_names, "leiden_L2"] = refined_labels

    # --------------------------------------------------
    # Global markers (L2)
    # --------------------------------------------------
    print("\nComputing global L2 markers...")

    sc.tl.rank_genes_groups(
        adata,
        "leiden_L2",
        method="wilcoxon",
        use_raw=True
    )

    result = adata.uns["rank_genes_groups"]
    marker_clusters = result["names"].dtype.names

    markers_file = TABLE_DIR / f"{input_file.stem}_L2_markers.xlsx"

    with pd.ExcelWriter(markers_file, engine="openpyxl") as writer:

        X = adata.raw.X if adata.raw is not None else adata.X
        var_names = adata.raw.var_names if adata.raw is not None else adata.var_names

        for cl in marker_clusters:

            genes = result["names"][cl]

            valid = [g for g in genes if g in var_names]
            gene_idx = [var_names.get_loc(g) for g in valid]

            expr = X[:, gene_idx]
            expr = expr.A if hasattr(expr, "A") else expr

            cluster_cells = adata.obs["leiden_L2"] == cl
            other_cells = ~cluster_cells

            df = pd.DataFrame({
                "gene": valid,
                "logfoldchange": result["logfoldchanges"][cl][:len(valid)],
                "score": result["scores"][cl][:len(valid)],
                "pvals": result["pvals"][cl][:len(valid)],
                "pvals_adj": result["pvals_adj"][cl][:len(valid)],
            })

            pct_in = (expr[cluster_cells.values] > 0).mean(axis=0)
            pct_out = (expr[other_cells.values] > 0).mean(axis=0)

            df["pct_in"] = pct_in
            df["pct_out"] = pct_out
            df["pct_diff"] = pct_in - pct_out

            df["mean_in"] = expr[cluster_cells.values].mean(axis=0)
            df["mean_out"] = expr[other_cells.values].mean(axis=0)

            df["cluster_size"] = cluster_cells.sum()

            df = df.sort_values(
                ["pvals_adj", "pct_diff", "logfoldchange"],
                ascending=[True, False, False]
            ).head(args.nmarkers)

            sheet_name = f"cluster_{cl}"[:31]
            df.to_excel(writer, sheet_name=sheet_name, index=False)

    # --------------------------------------------------
    # Finalise
    # --------------------------------------------------
    adata.obs["leiden_L2_id"] = (
        adata.obs["leiden_L2"].astype("category").cat.codes
    )

    output_file = input_file.with_name(f"{input_file.stem}_refined.h5ad")

    print(f"\nSaving → {output_file}")
    adata.write(output_file)
