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

    if "counts" not in adata.layers:
        raise ValueError("Input AnnData must contain layers['counts']")

    if adata.raw is None:
        raise ValueError("Input AnnData must have .raw (log-normalised full matrix)")

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

        # --------------------------------------------------
        # Subset once (clean)
        # --------------------------------------------------
        mask = adata.obs["leiden"].isin(subset_clusters)

        adata_sub = adata[mask].copy()

        print(f"Cells: {adata_sub.n_obs}")

        # Skip small groups
        if adata_sub.n_obs < args.min_cells_per_group:
            print("Skipping (too few cells)")
            continue

        # --------------------------------------------------
        # Preprocessing
        # --------------------------------------------------

        # Restore full gene space (raw counts)
        # NOTE: this replaces log-normalised data with raw counts
        adata_sub = restore_counts(adata_sub)
        # Store counts
        adata_sub.layers["counts"] = adata_sub.X.copy()

        # Gene filtering
        sc.pp.filter_genes(adata_sub, min_cells=args.min_cells_per_gene)

        # HVG identification
        sc.pp.highly_variable_genes(
            adata_sub,
            n_top_genes=args.hvgs,
            batch_key=args.batch,
            flavor="seurat_v3"
        )

        # Normalisation
        sc.pp.normalize_total(adata_sub)
        sc.pp.log1p(adata_sub)

        # Store raw
        adata_sub.raw = adata_sub.copy()

        # Subset on HVGs
        adata_sub = adata_sub[:, adata_sub.var.highly_variable].copy()

        # Optional scaling
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

                cur_series = pd.Series(current)

                aligned = cur_series.map(mapping).fillna(cur_series).values

                return aligned


            for i in range(args.n_repeats):
                sc.tl.leiden(
                    adata_sub,
                    resolution=res,
                    key_added=f"tmp_{i}",
                    random_state=i,
                    flavor="igraph",
                    directed=False,
                    n_iterations=2
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

        # ---- Filter clusters for DE ----
        cluster_sizes = adata_sub.obs["leiden_refined"].value_counts()
        valid_clusters = cluster_sizes[cluster_sizes >= 2].index

        adata_de = None
        result = None

        if len(valid_clusters) < 2:
            print("Skipping within-lineage DE: fewer than 2 valid clusters")
        else:
            adata_de = adata_sub[
                adata_sub.obs["leiden_refined"].isin(valid_clusters)
            ].copy()

            sc.tl.rank_genes_groups(
                adata_de,
                "leiden_refined",
                method="wilcoxon",
                use_raw=True
            )

            result = adata_de.uns["rank_genes_groups"]

        # ---- Safe guard ----
        if result is None:
            print(f"No within-lineage DE results for {group}, skipping export")
            return

        marker_clusters = result["names"].dtype.names

        if marker_clusters is None or len(marker_clusters) == 0:
            print(f"No marker clusters found for {group}, skipping export")
        else:
            markers_file = TABLE_DIR / f"{input_file.stem}_{group}_within_lineage_markers.xlsx"

            with pd.ExcelWriter(markers_file, engine="openpyxl") as writer:

                if adata_de.raw is None:
                    raise ValueError("Expected .raw for marker extraction")

                X = adata_de.raw.X
                var_names = adata_de.raw.var_names

                parent_series = adata.obs.loc[adata_sub.obs_names, "leiden"].astype(str)
                parent_label = "+".join(sorted(parent_series.unique(), key=int))

                sheets_written = 0

                for cl in marker_clusters:

                    genes = np.array(result["names"][cl])
                    logfc_all = np.array(result["logfoldchanges"][cl])
                    scores_all = np.array(result["scores"][cl])
                    pvals_all = np.array(result["pvals"][cl])
                    pvals_adj_all = np.array(result["pvals_adj"][cl])

                    valid_mask = np.isin(genes, var_names)

                    valid = genes[valid_mask]
                    if valid.size == 0:
                        continue

                    logfc = logfc_all[valid_mask]
                    scores = scores_all[valid_mask]
                    pvals = pvals_all[valid_mask]
                    pvals_adj = pvals_adj_all[valid_mask]

                    gene_idx = [var_names.get_loc(g) for g in valid]

                    expr = X[:, gene_idx]
                    expr = expr.toarray() if hasattr(expr, "toarray") else np.asarray(expr)

                    if expr.shape[1] != len(valid):
                        continue

                    cluster_cells = adata_de.obs["leiden_refined"] == cl
                    other_cells = ~cluster_cells

                    df = pd.DataFrame({
                        "gene": valid,
                        "logfoldchange": logfc,
                        "score": scores,
                        "pvals": pvals,
                        "pvals_adj": pvals_adj,
                    }).reset_index(drop=True)

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

                    if df.empty:
                        continue

                    sheet_name = f"cluster_{parent_label}-{cl}"[:31]
                    df.to_excel(writer, sheet_name=sheet_name, index=False)
                    sheets_written += 1

                if sheets_written == 0:
                    print(f"No valid marker sheets written for {group}, skipping file")

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

        cluster_sizes = adata.obs["leiden_L2"].value_counts()
        valid_clusters = cluster_sizes[cluster_sizes >= 2].index

        adata_de = None
        result = None

        if len(valid_clusters) < 2:
            print("Skipping global L2 DE: fewer than 2 valid clusters")
        else:
            adata_de = adata[
                adata.obs["leiden_L2"].isin(valid_clusters)
            ].copy()

            sc.tl.rank_genes_groups(
                adata_de,
                "leiden_L2",
                method="wilcoxon",
                use_raw=True
            )

            result = adata_de.uns["rank_genes_groups"]

        # ---- Safe guard ----
        if result is None:
            print("No global L2 DE results, skipping export")
            return

        marker_clusters = result["names"].dtype.names

        if marker_clusters is None or len(marker_clusters) == 0:
            print("No global L2 marker clusters found, skipping export")
        else:
            markers_file = TABLE_DIR / f"{input_file.stem}_L2_markers.xlsx"

            with pd.ExcelWriter(markers_file, engine="openpyxl") as writer:

                if adata_de.raw is None:
                    raise ValueError("Expected .raw for marker extraction")

                X = adata_de.raw.X
                var_names = adata_de.raw.var_names

                sheets_written = 0

                for cl in marker_clusters:

                    genes = np.array(result["names"][cl])
                    logfc_all = np.array(result["logfoldchanges"][cl])
                    scores_all = np.array(result["scores"][cl])
                    pvals_all = np.array(result["pvals"][cl])
                    pvals_adj_all = np.array(result["pvals_adj"][cl])

                    valid_mask = np.isin(genes, var_names)

                    valid = genes[valid_mask]
                    if valid.size == 0:
                        continue

                    logfc = logfc_all[valid_mask]
                    scores = scores_all[valid_mask]
                    pvals = pvals_all[valid_mask]
                    pvals_adj = pvals_adj_all[valid_mask]

                    gene_idx = [var_names.get_loc(g) for g in valid]

                    expr = X[:, gene_idx]
                    expr = expr.toarray() if hasattr(expr, "toarray") else np.asarray(expr)

                    if expr.shape[1] != len(valid):
                        continue

                    cluster_cells = adata_de.obs["leiden_L2"] == cl
                    other_cells = ~cluster_cells

                    df = pd.DataFrame({
                        "gene": valid,
                        "logfoldchange": logfc,
                        "score": scores,
                        "pvals": pvals,
                        "pvals_adj": pvals_adj,
                    }).reset_index(drop=True)

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

                    if df.empty:
                        continue

                    sheet_name = f"cluster_{cl}"[:31]
                    df.to_excel(writer, sheet_name=sheet_name, index=False)
                    sheets_written += 1

                if sheets_written == 0:
                    print("No valid marker sheets written, skipping file")

    # --------------------------------------------------
    # Finalise
    # --------------------------------------------------
    adata.obs["leiden_L2_id"] = (
        adata.obs["leiden_L2"].astype("category").cat.codes
    )

    output_file = input_file.with_name(f"{input_file.stem}_refined.h5ad")

    print(f"\nSaving → {output_file}")
    adata.write(output_file)
