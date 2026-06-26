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
    # Build refinement tasks (merge + subset)
    # --------------------------------------------------

    tasks = []

    # --- MERGE GROUPS ---
    merge_groups = [
        g for g in decisions["merge_group"].dropna().unique()
        if g != ""
    ]

    for group in merge_groups:
        clusters = decisions.loc[
            decisions["merge_group"] == group, "cluster"
        ].astype(str).tolist()

        tasks.append({
            "type": "merge",
            "name": f"merge_{group}",
            "clusters": clusters
        })

    # --- SUBSET CLUSTERS ---
    subset_cluster_list = decisions.loc[
        decisions["action"] == "subset", "cluster"
    ].astype(str).tolist()

    for cl in subset_cluster_list:
        tasks.append({
            "type": "subset",
            "name": f"subset_{cl}",
            "clusters": [cl]
        })

    # --------------------------------------------------
    # Validate tasks
    # --------------------------------------------------
    if len(tasks) == 0:
        print("No refinement tasks found (no merge or subset). Nothing to refine.")
        return

    print(f"Found {len(tasks)} refinement task(s)")


    # --------------------------------------------------
    # Process each group
    # --------------------------------------------------
    for task in tasks:

        print(f"\nProcessing {task['name']} ({task['type']})")

        subset_clusters = task["clusters"]

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
        # Reset neighbours before final clustering
        # --------------------------------------------------
        sc.pp.neighbors(
            adata_sub,
            use_rep="X_pca_harmony",
            n_pcs=args.npcs,
            n_neighbors=args.neighbors
        )

        # --------------------------------------------------
        # Clustering
        # --------------------------------------------------
        sc.tl.leiden(
            adata_sub,
            resolution=res,
            key_added="leiden_refined",
            flavor="igraph",
            directed=False,
            n_iterations=2,
            random_state=0
        )

        # --------------------------------------------------
        # Stability (fast, label-invariant)
        # --------------------------------------------------
        if args.n_repeats > 1:

            print(f"Computing stability ({args.n_repeats} repeats)")

            assignments = []

            # Collect clustering runs
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

                assignments.append(
                    adata_sub.obs[f"tmp_{i}"].to_numpy(dtype=str)
                )

            assignments = np.array(assignments)

            # reference for mapping
            ref_labels = assignments[0]
            aligned_assignments = [ref_labels]

            for i in range(1, assignments.shape[0]):
                cur = assignments[i]

                contingency = pd.crosstab(cur, ref_labels)
                mapping = contingency.idxmax(axis=1).to_dict()

                cur_series = pd.Series(cur)
                aligned = cur_series.map(mapping).fillna(cur_series).to_numpy()

                aligned_assignments.append(aligned)

            aligned_assignments = np.array(aligned_assignments)

            # vectorised stability
            df_assign = pd.DataFrame(aligned_assignments.T)

            stability = df_assign.apply(
                lambda row: row.value_counts().max(),
                axis=1
            ).values / args.n_repeats

            adata_sub.obs["cluster_stability"] = stability

            # cleanup
            for i in range(args.n_repeats):
                del adata_sub.obs[f"tmp_{i}"]

        # --------------------------------------------------
        # Within-lineage markers
        # --------------------------------------------------
        print(f"Computing within-lineage markers for {task['name']}")

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
            print(f"No within-lineage DE results for {task['name']}, skipping export")
            continue

        marker_clusters = result["names"].dtype.names

        if marker_clusters is None or len(marker_clusters) == 0:
            print(f"No marker clusters found for {task['name']}, skipping export")
        else:
            markers_file = TABLE_DIR / f"{input_file.stem}_{task['name']}_within_lineage_markers.xlsx"

            with pd.ExcelWriter(markers_file, engine="openpyxl") as writer:

                if adata_de.raw is None:
                    raise ValueError("Expected .raw for marker extraction")

                X = adata_de.raw.X
                var_names = adata_de.raw.var_names
                var_index = pd.Series(range(len(var_names)), index=var_names)

                parent_label = "+".join(sorted(subset_clusters, key=int))

                sheets_written = 0

                for cl in marker_clusters:

                    genes = np.array(result["names"][cl]).ravel()
                    logfc_all = np.array(result["logfoldchanges"][cl]).ravel()
                    scores_all = np.array(result["scores"][cl]).ravel()
                    pvals_all = np.array(result["pvals"][cl]).ravel()
                    pvals_adj_all = np.array(result["pvals_adj"][cl]).ravel()

                    valid_mask = np.isin(genes, var_names)

                    valid = genes[valid_mask]
                    if valid.size == 0:
                        continue

                    logfc = logfc_all[valid_mask]
                    scores = scores_all[valid_mask]
                    pvals = pvals_all[valid_mask]
                    pvals_adj = pvals_adj_all[valid_mask]

                    if not (len(valid) == len(logfc) == len(scores) == len(pvals) == len(pvals_adj)):
                        continue

                    gene_idx = var_index.loc[valid].values

                    cluster_cells = adata_de.obs["leiden_refined"] == cl
                    other_cells = ~cluster_cells

                    expr = X[:, gene_idx]
                    expr = expr.tocsr()

                    pct_in = (expr[cluster_cells.values] > 0).mean(axis=0)
                    pct_out = (expr[other_cells.values] > 0).mean(axis=0)

                    pct_in = np.asarray(pct_in).ravel()
                    pct_out = np.asarray(pct_out).ravel()

                    mean_in = expr[cluster_cells.values].mean(axis=0)
                    mean_out = expr[other_cells.values].mean(axis=0)

                    mean_in = np.asarray(mean_in).ravel()
                    mean_out = np.asarray(mean_out).ravel()

                    if expr.shape[1] != len(valid):
                        continue

                    df = pd.DataFrame({
                        "gene": valid,
                        "logfoldchange": logfc,
                        "score": scores,
                        "pvals": pvals,
                        "pvals_adj": pvals_adj,
                    }).reset_index(drop=True)

                    df["pct_in"] = pct_in
                    df["pct_out"] = pct_out
                    df["pct_diff"] = pct_in - pct_out

                    df["mean_in"] = mean_in
                    df["mean_out"] = mean_out

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
                    print(f"No valid marker sheets written for {task['name']}, skipping file")

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
            continue

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
                var_index = pd.Series(range(len(var_names)), index=var_names)

                sheets_written = 0

                for cl in marker_clusters:

                    genes = np.array(result["names"][cl]).ravel()
                    logfc_all = np.array(result["logfoldchanges"][cl]).ravel()
                    scores_all = np.array(result["scores"][cl]).ravel()
                    pvals_all = np.array(result["pvals"][cl]).ravel()
                    pvals_adj_all = np.array(result["pvals_adj"][cl]).ravel()

                    valid_mask = np.isin(genes, var_names)

                    valid = genes[valid_mask]
                    if valid.size == 0:
                        continue

                    logfc = logfc_all[valid_mask]
                    scores = scores_all[valid_mask]
                    pvals = pvals_all[valid_mask]
                    pvals_adj = pvals_adj_all[valid_mask]

                    if not (len(valid) == len(logfc) == len(scores) == len(pvals) == len(pvals_adj)):
                        continue

                    gene_idx = var_index.loc[valid].values

                    cluster_cells = adata_de.obs["leiden_L2"] == cl
                    other_cells = ~cluster_cells

                    expr = X[:, gene_idx]
                    expr = expr.tocsr()

                    pct_in = (expr[cluster_cells.values] > 0).mean(axis=0)
                    pct_out = (expr[other_cells.values] > 0).mean(axis=0)

                    pct_in = np.asarray(pct_in).ravel()
                    pct_out = np.asarray(pct_out).ravel()

                    mean_in = expr[cluster_cells.values].mean(axis=0)
                    mean_out = expr[other_cells.values].mean(axis=0)

                    mean_in = np.asarray(mean_in).ravel()
                    mean_out = np.asarray(mean_out).ravel()

                    if expr.shape[1] != len(valid):
                        continue

                    df = pd.DataFrame({
                        "gene": valid,
                        "logfoldchange": logfc,
                        "score": scores,
                        "pvals": pvals,
                        "pvals_adj": pvals_adj,
                    }).reset_index(drop=True)

                    df["pct_in"] = pct_in
                    df["pct_out"] = pct_out
                    df["pct_diff"] = pct_in - pct_out

                    df["mean_in"] = mean_in
                    df["mean_out"] = mean_out

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
