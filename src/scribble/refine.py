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



def robust_hvg(sc, adata_sub, args):
    """
    Single-attempt HVG selection.
    If Seurat v3 HVG fails, treat cluster as terminal (no refinement).
    """

    try:
        if args.skip_integration:

            sc.pp.highly_variable_genes(
                adata_sub,
                n_top_genes=args.hvgs,
                flavor="seurat_v3",
                layer="counts"
            )

        else:

            sc.pp.highly_variable_genes(
                adata_sub,
                n_top_genes=args.hvgs,
                batch_key=args.batch,
                flavor="seurat_v3",
                layer="counts"
            )

        return True

    except Exception as e:

        print(
            f"[HVG FAILURE] {e}\n"
            "Refinement abandoned for this cluster."
        )

        return False


def is_marker_strong(marker_df, threshold):
    """
    Quick proxy: are top markers clearly separable?
    Uses logfoldchange as proxy for biological signal.
    """
    if marker_df is None or marker_df.empty:
        return False

    return marker_df["logfoldchange"].head(5).mean() >= threshold


def compute_markers(sc, pd, np, adata_de, args, groupby, min_cluster_size=None):

    if min_cluster_size is None:
        min_cluster_size = args.min_cells_per_cluster

    cluster_sizes = (
        adata_de.obs[groupby]
        .value_counts()
    )

    valid = cluster_sizes[
        cluster_sizes >= min_cluster_size
    ].index

    if len(valid) < 2:

        print(
            f"[Markers] Skipping {groupby}: "
            f"fewer than two groups with >={min_cluster_size} cells"
        )

        return {}

    adata_de = adata_de[adata_de.obs[groupby].isin(valid)].copy()

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


def run_clustering(np, sc, pd, adata_sub, args):
    from scribble.cluster import optimise_resolution

    sc.pp.filter_genes(
        adata_sub,
        min_cells=args.min_cells_per_gene
    )

    print(
        f"[HVG input] cells={adata_sub.n_obs}, "
        f"genes={adata_sub.n_vars}"
    )

    if args.batch in adata_sub.obs:
        print(
            "[HVG batches]",
            adata_sub.obs[args.batch].value_counts().to_dict()
        )

    gene_ncells = np.asarray(
        (adata_sub.X > 0).sum(axis=0)
    ).ravel()

    print(
        f"[HVG genes] "
        f"n=3: {(gene_ncells == 3).sum()}, "
        f"n<=5: {(gene_ncells <= 5).sum()}, "
        f"n<=10: {(gene_ncells <= 10).sum()}"
    )

    # preserve raw counts
    adata_sub.layers["counts"] = adata_sub.X.copy()

    # ------------------------------------------
    # HVGs must be selected on counts for
    # flavor="seurat_v3"
    # ------------------------------------------
    hvg_ok = robust_hvg(sc, adata_sub, args)

    if not hvg_ok:
        return None

    # ------------------------------------------
    # Normalisation for PCA / Harmony / DE
    # ------------------------------------------
    sc.pp.normalize_total(adata_sub)
    sc.pp.log1p(adata_sub)

    # store log-normalised matrix for DE
    adata_sub.raw = adata_sub.copy()

    # restrict to HVGs after they have been selected
    adata_sub = adata_sub[:, adata_sub.var.highly_variable].copy()

    if not args.no_scale:
        sc.pp.scale(adata_sub, max_value=10)

    sc.tl.pca(adata_sub, n_comps=args.npcs)

    if args.skip_integration:

        use_rep = "X_pca"

    else:

        ho = hm.run_harmony(
            adata_sub.obsm["X_pca"],
            adata_sub.obs,
            args.batch,
            theta=args.theta
        )

        adata_sub.obsm["X_pca_harmony"] = ho.Z_corr
        use_rep = "X_pca_harmony"

    sc.pp.neighbors(
        adata_sub,
        use_rep=use_rep,
        n_neighbors=args.neighbors
    )

    if args.auto_resolution:
        res, _, _ = optimise_resolution(
            np,
            pd,
            sc,
            adata_sub,
            use_rep,
            args.neighbors,
            (args.res_min, args.res_max),
            args.fine_width,
            args.res_steps
        )
    else:
        res = args.resolution

    sc.pp.neighbors(
        adata_sub,
        use_rep=use_rep,
        n_neighbors=args.neighbors
    )

    sc.tl.leiden(
        adata_sub,
        resolution=res,
        key_added="leiden_refined",
        flavor="igraph",
        directed=False,
        n_iterations=2,
        random_state=0
    )

    if args.n_repeats > 1:

        stability_results = compute_cluster_stability(
            adata=adata_sub,
            resolution=res,
            n_repeats=args.n_repeats,
        )

        adata_sub.obs["cluster_stability"] = (
            stability_results["cell_stability"]
        )

        adata_sub.uns["stability_aligned_assignments"] = (
            stability_results["aligned_assignments"]
        )

        reassignment_df = (
            summarise_cluster_reassignments(
                final_labels=(
                    adata_sub.obs["leiden_refined"]
                    .astype(str)
                    .values
                ),
                aligned_assignments=(
                    stability_results["aligned_assignments"]
                ),
            )
        )

        adata_sub.uns["cluster_reassignments"] = (
            reassignment_df
        )

    return adata_sub


def run_refine(args):
    import scanpy as sc
    import pandas as pd
    import numpy as np
    from scipy import sparse
    import random
    import harmonypy as hm
    from pathlib import Path
    import matplotlib.pyplot as plt

    from scribble.import_data import setup_environment
    from scipy.stats import entropy
    from scribble.utils import (
        write_run_manifest,
        compute_cluster_stability,
        summarise_cluster_reassignments,
        compute_cluster_summary,
    )

    PROJECT_DIR = Path(args.project_dir)
    MANIFEST_DIR = PROJECT_DIR / "scribble/manifests"
    PLOT_DIR = PROJECT_DIR / "scribble/plots"
    TABLE_DIR = PROJECT_DIR / "scribble/tables"

    setup_environment(sc, np, random, PLOT_DIR)

    input_file = Path(args.input)
    manifest_file = (MANIFEST_DIR / f"{input_file.stem}_refine_manifest.json")
    reassignments_file = TABLE_DIR / "refine_cluster_reassignments.tsv"
    all_markers_file = TABLE_DIR / "all_markers.tsv"
    refinement_report_file = TABLE_DIR / "refinement_report.tsv"
    umap_png_file = PLOT_DIR / "UMAP_refine_cluster.png"

    decisions = pd.read_csv(args.decisions, sep="\t")

    print(f"Loading AnnData → {input_file}")
    adata = sc.read(input_file)

    if (
        not args.ignore_saved_params and
        "scribble_cluster_params" in adata.uns
    ):

        saved = adata.uns["scribble_cluster_params"]

        print(
            "\nUsing clustering parameters "
            "stored in AnnData:"
        )

        mapping = [
            "neighbors",
            "n_repeats",
            "res_min",
            "res_max",
            "res_steps",
            "fine_width",
            "dominant_reassignment_threshold",
            "reassignment_display_threshold",
            "high_reassignment",
        ]

        for param in mapping:

            if param in saved:

                setattr(
                    args,
                    param,
                    saved[param]
                )

                print(
                    f"  {param} = "
                    f"{saved[param]}"
                )

    # --------------------------------------------------
    # Init labels + refinement tree
    # --------------------------------------------------
    adata.obs["refine_label"] = adata.obs["leiden"].astype(str)

    # Each refine run constructs a new tree from scratch.
    adata.uns["refinement_tree"] = {}

    # Refine cluster registry
    next_cluster_id = 1
    cluster_registry = []

    # Global marker registry
    master_markers = {}
    master_marker_rows = []
    master_reassignments = []

    # --------------------------------------------------
    # Build initial tasks
    # --------------------------------------------------

    tasks = []

    refine_df = decisions.copy()

    refine_df["cluster"] = refine_df["cluster"].astype(str)

    # ----------------------------------------
    # 1. Handle merge groups
    # ----------------------------------------

    processed = set()

    if "merge_group" in refine_df.columns:

        merge_df = decisions[
            decisions["merge_group"].notna()
            & (decisions["merge_group"].astype(str).str.strip() != "")
        ]

        print("\n[MERGE DF]")
        print(
            merge_df[
                [
                    "cluster",
                    "action",
                    "merge_group",
                    "suggested_merge_target",
                ]
            ]
        )

        print("\n[MERGE GROUP CONTENTS]")

        for mg, grp in merge_df.groupby("merge_group"):

            print(
                mg,
                grp["cluster"]
                .astype(str)
                .tolist()
            )

        for merge_group, grp in merge_df.groupby("merge_group"):

            clusters = (
                grp["cluster"]
                .astype(str)
                .tolist()
            )

            tasks.append({
                "clusters": clusters,
                "level": 2,
            })

            print(
                f"[TASK MERGE] {merge_group}: "
                f"clusters={clusters}"
            )


            processed.update(clusters)

            print(
                f"[Refine] Merge group '{merge_group}' "
                f"→ clusters {clusters}"
            )

    # ----------------------------------------
    # 2. Handle standalone subset clusters
    # ----------------------------------------

    for cl in refine_df.loc[
        refine_df["action"] == "subset",
        "cluster"
    ]:

        if cl in processed:
            continue

        tasks.append({
            "clusters": [cl],
            "level": 2,
        })

    print("\n[INITIAL TASKS]")
    for t in tasks:
        print(t)
    print()

    if len(tasks) == 0:
        print("No refinement tasks.")

    # --------------------------------------------------
    # INTERNAL FUNCTIONS
    # --------------------------------------------------

    def _entropy_row(row):
        vals = row.to_numpy(dtype=float)
        probs = vals / vals.sum()
        return entropy(probs)


    def _refine_task(task):
        nonlocal next_cluster_id

        clusters = task["clusters"]
        level = task["level"]

        print(
            f"\n[REFINE TASK START] "
            f"level={level} "
            f"clusters={clusters}"
        )

        mask = (
            adata.obs["refine_label"]
            .astype(str)
            .isin([str(x) for x in clusters])
        )

        print(
            f"[TASK SELECT] clusters={clusters} "
            f"matched_cells={mask.sum()}"
        )
        adata_sub = adata[mask].copy()

        print(f"[Refine] Level {level} | clusters = {clusters} | cells = {adata_sub.n_obs}")

        if adata_sub.n_obs < args.min_cells_per_group:
            return []

        adata_sub = restore_counts(adata_sub)

        adata_sub = run_clustering(np, sc, pd, adata_sub, args)

        if adata_sub is not None:

            print(
                "[REFINE RESULT] "
                f"clusters={clusters} "
                f"resolved_clusters="
                f"{adata_sub.obs['leiden_refined'].nunique()}"
            )

            print(
                adata_sub.obs["leiden_refined"]
                .value_counts()
                .sort_index()
            )

        if adata_sub is None:
            parent_label = (
                clusters[0]
                if len(clusters) == 1
                else "_".join(sorted(clusters))
            )

            # Recover existing ancestry if this is a recursive refinement
            parent_info = adata.uns["refinement_tree"].get(parent_label, {})

            if isinstance(parent_info, dict):
                parent_path = parent_info.get("path", [parent_label])
            else:
                parent_path = [parent_label]

            adata.uns["refinement_tree"][parent_label] = {
                "parent": parent_info.get("parent")
                if isinstance(parent_info, dict) else None,
                "level": level,
                "path": parent_path,
                "node_type": "refinement_failed"
            }

            print(
                f"[Refine] Level {level} | clusters={clusters} "
                "→ TERMINAL (refinement failed)"
            )
            return []


        # --------------------------------------------------
        # Parent node markers
        # --------------------------------------------------

        if len(clusters) == 1:

            parent_name = (
                f"L{level-1}_"
                + clusters[0]
            )

        else:

            parent_name = None

        if parent_name is not None:

            try:

                parent_df = pd.DataFrame()

                expr = adata_sub.X

                mean_expr = np.asarray(
                    expr.mean(axis=0)
                ).ravel()

                parent_df["gene"] = adata_sub.var_names
                parent_df["mean_expr"] = mean_expr

                parent_df = (
                    parent_df
                    .sort_values(
                        "mean_expr",
                        ascending=False
                    )
                    .head(args.nmarkers)
                    .copy()
                )


                tmp = parent_df.copy()

                tmp["refine_cluster"] = None
                tmp["parent_refine_cluster"] = (
                    clusters[0]
                    if len(clusters) == 1
                    else "_".join(sorted(clusters))
                )
                tmp["level"] = level - 1
                tmp["node_type"] = "parent"
                tmp["n_parent_clusters"] = len(clusters)
                tmp["marker_type"] = "parent_mean_expression"

                master_marker_rows.append(tmp)

            except Exception as e:

                print(
                    f"[Parent markers] "
                    f"{parent_name}: {e}"
                )



        # -----------------------
        # Map refined labels
        # -----------------------
        cluster_map = {}

        cluster_sizes = (
            adata_sub.obs["leiden_refined"]
            .astype(str)
            .value_counts()
        )

        for cl in sorted(
            adata_sub.obs["leiden_refined"].astype(str).unique()
        ):

            n_cells = int(cluster_sizes[cl])

            # too small -> collapse into parent
            if n_cells < args.min_cells_per_cluster:

                cluster_map[cl] = None

                print(
                    f"[Collapse] refined cluster {cl} "
                    f"({n_cells} cells) "
                    f"returned to parent"
                )

                continue

            cluster_id = f"R{next_cluster_id}"

            cluster_map[cl] = cluster_id

            cluster_registry.append({
                "refine_cluster": cluster_id,
                "parent_refine_cluster": (
                    clusters[0]
                    if len(clusters) == 1
                    else "_".join(sorted(clusters))
                ),
                "local_cluster": cl,
                "level": level,
                "n_cells": n_cells
            })

            next_cluster_id += 1

        refined = (
            adata_sub.obs["leiden_refined"]
            .astype(str)
            .map(cluster_map)
        )

        print(
            "[REFINE MAP]",
            cluster_map
        )

        if "cluster_reassignments" in adata_sub.uns:

            tmp = (
                adata_sub.uns[
                    "cluster_reassignments"
                ]
                .copy()
            )

            # local Leiden -> global refine label
            global_map = {
                str(k): v
                for k, v in cluster_map.items()
                if v is not None
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

            tmp = tmp.dropna(
                subset=["cluster", "destination"]
            )

            tmp["level"] = level

            tmp["parent_clusters"] = (
                "|".join(clusters)
            )

            master_reassignments.append(tmp)

        # any cluster collapsed for being too small
        collapsed = refined.isna()

        if collapsed.any():

            largest_refined = (
                pd.Series(cluster_map)
                .dropna()
                .iloc[0]
            )

            refined.loc[collapsed] = largest_refined

            print(
                f"[COLLAPSE MERGE] "
                f"{collapsed.sum()} cells -> "
                f"{largest_refined}"
            )

        adata.obs.loc[
            adata_sub.obs_names,
            "refine_label"
        ] = refined.values

        print(
            f"[ASSIGN] parent={clusters}"
        )

        print(
            pd.Series(refined)
            .value_counts()
            .sort_index()
        )

        if "cluster_stability" in adata_sub.obs:
            adata.obs.loc[
                adata_sub.obs_names,
                "cluster_stability"
            ] = (
                adata_sub.obs["cluster_stability"]
            )

        # -----------------------
        # Refinement tracking
        # -----------------------
        children = refined.unique().tolist()

        for child in children:

            # Which current parent(s) produced this child?
            parent_label = (
                clusters[0]
                if len(clusters) == 1
                else "_".join(sorted(clusters))
            )

            # Recover ancestry already associated with the parent
            parent_info = adata.uns["refinement_tree"].get(
                parent_label,
                {}
            )

            if isinstance(parent_info, dict):
                parent_path = parent_info.get(
                    "path",
                    [parent_label]
                )
            else:
                # Backward compatibility with old tree structure
                parent_path = [parent_label]

            # Do not append parent again when a collapsed cluster
            # retains the parent label
            if child == parent_label:
                child_path = parent_path.copy()
            else:
                child_path = parent_path + [child]

            # A collapsed group that retains the parent label is a root/terminal
            # node, not its own child.
            if child == parent_label:
                adata.uns["refinement_tree"][child] = {
                    "parent": None,
                    "level": level,
                    "path": child_path,
                    "node_type": "retained_parent",
                    "parent_clusters": clusters.copy(),
                    "parent_cluster_type": (
                        "merge_group"
                        if len(clusters) > 1
                        else "single_cluster"
                    ),
                    "refinement_source": (
                        "merge_group"
                        if len(clusters) > 1
                        else "single_cluster"
                    ),
                }

            else:
                adata.uns["refinement_tree"][child] = {
                    "parent": parent_label,
                    "level": level,
                    "path": child_path,
                    "node_type": "refined_child",
                    "parent_clusters": clusters.copy(),
                    "parent_cluster_type": (
                        "merge_group"
                        if len(clusters) > 1
                        else "single_cluster"
                    ),
                    "refinement_source": (
                        "merge_group"
                        if len(clusters) > 1
                        else "single_cluster"
                    ),
                }

        # -----------------------
        # Markers
        # -----------------------
        cluster_sizes = adata_sub.obs["leiden_refined"].value_counts()
        valid = cluster_sizes[cluster_sizes >= 2].index
        marker_tables = {}

        if len(valid) >= 2:

            adata_de = adata_sub[
                adata_sub.obs["leiden_refined"].isin(valid)
            ].copy()

            marker_tables = compute_markers(
                sc, pd, np, adata_de, args,
                groupby = "leiden_refined"
            )

            clusters_str = "+".join(task["clusters"])

            # Save workbook for this refinement step
            out_file = (
                TABLE_DIR
                / f"L{task['level']}_{clusters_str}_markers.xlsx"
            )

            with pd.ExcelWriter(out_file) as writer:

                for cl, df in marker_tables.items():

                    # Refine node label
                    node_label = cluster_map.get(str(cl))

                    if node_label is None:
                        continue

                    # Workbook storage
                    master_markers[node_label] = df.copy()

                    df.to_excel(
                        writer,
                        sheet_name=node_label,
                        index=False
                    )

                    # Flat table storage
                    tmp = df.copy()

                    tmp["refine_cluster"] = node_label

                    tmp["parent_clusters"] = (
                        "|".join(clusters)
                    )

                    tmp["level"] = task["level"]

                    tmp["node_type"] = "child"
                    tmp["n_parent_clusters"] = len(clusters)
                    tmp["marker_type"] = "differential"

                    master_marker_rows.append(tmp)

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
                    adata_sub.obs
                    .groupby("leiden_refined",observed=False)["cluster_stability"]
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
                    if is_marker_strong(marker_df, args.marker_strength_threshold):
                        print("    → STOP (already well separated)")
                        continue

                # ------------------------------
                # OTHERWISE → REFINE
                # ------------------------------
                print("    → REFINE further")

                child_labels = refined[
                    adata_sub.obs["leiden_refined"].astype(str) == cluster_label
                ].unique().tolist()

                print(
                    f"[RECURSION CHECK] local={cluster_label} "
                    f"global={child_labels}"
                )

                new_tasks.append({
                    "clusters": child_labels,
                    "level": level + 1
                })

        return new_tasks


    # --------------------------------------------------
    # MAIN LOOP (safe recursion)
    # --------------------------------------------------


    if len(tasks) > 0:

        i = 0
        while i < len(tasks):

            print(tasks[i])
            new_tasks = _refine_task(tasks[i])
            tasks.extend(new_tasks)

            i += 1

    # --------------------------------------------------
    # Store final refine IDs
    # --------------------------------------------------
    print("\n[FINAL LABEL COUNTS]")
    print(
        adata.obs["refine_label"]
        .astype(str)
        .value_counts()
        .sort_index()
    )

    final_labels = sorted(
        adata.obs["refine_label"]
        .astype(str)
        .unique(),
        key=str
    )

    final_registry = []

    for label in final_labels:

        n_cells = int(
            (adata.obs["refine_label"].astype(str) == label).sum()
        )

        final_registry.append({
            "refine_cluster": label,
            "n_cells": n_cells
        })

    cluster_registry = final_registry


    adata.obs["leiden_L2"] = (
        adata.obs["refine_label"]
        .astype(str)
    )

    print(
        "Final refine clusters:",
        adata.obs["refine_label"].nunique()
    )


    # Identify global markers
    adata_global = restore_counts(adata)
    if "X_umap" in adata.obsm:
        adata_global.obsm["X_umap"] = adata.obsm["X_umap"].copy()
    sc.pp.normalize_total(adata_global)
    sc.pp.log1p(adata_global)
    adata_global.raw = adata_global.copy()
    adata_global.obs["refine_cluster"] = (
        adata.obs["refine_label"]
        .astype(str)
    )

    refine_reassignments = None

    if len(master_reassignments) > 0:

        refine_reassignments = pd.concat(
            master_reassignments,
            ignore_index=True,
        )

        refine_reassignments.to_csv(
            reassignments_file,
            sep="\t",
            index=False,
        )

    refine_summary = compute_cluster_summary(
        adata=adata_global,
        cluster_key="refine_cluster",
        stability_key="cluster_stability",
        reassignment_df=refine_reassignments,
        dominant_reassignment_threshold=args.dominant_reassignment_threshold,
        reassignment_display_threshold=args.reassignment_display_threshold,
        merge_candidate_threshold=args.high_reassignment
    )

    global_markers = compute_markers(
        sc, pd, np, adata_global, args,
        groupby = "refine_cluster",
        min_cluster_size=2
    )

    print(
        "Global marker clusters:",
        len(global_markers)
    )

    print(
        "Registry clusters:",
        len(cluster_registry)
    )


    # --------------------------------------------------
    # Final global markers
    # --------------------------------------------------
    master_file = (
        TABLE_DIR
        / "L2_markers.xlsx"
    )

    with pd.ExcelWriter(master_file) as writer:

        # INDEX worksheet first
        pd.DataFrame(
            cluster_registry
        ).to_excel(
            writer,
            sheet_name="INDEX",
            index=False
        )

        # Marker sheets
        def _cluster_sort_key(row):
            label = str(row["refine_cluster"])

            # Original Leiden cluster
            if label.isdigit():
                return (0, int(label))

            # Refinement-derived cluster
            if label.startswith("R") and label[1:].isdigit():
                return (1, int(label[1:]))

            # Composite/retained parent such as 5|8|15|6|16
            return (2, label)

        all_clusters = sorted(
            cluster_registry,
            key=_cluster_sort_key
        )

        for row in all_clusters:

            label = str(row["refine_cluster"])

            if label in global_markers:

                df = global_markers[label]

            else:

                df = pd.DataFrame({
                    "note": [
                        "No markers available"
                    ]
                })

            df.to_excel(
                writer,
                sheet_name=label,
                index=False
            )



    # --------------------------------------------------
    # Flat marker table
    # --------------------------------------------------

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

        print(f"Saved {len(all_markers):,} marker rows to {all_markers_file}")

    # --------------------------------------------------
    # Save refinement summary / tree
    # --------------------------------------------------

    refinement_rows = []

    for cluster, info in adata.uns["refinement_tree"].items():

        if not isinstance(info, dict):
            continue

        path = info.get("path", [cluster])

        refinement_rows.append({
            "cluster": cluster,
            "parent": info.get("parent"),
            "level": info.get("level"),
            "refinement_path": " -> ".join(map(str, path)),
            "refinement_source": info.get(
                "refinement_source"
            ),
            "parent_clusters": "|".join(
                map(
                    str,
                    info.get(
                        "parent_clusters",
                        []
                    )
                )
            ),
            "n_parent_clusters": len(
                info.get(
                    "parent_clusters",
                    []
                )
            ),
        })

    if len(refinement_rows) > 0:

        refinement_df = (
            pd.DataFrame(refinement_rows)
            .drop_duplicates()
        )

        cluster_sizes = (
            adata.obs["refine_label"]
            .astype(str)
            .value_counts()
        )

        refinement_df["n_cells"] = (
            refinement_df["cluster"]
            .map(cluster_sizes)
        )

        # ------------------------------
        # Merge with refinement summary
        # ------------------------------
        refinement_report = refinement_df.merge(
            refine_summary.drop(columns=["n_cells"]),
            on="cluster",
            how="left",
        )

        refinement_report = refinement_report[
            [
                "cluster",
                "n_cells",

                "parent",
                "level",
                "refinement_path",

                "refinement_source",
                "parent_clusters",
                "n_parent_clusters",

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

        refinement_report.to_csv(refinement_report_file, sep="\t", index=False)

        print(f"Saved {len(refinement_report):,} refinement records to {refinement_report_file}")

    # --------------------------------------------------
    # UMAP final cluster labels
    # --------------------------------------------------
    sc.pl.umap(
        adata_global,
        color="refine_cluster",
        legend_loc="right margin",
        legend_fontsize=7,
        show=False
    )
    plt.savefig(umap_png_file, dpi=300, bbox_inches="tight")
    plt.close()


    # --------------------------------------------------
    # Save
    # --------------------------------------------------
    adata.uns["scribble_manifest"] = str(manifest_file)
    adata.uns["scribble_refine_params"] = (vars(args).copy())
    output = input_file.with_name(f"{input_file.stem}_refined.h5ad")
    print(f"Saving → {output}")
    adata.write(output)

    print(f"Writing run manifest to {manifest_file}")
    write_run_manifest(
        manifest_file=manifest_file,
        module="refine",
        args=args,
        inputs={
            "input_file": str(input_file),
            "tsv": str(args.decisions),
        },
        outputs={
            "adata": str(output),
            "png": str(umap_png_file),
            "tsv": str(reassignments_file),
            "tsv": str(all_markers_file),
            "tsv": str(refinement_report_file),
        },
    )
