#!/usr/bin/env python

from pathlib import Path

from seaborn.relational import _Markers

def optimise_resolution(np, pd, sc, adata, embedding, neighbors, coarse_range, fine_width, n_steps):
    from sklearn.metrics import silhouette_score

    # --------------------------------------------------
    # COARSE PASS
    # --------------------------------------------------
    print("Running coarse resolution scan...")

    coarse_resolutions = np.linspace(coarse_range[0], coarse_range[1], n_steps)
    coarse_results = []

    for res in coarse_resolutions:
        sc.tl.leiden(adata, resolution=res, key_added="leiden_tmp",
            flavor="igraph", directed=False, n_iterations=2)

        labels = adata.obs["leiden_tmp"].astype(int)
        n_clusters = labels.nunique()

        if n_clusters < 2:
            sil = -1
        else:
            sil = silhouette_score(adata.obsm[embedding], labels)

        coarse_results.append({
            "resolution": res,
            "silhouette": sil,
            "n_clusters": n_clusters
        })

    coarse_df = pd.DataFrame(coarse_results)

    best_coarse = coarse_df.loc[coarse_df["silhouette"].idxmax(), "resolution"]
    print(f"Best coarse resolution: {best_coarse:.3f}")

    # --------------------------------------------------
    # FINE PASS
    # --------------------------------------------------
    print("Running fine resolution scan...")

    fine_min = max(0.01, best_coarse - fine_width)
    fine_max = best_coarse + fine_width
    fine_resolutions = np.linspace(fine_min, fine_max, n_steps)

    fine_results = []

    for res in fine_resolutions:
        sc.tl.leiden(adata, resolution=res, key_added="leiden_tmp",
            flavor="igraph", directed=False, n_iterations=2)

        labels = adata.obs["leiden_tmp"].astype(int)
        n_clusters = labels.nunique()

        if n_clusters < 2:
            sil = -1
        else:
            sil = silhouette_score(adata.obsm[embedding], labels)

        fine_results.append({
            "resolution": res,
            "silhouette": sil,
            "n_clusters": n_clusters
        })

    fine_df = pd.DataFrame(fine_results)

    best_fine = fine_df.loc[fine_df["silhouette"].idxmax(), "resolution"]
    print(f"Optimal resolution: {best_fine:.3f}")

    return best_fine, coarse_df, fine_df


def run_cluster(args):
    import scanpy as sc
    import numpy as np
    import pandas as pd
    import matplotlib.pyplot as plt
    from scipy.optimize import linear_sum_assignment
    from scipy.stats import entropy, mode
    import random
    from scribble.import_data import setup_environment

    PROJECT_DIR = Path(args.project_dir)
    PLOT_DIR = PROJECT_DIR / "sc_plots"
    setup_environment(sc, np, random, PLOT_DIR)

    input_file = Path(args.input)
    output_file = input_file.with_name(f"{input_file.stem}_clustered.h5ad")
    markers_file = input_file.with_name(f"{input_file.stem}_clusters.xlsx")

    print(f"Loading {input_file}")
    adata = sc.read(input_file)

    # --------------------------------------------------
    # Validation
    # --------------------------------------------------
    if args.embedding not in adata.obsm:
        raise ValueError(f"{args.embedding} not found in adata.obsm")

    if "sample" not in adata.obs.columns:
        raise ValueError("sample column required for cluster diagnostics")

    # --------------------------------------------------
    # Build neighbors
    # --------------------------------------------------
    print(f"Building neighbors (embedding={args.embedding})")

    sc.pp.neighbors(
        adata,
        use_rep=args.embedding,
        n_neighbors=args.neighbors
    )

    # --------------------------------------------------
    # Resolution selection
    # --------------------------------------------------
    if args.auto_resolution:
        best_res, coarse_df, fine_df = optimise_resolution(
            np, pd, sc,
            adata,
            args.embedding,
            args.neighbors,
            (args.res_min, args.res_max),
            args.fine_width,
            args.res_steps
        )

        args.resolution = best_res

        # Save diagnostics
        coarse_file = input_file.with_name(f"{input_file.stem}_res_coarse.tsv")
        fine_file = input_file.with_name(f"{input_file.stem}_res_fine.tsv")

        coarse_df.to_csv(coarse_file, sep="\t", index=False)
        fine_df.to_csv(fine_file, sep="\t", index=False)

        # --------------------------------------------------
        # Plot optimisation diagnostics
        # --------------------------------------------------
        plot_file = PLOT_DIR / f"{input_file.stem}_resolution_optimisation.png"

        fig, ax1 = plt.subplots(figsize=(7, 4))

        # --- silhouette axis ---
        ax1.set_xlabel("Resolution")
        ax1.set_ylabel("Silhouette score")
        ax1.plot(
            coarse_df["resolution"],
            coarse_df["silhouette"],
            marker="o",
            label="coarse (silhouette)"
        )
        ax1.plot(
            fine_df["resolution"],
            fine_df["silhouette"],
            marker="o",
            label="fine (silhouette)"
        )

        ax1.axvline(args.resolution, linestyle="--", label=f"selected={args.resolution:.2f}")
        ax1.grid(alpha=0.3)

        # --- cluster count axis ---
        ax2 = ax1.twinx()
        ax2.set_ylabel("Number of clusters")

        ax2.plot(
            coarse_df["resolution"],
            coarse_df["n_clusters"],
            linestyle="--",
            label="coarse (clusters)"
        )
        ax2.plot(
            fine_df["resolution"],
            fine_df["n_clusters"],
            linestyle="--",
            label="fine (clusters)"
        )

        # --- legend handling ---
        lines_1, labels_1 = ax1.get_legend_handles_labels()
        lines_2, labels_2 = ax2.get_legend_handles_labels()

        ax1.legend(lines_1 + lines_2, labels_1 + labels_2, loc="best")

        plt.title("Leiden Resolution Optimisation")
        plt.savefig(plot_file, dpi=300, bbox_inches="tight")
        plt.close()

        print(f"Saved optimisation plot → {plot_file}")

        if "leiden_tmp" in adata.obs:
            del adata.obs["leiden_tmp"]


    # --------------------------------------------------
    # Final clustering
    # --------------------------------------------------
    print(f"Running Leiden clustering (resolution={args.resolution})")

    sc.tl.leiden(
        adata,
        resolution=args.resolution,
        key_added="leiden",
        flavor="igraph", directed=False, n_iterations=2
    )

    print("Cluster sizes:")
    print(adata.obs["leiden"].value_counts())

    # --------------------------------------------------
    # Stability via Hungarian alignment
    # --------------------------------------------------
    if args.n_repeats > 1:
        print(f"Computing cluster stability ({args.n_repeats} repeats)...")

        # reference clustering
        ref_labels = adata.obs["leiden"].astype(str).values

        all_assignments = []

        def align_labels(reference, current):
            ref = pd.Series(reference)
            cur = pd.Series(current)

            contingency = pd.crosstab(cur, ref)

            # convert to cost matrix (maximize overlap → minimize negative)
            cost_matrix = -contingency.values

            row_ind, col_ind = linear_sum_assignment(cost_matrix)

            mapping = {
                contingency.index[row]: contingency.columns[col]
                for row, col in zip(row_ind, col_ind)
            }

            aligned = cur.map(mapping)

            # fallback: keep original label if not mapped
            aligned = aligned.fillna(cur)

            return aligned.values

        for i in range(args.n_repeats):
            sc.tl.leiden(
                adata,
                resolution=args.resolution,
                key_added=f"leiden_tmp_{i}",
                flavor="igraph", directed=False, n_iterations=2
            )

            raw_labels = adata.obs[f"leiden_tmp_{i}"].astype(str).values
            aligned = align_labels(ref_labels, raw_labels)

            all_assignments.append(aligned)

        all_assignments = np.array(all_assignments, dtype=str)  # (runs, cells)

        # compute stability
        stability = np.array([
            np.unique(all_assignments[:, i], return_counts=True)[1].max() / args.n_repeats
            for i in range(all_assignments.shape[1])
        ])

        adata.obs["cluster_stability"] = np.array(stability, dtype=float)

        print(f"Mean stability: {np.mean(stability):.3f}")

        # cleanup temp columns
        for i in range(args.n_repeats):
            del adata.obs[f"leiden_tmp_{i}"]

    # --------------------------------------------------
    # Export stats
    # --------------------------------------------------
    if args.n_repeats > 1:

        stats_file = input_file.with_name(f"{input_file.stem}_cluster_summary.tsv")
        cell_file = input_file.with_name(f"{input_file.stem}_cluster_cells.tsv")

        print(f"Exporting cluster stats → {stats_file}")

        adata.obs[["leiden", "cluster_stability"]].to_csv(cell_file, sep="\t")

        cluster_summary = (
            adata.obs
            .groupby("leiden", observed=True)
            .agg(
                n_cells=("leiden", "size"),
                mean_stability=("cluster_stability", "mean"),
                median_stability=("cluster_stability", "median")
            )
        )

        cluster_summary["fraction"] = cluster_summary["n_cells"] / adata.n_obs
        cluster_summary.index = cluster_summary.index.astype(str)

        sample_counts = (
            adata.obs
            .groupby(["leiden", "sample"], observed=True)
            .size()
            .unstack(fill_value=0)
        )

        sample_counts.index = sample_counts.index.astype(str)
        sample_counts.columns = sample_counts.columns.astype(str)

        # reset indices → convert to columns
        cluster_summary = cluster_summary.reset_index()
        sample_counts = sample_counts.reset_index()

        # ensure string keys
        cluster_summary["leiden"] = cluster_summary["leiden"].astype(str)
        sample_counts["leiden"] = sample_counts["leiden"].astype(str)

        # merge
        cluster_summary = pd.merge(
            cluster_summary,
            sample_counts,
            on="leiden",
            how="outer"
        ).fillna(0)

        cluster_summary.rename(columns={"leiden": "cluster"}, inplace=True)

        # --------------------------------------------------
        # Compute sample entropy (mixing metric)
        # --------------------------------------------------

        # Identify sample columns (all columns that came from sample_counts)
        exclude_cols = {
            "cluster",
            "n_cells",
            "mean_stability",
            "median_stability",
            "fraction",
            "resolution",
            "embedding",
            "n_repeats",
            "sample_entropy",
            "low_stability",
            "low_mixing"
        }

        sample_cols = [col for col in cluster_summary.columns if col not in exclude_cols]

        def compute_entropy_row(row):
            values = row[sample_cols].to_numpy(dtype=float)

            total = values.sum()
            if total == 0:
                return 0.0

            probs = values / total
            return entropy(probs)

        cluster_summary["sample_entropy"] = cluster_summary.apply(compute_entropy_row, axis=1)

        # flags
        cluster_summary["low_stability"] = cluster_summary["mean_stability"] < 0.7
        cluster_summary["low_mixing"] = cluster_summary["sample_entropy"] < 0.5

        # metadata
        cluster_summary["resolution"] = args.resolution
        cluster_summary["embedding"] = args.embedding
        cluster_summary["n_repeats"] = args.n_repeats

        cluster_summary.reset_index(drop=True, inplace=True)

        cluster_summary = cluster_summary.sort_values("n_cells", ascending=False)

        cluster_summary.to_csv(stats_file, sep="\t", index=False)

    # --------------------------------------------------
    # UMAP
    # --------------------------------------------------
    print("Computing UMAP...")
    sc.tl.umap(adata)

    # --------------------------------------------------
    # Plotting
    # --------------------------------------------------
    print("Generating plots...")

    cluster_file = PLOT_DIR / f"{input_file.stem}_clusters.png"

    sc.pl.umap(
        adata,
        color=["leiden"] + args.vars,
        wspace=0.4,
        show=False
    )
    plt.savefig(cluster_file, dpi=300, bbox_inches="tight")
    plt.close()

    if "cluster_stability" in adata.obs:
        stability_file = PLOT_DIR / f"{input_file.stem}_stability.png"
        sc.pl.umap(
            adata,
            color="cluster_stability",
            cmap="viridis",
            show=False
        )
        plt.savefig(stability_file, dpi=300, bbox_inches="tight")
        plt.close()

    # --------------------------------------------------
    # Extract marker genes
    # --------------------------------------------------
    print("Extracting markers...")
    # Identify marker genes
    sc.tl.rank_genes_groups(
        adata,
        "leiden",
        method="wilcoxon",
        use_raw=True
    )

    result = adata.uns["rank_genes_groups"]
    clusters = result["names"].dtype.names

    print(f"Writing markers → {markers_file}")
    with pd.ExcelWriter(_Markers_file, engine="openpyxl") as writer:
        for cl in clusters:

            genes = result["names"][cl]

            df = pd.DataFrame({
                "gene": genes,
                "logfoldchange": result["logfoldchanges"][cl],
                "score": result["scores"][cl],
                "pvals": result["pvals"][cl],
                "pvals_adj": result["pvals_adj"][cl],
            })

            # --------------------------------------------------
            # Expression statistics
            # --------------------------------------------------
            cluster_cells = adata.obs["leiden"] == cl
            other_cells = ~cluster_cells

            # Use raw if available
            X = adata.raw.X if adata.raw is not None else adata.X
            var_names = adata.raw.var_names if adata.raw is not None else adata.var_names

            # Convert to dense if needed
            X = X.toarray() if hasattr(X, "toarray") else X

            gene_idx = [var_names.get_loc(g) for g in genes]
            expr = X[:, gene_idx]

            # Expression fractions
            pct_in = (expr[cluster_cells.values] > 0).mean(axis=0)
            pct_out = (expr[other_cells.values] > 0).mean(axis=0)

            # Mean expression
            mean_in = expr[cluster_cells.values].mean(axis=0)
            mean_out = expr[other_cells.values].mean(axis=0)

            df["pct_in"] = pct_in
            df["pct_out"] = pct_out
            df["pct_diff"] = pct_in - pct_out

            df["mean_in"] = mean_in
            df["mean_out"] = mean_out

            # --------------------------------------------------
            # Cluster-level metadata
            # --------------------------------------------------
            cluster_size = cluster_cells.sum()
            df["cluster_size"] = cluster_size

            # --------------------------------------------------
            # Sorting (no filtering here)
            # --------------------------------------------------
            df = df.sort_values(
                ["pct_diff", "logfoldchange"],
                ascending=False
            )

            df.to_excel(writer, sheet_name=f"cluster_{cl}", index=False)

    # --------------------------------------------------
    # Save
    # --------------------------------------------------
    print(f"Saving updated AnnData → {output_file}")
    adata.write(output_file)
