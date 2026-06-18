#!/usr/bin/env python

from pathlib import Path


def run_cluster(args):
    import scanpy as sc
    import numpy as np
    import pandas as pd
    import matplotlib.pyplot as plt
    import random
    from pathlib import Path
    from scribble.import_data import setup_environment

    PROJECT_DIR = Path(args.project_dir)
    PLOT_DIR = PROJECT_DIR / "sc_plots"
    setup_environment(sc, np, random, PLOT_DIR)

    input_file = Path(args.input)
    output_file = input_file.with_name(f"{input_file.stem}_clustered.h5ad")

    print(f"Loading {input_file}")
    adata = sc.read(input_file)

    # --------------------------------------------------
    # Validation
    # --------------------------------------------------
    if args.embedding not in adata.obsm:
        raise ValueError(f"{args.embedding} not found in adata.obsm")

    # --------------------------------------------------
    # Build neighbors graph
    # --------------------------------------------------
    print(f"Building neighbors (embedding={args.embedding})")

    sc.pp.neighbors(
        adata,
        use_rep=args.embedding,
        n_neighbors=args.neighbors
    )

    # --------------------------------------------------
    # Run Leiden clustering
    # --------------------------------------------------
    print(f"Running Leiden clustering (resolution={args.resolution})")

    sc.tl.leiden(
        adata,
        resolution=args.resolution,
        key_added="leiden"
    )

    print(f"Number of clusters: {adata.obs['leiden'].nunique()}")

    # --------------------------------------------------
    # OPTIONAL: cluster stability (repeat runs)
    # --------------------------------------------------
    if args.n_repeats > 1:
        print(f"Computing cluster stability ({args.n_repeats} repeats)...")

        cluster_matrix = np.zeros((adata.n_obs, args.n_repeats))

        for i in range(args.n_repeats):
            sc.tl.leiden(
                adata,
                resolution=args.resolution,
                key_added=f"leiden_tmp_{i}"
            )
            cluster_matrix[:, i] = adata.obs[f"leiden_tmp_{i}"].astype(int)

        # Stability = most frequent assignment frequency
        stability = []

        for row in cluster_matrix:
            counts = np.bincount(row.astype(int))
            stability.append(counts.max() / args.n_repeats)

        adata.obs["cluster_stability"] = stability

        print(f"Mean stability: {np.mean(stability):.3f}")

    # --------------------------------------------------
    # Export cluster stats (if repeats > 1)
    # --------------------------------------------------
    if args.n_repeats > 1:

        stats_file = input_file.with_name(f"{input_file.stem}_cluster_stats.tsv")
        cell_file = input_file.with_name(f"{input_file.stem}_cluster_cells.tsv")

        print(f"Exporting cluster stats → {stats_file}")
        print(f"Exporting per-cell data → {cell_file}")

        # --------------------------------------------------
        # Per-cell data
        # --------------------------------------------------
        cell_df = adata.obs[["leiden", "cluster_stability"]].copy()
        cell_df.to_csv(cell_file, sep="\t")

        # --------------------------------------------------
        # Per-cluster summary
        # --------------------------------------------------
        cluster_summary = (
            adata.obs
            .groupby("leiden")
            .agg(
                n_cells=("leiden", "size"),
                mean_stability=("cluster_stability", "mean"),
                median_stability=("cluster_stability", "median")
            )
        )

        # Add sample composition
        sample_counts = (
            adata.obs
            .groupby(["leiden", "sample"])
            .size()
            .unstack(fill_value=0)
        )

        cluster_summary = cluster_summary.join(sample_counts)

        # Add fraction of total cells
        cluster_summary["fraction"] = cluster_summary["n_cells"] / adata.n_obs

        # Add parameters
        cluster_summary["resolution"] = args.resolution
        cluster_summary["embedding"] = args.embedding
        cluster_summary["n_repeats"] = args.n_repeats

        # Sort by size
        cluster_summary = cluster_summary.sort_values("n_cells", ascending=False)

        cluster_summary.to_csv(stats_file, sep="\t")


    # --------------------------------------------------
    # Ensure UMAP exists
    # --------------------------------------------------
    if "X_umap" not in adata.obsm:
        print("UMAP not found — computing UMAP...")
        sc.tl.umap(adata)

    # --------------------------------------------------
    # Plotting
    # --------------------------------------------------
    print("Generating clustering plots...")

    # Cluster plot
    cluster_file = PLOT_DIR / f"{input_file.stem}_clusters.png"
    sc.pl.umap(
        adata,
        color=["leiden"] + args.vars,
        wspace=0.4,
        show=False
    )
    plt.savefig(cluster_file, dpi=300, bbox_inches="tight")
    plt.close()

    # Stability plot
    if "cluster_stability" in adata.obs:
        stability_file = PLOT_DIR / f"{input_file.stem}_cluster_stability.png"
        sc.pl.umap(
            adata,
            color="cluster_stability",
            cmap="viridis",
            show=False
        )
        plt.savefig(stability_file, dpi=300, bbox_inches="tight")
        plt.close()

    # --------------------------------------------------
    # Save
    # --------------------------------------------------
    print(f"Saving updated AnnData → {output_file}")
    adata.write(output_file)
