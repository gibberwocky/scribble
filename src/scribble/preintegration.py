#!/usr/bin/env python

from pathlib import Path


def run_preintegration(args):
    import scanpy as sc
    import numpy as np
    import matplotlib.pyplot as plt
    import random
    from scribble.import_data import setup_environment

    PROJECT_DIR = Path(args.project_dir)
    PLOT_DIR = PROJECT_DIR / "sc_plots"
    setup_environment(sc, np, random, PLOT_DIR)

    input_file = Path(args.input)
    output_file = input_file.with_name(f"{input_file.stem}_preintegration.h5ad")
    umap_file = PLOT_DIR / f"{input_file.stem}_preintegration_umap.png"
    pca_file = PLOT_DIR / f"{input_file.stem}_preintegration_pca.png"

    print(f"Loading {input_file}")
    adata = sc.read(input_file)

    # Check requested args are present before proceeding
    missing = [v for v in args.vars if v not in adata.obs.columns]
    if missing:
        raise ValueError(f"Variables not found in adata.obs: {missing}")

    # Highly variable genes
    sc.pp.highly_variable_genes(
        adata,
        n_top_genes=args.hvgs, # 3000
        batch_key=args.batch, # "sample"
        flavor="seurat_v3"
    )

    # Normalisation
    sc.pp.normalize_total(adata)
    sc.pp.log1p(adata)
    adata.raw = adata.copy()
    adata = adata[:, adata.var.highly_variable].copy()

    print(f"Cells: {adata.n_obs}")
    print(f"HVGs: {adata.n_vars}")

    if not args.no_scale:
        sc.pp.scale(adata, max_value=10)
    else:
        print("Skipping scaling...")

    sc.tl.pca(adata, n_comps=args.npcs, mask_var="highly_variable")
    print("Top 10 PC variance ratios:")
    print(adata.uns["pca"]["variance_ratio"][:10])

    print("Generating PCA plot(s) ...")
    sc.pl.pca(
        adata,
        color=args.vars, # ["sample", "batchinfo", "sex"]
        wspace=0.4,
        show=False
    )
    plt.savefig(pca_file, dpi=300, bbox_inches="tight")
    plt.close()

    sc.pp.neighbors(adata, n_pcs=args.npcs, n_neighbors=args.neighbors)
    sc.tl.umap(adata)

    print("Generating UMAP plot(s) ...")
    sc.pl.umap(
        adata,
        color=args.vars, # ["sample", "batchinfo", "sex"]
        wspace=0.4,
        show=False
    )
    plt.savefig(umap_file, dpi=300, bbox_inches="tight")
    plt.close()

    # --------------------------------------------------
    # Save output
    # --------------------------------------------------
    print(f"Saving updated AnnData → {output_file}")
    adata.write(output_file)
