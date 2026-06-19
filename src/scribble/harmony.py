#!/usr/bin/env python

from pathlib import Path


def run_harmony(args):
    import scanpy as sc
    import numpy as np
    import matplotlib.pyplot as plt
    import harmonypy as hm
    import random
    from scribble.import_data import setup_environment

    PROJECT_DIR = Path(args.project_dir)
    PLOT_DIR = PROJECT_DIR / "sc_plots"
    setup_environment(sc, np, random, PLOT_DIR)

    input_file = Path(args.input)
    output_file = input_file.with_name(f"{input_file.stem}_harmony_theta-{args.theta}.h5ad")
    umap_file = PLOT_DIR / f"{input_file.stem}_harmony_theta-{args.theta}_umap.png"

    print(f"Loading {input_file}")
    adata = sc.read(input_file)

    if "X_pca" not in adata.obsm:
        raise ValueError("X_pca not found. Run preintegration first.")

    # --------------------------------------------------
    # Harmony
    # --------------------------------------------------
    ho = hm.run_harmony(
        adata.obsm['X_pca'],
        adata.obs,
        args.batch,
        theta=args.theta
    )
    adata.obsm['X_pca_harmony'] = ho.Z_corr

    sc.pp.neighbors(adata, use_rep="X_pca_harmony", n_pcs=args.npcs, n_neighbors=args.neighbors)
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
