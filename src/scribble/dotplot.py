#!/usr/bin/env python

from pathlib import Path


def run_dotplot(args):

    import scanpy as sc
    import numpy as np
    import random
    import matplotlib.pyplot as plt

    from scribble.import_data import setup_environment
    from scribble.refine import restore_counts

    PROJECT_DIR = Path(args.project_dir)
    PLOT_DIR = (PROJECT_DIR / "scribble/plots")

    setup_environment(sc, np, random, PLOT_DIR)

    input_file = Path(args.input)
    output_plot = Path(args.output)

    print(f"Loading {input_file}")
    adata = sc.read(input_file)

    try:

        adata_plot = restore_counts(adata)
        adata_plot.obs = (adata.obs.copy())
        adata_plot.obsm = (adata.obsm.copy())
        adata = adata_plot
        sc.pp.normalize_total(adata)
        sc.pp.log1p(adata)
        adata.raw = adata

    except Exception as e:
        print(f"Could not restore counts: {e}")

    if args.samples:

        adata = adata[adata.obs["sample"].isin(args.samples)].copy()
        print(f"Cells after sample filter: {adata.n_obs:,}")

    if args.annotation not in adata.obs:
        raise ValueError(f"{args.annotation} not found in adata.obs")

    valid_genes = [
        gene
        for gene in args.genes
        if gene in adata.var_names
    ]

    missing_genes = sorted(
        set(args.genes)
        - set(valid_genes)
    )

    if missing_genes:

        print(
            f"Missing genes: "
            f"{missing_genes}"
        )

    if len(valid_genes) == 0:

        raise ValueError(
            "None of the requested "
            "genes are present."
        )

    print(
        f"Plotting "
        f"{len(valid_genes)} genes"
    )

    sc.pl.dotplot(
        adata,
        var_names=valid_genes,
        groupby=args.annotation,
        show=False,
    )

    plt.savefig(
        output_plot,
        dpi=300,
        bbox_inches="tight",
    )

    plt.close()
