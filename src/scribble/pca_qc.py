#!/usr/bin/env python

from pathlib import Path


def run_pca_qc(args):
    import scanpy as sc
    import numpy as np
    import matplotlib.pyplot as plt
    import random
    from scribble.import_data import setup_environment
    from scribble.plots import pca_qc_panel

    PROJECT_DIR = Path(args.project_dir)
    PLOT_DIR = PROJECT_DIR / "plots"
    setup_environment(sc, np, random, PLOT_DIR)

    input_file = Path(args.input)
    plot_file = PLOT_DIR / f"{input_file.stem}_pca.png"

    print(f"Loading {input_file}")
    adata = sc.read(input_file)

    print("Generating PCA before/after filtering plot...")
    pca_qc_panel(np, sc, plt, adata, plot_file,
        args.mingenes, args.maxgenes, args.hvgs, args.vmax)
