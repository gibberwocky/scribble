#!/usr/bin/env python

from pathlib import Path


def compute_mt_outliers(np, adata, nmads):

    adata.obs["mt_outlier"] = False
    mt_stats = {}

    for sample in adata.obs["sample"].unique():
        mask = adata.obs["sample"] == sample
        mt = adata.obs.loc[mask, "pct_counts_mt"]

        median = np.median(mt)
        mad = np.median(np.abs(mt - median))

        threshold = median + nmads * mad
        outliers = (mt - median) > nmads * mad

        adata.obs.loc[mask, "mt_outlier"] = outliers

        n_cells = int(mask.sum())
        n_outliers = int(outliers.sum())
        frac = float(n_outliers / n_cells)

        mt_stats[sample] = {
            "n_cells_pre": n_cells,
            "n_outliers": n_outliers,
            "fraction_outliers": frac,
            "threshold": float(threshold),
            "median": float(median),
            "mad": float(mad),
            "nmads": nmads,
        }

        print(f"{sample}: {frac:.3f}")

    adata.uns["mt_qc"] = mt_stats

    return adata




def run_mt_qc(args):
    import scanpy as sc
    import numpy as np
    import matplotlib.pyplot as plt
    import random
    from scribble.import_data import setup_environment
    from scribble.plots import mt_qc_panel

    PROJECT_DIR = Path(args.project_dir)
    PLOT_DIR = PROJECT_DIR / "sc_plots"
    setup_environment(sc, np, random, PLOT_DIR)

    input_file = Path(args.input)
    output_file = input_file.with_name(f"{input_file.stem}_mtqc_nMADs-{args.nmads}.h5ad")
    plot_file = PLOT_DIR / f"_mtqc_nMADs-{args.nmads}.png"

    print(f"Loading {input_file}")
    adata = sc.read(input_file)

    adata = compute_mt_outliers(np, adata, args.nmads)

    print("Generating MT QC plot...")
    mt_qc_panel(np, plt, adata, plot_file, args.nmads)

    print(f"Saving updated AnnData → {output_file}")
    adata.write(output_file)
