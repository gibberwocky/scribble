#!/usr/bin/env python

from pathlib import Path





def run_dbl_qc(args):
    import scanpy as sc
    import numpy as np
    import matplotlib.pyplot as plt
    import random
    from scribble.import_data import setup_environment
    from scribble.plots import doublet_qc_panel

    PROJECT_DIR = Path(args.project_dir)
    PLOT_DIR = PROJECT_DIR / "sc_plots"
    setup_environment(sc, np, random, PLOT_DIR)

    input_file = Path(args.input)
    output_file = input_file.with_name(f"{input_file.stem}_dblqc_exp-{args.expected}.h5ad")
    plot_file = PLOT_DIR / f"{input_file.stem}_dblqc_exp-{args.expected}.png"

    print(f"Loading {input_file}")
    adata = sc.read(input_file)

    adata.obs["doublet_outlier"] = False
    doublet_stats = {}
    for sample in adata.obs["sample"].unique():
        mask = adata.obs["sample"] == sample
        scores = adata.obs.loc[mask, "doublet_score"]
        threshold = np.quantile(scores, 1-args.expected)
        outliers = scores > threshold
        adata.obs.loc[mask, "doublet_outlier"] = outliers
        n_cells = int(mask.sum())
        n_outliers = int(outliers.sum())
        frac = float(n_outliers / n_cells)
        doublet_stats[sample] = {
            "n_cells_pre": n_cells,
            "n_outliers": n_outliers,
            "fraction_outliers": frac,
            "threshold": float(threshold)
        }
        print(f"{sample}: {frac:.3f}")

    adata.uns["doublet_qc"] = doublet_stats

    # QC summary
    adata.uns["qc_summary"] = {
        "n_cells_total": int(adata.n_obs),
        "doublet_fraction_total": float(np.mean(adata.obs["doublet_outlier"])),
        "mt_fraction_total": float(np.mean(adata.obs["mt_outlier"]))
    }

    print("Generating doublet QC plot...")
    doublet_qc_panel(np, plt, adata, outfile=plot_file, threshold=1-args.expected)

    retention_stats = {}
    for sample in adata.obs["sample"].unique():
        mask = adata.obs["sample"] == sample
        n_total = int(mask.sum())
        keep_mask = (
            ~adata.obs["mt_outlier"] &
            ~adata.obs["doublet_outlier"]
        )
        mask_keep = mask & keep_mask
        n_keep = int(mask_keep.sum())
        frac_keep = float(n_keep / n_total)
        retention_stats[sample] = {
            "n_cells_pre": n_total,
            "n_cells_post": n_keep,
            "fraction_kept": frac_keep,
            "fraction_removed": 1 - frac_keep
        }

    adata.uns["qc_retention"] = retention_stats

    adata.uns["qc_retention_summary"] = {
        "n_cells_pre": int(adata.n_obs),
        "n_cells_post": int(sum(v["n_cells_post"] for v in retention_stats.values())),
        "fraction_kept": float(
            sum(v["n_cells_post"] for v in retention_stats.values()) / adata.n_obs
        )
    }

    print(f"Saving updated AnnData → {output_file}")
    adata.write(output_file, compression="gzip")
