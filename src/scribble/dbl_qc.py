#!/usr/bin/env python

from pathlib import Path


def run_dbl_qc(args):
    import scanpy as sc
    import numpy as np
    import pandas as pd
    import matplotlib.pyplot as plt
    import seaborn as sns
    import random
    from pathlib import Path
    from scribble.import_data import setup_environment
    from scribble.plots import doublet_qc_panel

    PROJECT_DIR = Path(args.project_dir)
    PLOT_DIR = PROJECT_DIR / "scribble/plots"
    TABLE_DIR = PROJECT_DIR / "scribble/tables"

    setup_environment(sc, np, random, PLOT_DIR)

    input_file = Path(args.input)
    output_file = input_file.with_name(f"{input_file.stem}_dblqc_exp-{args.expected}.h5ad")
    plot_file = PLOT_DIR / f"{input_file.stem}_dblqc_exp-{args.expected}.png"

    print(f"Loading {input_file}")
    adata = sc.read(input_file)

    # initialise
    adata.obs["doublet_outlier"] = False
    adata.obs["doublet_score"] = np.nan
    doublet_stats = {}

    print("Running Scrublet per sample...")

    for sample in adata.obs["sample"].unique():

        mask = adata.obs["sample"] == sample
        n_cells = int(mask.sum())

        # ✅ use CLI min_cells
        if n_cells < args.min_cells:
            print(f"WARNING: {sample} has very few cells ({n_cells}); skipping Scrublet")

            # mark all as non-doublets, record stats
            doublet_stats[sample] = {
                "n_cells_pre": n_cells,
                "n_outliers": 0,
                "fraction_outliers": 0.0,
                "quantile_threshold": np.nan,
                "scrublet_threshold": np.nan,
                "expected_rate": args.expected,
                "skipped": True
            }
            continue

        # -------- build per-sample AnnData --------
        mask_array = mask.values
        counts = adata.layers["counts"][mask_array]

        ad_sample = sc.AnnData(X=counts.copy())
        ad_sample.obs_names = adata.obs_names[mask_array]
        ad_sample.var_names = adata.var_names

        # -------- normalisation for Scrublet (optional) --------
        if getattr(args, "normalize", False):
            sc.pp.normalize_total(ad_sample, target_sum=1e4)
            sc.pp.log1p(ad_sample)

        # -------- run Scrublet --------
        sc.external.pp.scrublet(
            ad_sample,
            expected_doublet_rate=args.expected,
            threshold=None
        )

        scores = ad_sample.obs["doublet_score"].values
        scrublet_pred = ad_sample.obs["predicted_doublet"].values

        # store scores
        adata.obs.loc[mask, "doublet_score"] = scores

        # -------- quantile threshold --------
        if len(scores) < 10:
            quantile_threshold = np.nan
            quantile_mask = np.zeros_like(scores, dtype=bool)
        else:
            quantile_threshold = np.quantile(scores, 1 - args.expected)
            quantile_mask = scores > quantile_threshold

        # -------- thresholding mode --------
        mode = getattr(args, "mode", "hybrid")

        if mode == "hybrid":
            outliers = scrublet_pred | quantile_mask
        elif mode == "scrublet":
            outliers = scrublet_pred
        elif mode == "quantile":
            outliers = quantile_mask
        else:
            raise ValueError(f"Unknown mode: {mode}")

        adata.obs.loc[mask, "doublet_outlier"] = outliers

        n_outliers = int(outliers.sum())
        frac = float(n_outliers / n_cells)

        doublet_stats[sample] = {
            "n_cells_pre": n_cells,
            "n_outliers": n_outliers,
            "fraction_outliers": frac,
            "quantile_threshold": float(quantile_threshold) if not np.isnan(quantile_threshold) else np.nan,
            "scrublet_threshold": float(
                ad_sample.uns.get("scrublet", {}).get("threshold", np.nan)
            ),
            "expected_rate": args.expected,
            "skipped": False
        }

        print(f"{sample}: {frac:.3f}")

        # -------- diagnostic plot --------
        fig, ax = plt.subplots(figsize=(5, 4))

        if np.unique(scores).size > 1:
            sns.histplot(scores, bins=50, ax=ax, color="steelblue")
        else:
            ax.text(0.5, 0.5, "Constant scores", ha="center")

        if not np.isnan(quantile_threshold):
            ax.axvline(quantile_threshold, color="red", linestyle="--", label="Quantile")

        scrub_thresh = ad_sample.uns.get("scrublet", {}).get("threshold", np.nan)
        if scrub_thresh is not None:
            ax.axvline(scrub_thresh, color="green", linestyle=":", label="Scrublet")

        ax.set_title(f"{sample} doublet score")
        ax.set_xlabel("Score")
        ax.set_ylabel("Cells")
        ax.legend()

        fig_path = PLOT_DIR / f"{input_file.stem}_{sample}_doublet_hist.png"
        plt.savefig(fig_path, dpi=150, bbox_inches="tight")
        plt.close()

    adata.uns["doublet_qc"] = doublet_stats

    # -------- QC summary --------
    adata.uns["qc_summary"] = {
        "n_cells_total": int(adata.n_obs),
        "doublet_fraction_total": float(np.mean(adata.obs["doublet_outlier"])),
        "mt_fraction_total": float(np.mean(adata.obs["mt_outlier"]))
    }

    print("Generating global doublet QC plot...")
    doublet_qc_panel(
        np,
        plt,
        sns,
        adata,
        outfile=plot_file,
        threshold=1 - args.expected
    )

    # -------- retention stats --------
    retention_stats = {}

    for sample in adata.obs["sample"].unique():
        mask = adata.obs["sample"] == sample
        n_total = int(mask.sum())

        if n_total == 0:
            continue

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

    total_post = sum(v["n_cells_post"] for v in retention_stats.values())

    adata.uns["qc_retention_summary"] = {
        "n_cells_pre": int(adata.n_obs),
        "n_cells_post": int(total_post),
        "fraction_kept": float(total_post / adata.n_obs) if adata.n_obs > 0 else 0.0
    }

    # Generate QC summary
    rows = []

    for sample in adata.obs["sample"].unique():

        dbl = adata.uns["doublet_qc"].get(sample, {})
        ret = adata.uns["qc_retention"].get(sample, {})

        rows.append({
            "sample": sample,

            # doublet stats
            "n_cells_pre": dbl.get("n_cells_pre"),
            "n_doublets": dbl.get("n_outliers"),
            "doublet_fraction": dbl.get("fraction_outliers"),
            "expected_doublet_rate": dbl.get("expected_rate"),
            "quantile_threshold": dbl.get("quantile_threshold"),
            "scrublet_threshold": dbl.get("scrublet_threshold"),
            "scrublet_skipped": dbl.get("skipped", False),

            # retention
            "n_cells_post": ret.get("n_cells_post"),
            "fraction_kept": ret.get("fraction_kept"),
            "fraction_removed": ret.get("fraction_removed"),
        })

    qc_df = pd.DataFrame(rows)

    # -------- save as CSV --------
    xlsx_file = TABLE_DIR / f"{input_file.stem}_doublet_qc_summary.xlsx"
    qc_df.to_excel(xlsx_file, index=False)
    print(f"Saved QC summary → {xlsx_file}")

    print(f"Saving updated AnnData → {output_file}")
    adata.write(output_file, compression="gzip")
