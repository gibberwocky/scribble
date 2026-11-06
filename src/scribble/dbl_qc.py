#!/usr/bin/env python

from pathlib import Path


def run_dbl_qc(args):
    import os, sys
    os.environ['R_HOME'] = sys.exec_prefix+"/lib/R/"

    import scanpy as sc
    import numpy as np
    import pandas as pd
    import matplotlib.pyplot as plt
    import seaborn as sns
    import random
    import subprocess
    from pathlib import Path
    from scribble.import_data import setup_environment
    from scribble.utils import write_run_manifest
    from scribble.plots import doublet_qc_panel

    PROJECT_DIR = Path(args.project_dir)
    MANIFEST_DIR = PROJECT_DIR / "scribble/manifests"
    PLOT_DIR = PROJECT_DIR / "scribble/plots"
    TABLE_DIR = PROJECT_DIR / "scribble/tables"
    ADATA_DIR = PROJECT_DIR / "scribble/adata"

    setup_environment(sc, np, random, PLOT_DIR)

    input_file = Path(args.input)
    output_file = input_file.with_name(f"{input_file.stem}_dblqc.h5ad")
    plot_file = PLOT_DIR / f"{input_file.stem}_dblqc.png"
    xlsx_file = TABLE_DIR / f"{input_file.stem}_doublet_qc_summary.xlsx"
    manifest_file = (MANIFEST_DIR / f"{input_file.stem}_dbl_manifest.json")

    if not any([
        args.scrublet,
        args.scdblfinder,
    ]):
        raise ValueError("At least one doublet detection method must be selected.")

    # Import data
    print(f"Loading {input_file}")
    adata = sc.read(input_file)

    # initialise
    if args.scrublet:
        adata.obs["scrublet_doublet"] = False
        adata.obs["scrublet_score"] = np.nan

    if args.scdblfinder:
        adata.obs["scdblfinder_doublet"] = False

    adata.obs["doublet_outlier"] = False

    doublet_stats = {}
    fig_paths = []
    for sample in adata.obs["sample"].unique():

        mask = adata.obs["sample"] == sample
        n_cells = int(mask.sum())

        # ✅ use CLI min_cells
        if n_cells < args.min_cells:
            print(f"WARNING: {sample} has very few cells ({n_cells}); skipping doublet detection")

            # mark all as non-doublets, record stats
            stats = {
                "n_cells_pre": n_cells,
                "n_outliers": 0,
                "fraction_outliers": 0.0,
                "skipped": True,
            }

            if args.scrublet:
                stats["expected_rate"] = args.scrublet_expected
                stats["n_scrublet"] = 0
                stats["scrublet_threshold"] = np.nan

            if args.scdblfinder:
                stats["n_scdblfinder"] = 0

            doublet_stats[sample] = stats

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

        scrublet_pred = None
        scores = None
        scrub_thresh = np.nan

        # -------- run Scrublet --------
        if args.scrublet:
            sc.external.pp.scrublet(ad_sample, expected_doublet_rate=args.scrublet_expected, threshold=None)
            scores = ad_sample.obs["doublet_score"].values
            scrublet_pred = ad_sample.obs["predicted_doublet"].values
            adata.obs.loc[mask, "scrublet_score"] = scores
            adata.obs.loc[mask, "scrublet_doublet"] = scrublet_pred
            scrub_thresh = float(ad_sample.uns.get("scrublet", {}).get("threshold", np.nan))

        # -------- run scDblFinder --------
        if args.scdblfinder:
            scdblfinder_pred = np.zeros(ad_sample.n_obs, dtype=bool,)
            tmp_h5ad = (ADATA_DIR / f"{sample}_scdblfinder.h5ad")
            tmp_csv = (ADATA_DIR / f"{sample}_scdblfinder.csv")

            # save temporary AnnData
            ad_sample.write(tmp_h5ad, compression="gzip",)

            r_script = f"""
            suppressPackageStartupMessages(library(scDblFinder))
            suppressPackageStartupMessages(library(anndataR))
            suppressPackageStartupMessages(library(SingleCellExperiment))

            ad <- read_h5ad("{tmp_h5ad}")
            class(ad)
            sce <- ad$as_SingleCellExperiment()
            class(sce)
            sce <- scDblFinder(sce)

            write.csv(
                data.frame(
                    cell = colnames(sce),
                    class = sce$scDblFinder.class
                ),
                "{tmp_csv}",
                row.names = FALSE
            )
            """

            result = subprocess.run(
                ["Rscript", "-e", r_script],
                capture_output=True,
                text=True,
            )

            if result.returncode != 0:
                print(result.stdout)
                print(result.stderr)
                raise RuntimeError(f"scDblFinder failed for {sample}")

            sdf = pd.read_csv(tmp_csv)
            scdblfinder_pred = (sdf["class"].eq("doublet").values)
            adata.obs.loc[mask, "scdblfinder_doublet"] = scdblfinder_pred

            # cleanup
            tmp_h5ad.unlink(missing_ok=True)
            tmp_csv.unlink(missing_ok=True)


        outliers = np.zeros(ad_sample.n_obs, dtype=bool,)

        if args.scrublet:
            outliers |= scrublet_pred

        if args.scdblfinder:
            outliers |= scdblfinder_pred

        adata.obs.loc[mask, "doublet_outlier"] = outliers

        n_outliers = int(outliers.sum())
        frac = float(n_outliers / n_cells)

        stats = {
            "n_cells_pre": n_cells,
            "n_outliers": n_outliers,
            "fraction_outliers": frac,
            "skipped": False,
        }

        if args.scrublet:
            stats["expected_rate"] = args.scrublet_expected
            stats["n_scrublet"] = int(scrublet_pred.sum())
            stats["scrublet_threshold"] = (scrub_thresh)

        if args.scdblfinder:
            stats["n_scdblfinder"] = int(scdblfinder_pred.sum())

        doublet_stats[sample] = stats


        print(f"{sample}: {frac:.3f}")

        # -------- diagnostic plot --------
        fig, ax = plt.subplots(figsize=(5, 4))

        if (
            args.scrublet
            and scores is not None
            and np.unique(scores).size > 1
        ):
            sns.histplot(scores, bins=50, ax=ax, color="steelblue")
        elif args.scrublet:
            ax.text(0.5, 0.5, "Constant scores", ha="center")
        else:
            ax.text(0.5, 0.5, "No Scrublet scores", ha="center")


        if (
            args.scrublet
            and np.isfinite(scrub_thresh)
        ):
            ax.axvline(scrub_thresh, color="green", linestyle=":", label="Scrublet")

        ax.set_title(f"{sample} doublet score")
        ax.set_xlabel("Score")
        ax.set_ylabel("Cells")

        handles, labels = ax.get_legend_handles_labels()

        if len(handles) > 0:
            ax.legend()

        fig_path = PLOT_DIR / f"{input_file.stem}_{sample}_doublet_hist.png"
        fig_paths.append(fig_path)
        plt.savefig(fig_path, dpi=150, bbox_inches="tight")
        plt.close()

    adata.uns["doublet_methods"] = {
        "scrublet": args.scrublet,
        "scdblfinder": args.scdblfinder,
    }
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
        threshold=1 - args.scrublet_expected
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

        row = {
            "sample": sample,
            "n_cells_pre": dbl.get("n_cells_pre"),
            "n_doublets": dbl.get("n_outliers"),
            "doublet_fraction": dbl.get("fraction_outliers"),
            "n_cells_post": ret.get("n_cells_post"),
            "fraction_kept": ret.get("fraction_kept"),
            "fraction_removed": ret.get("fraction_removed"),
        }

        if args.scrublet:
            row["scrublet_expected_doublet_rate"] = dbl.get("expected_rate")
            row["n_scrublet"] = (dbl.get("n_scrublet"))
            row["scrublet_threshold"] = (dbl.get("scrublet_threshold"))

        if args.scdblfinder:
            row["n_scdblfinder"] = (dbl.get("n_scdblfinder"))

        rows.append(row)

    qc_df = pd.DataFrame(rows)

    # -------- Save as Excel file --------
    qc_df.to_excel(xlsx_file, index=False)
    print(f"Saved QC summary → {xlsx_file}")

    adata.uns["scribble_manifest"] = str(manifest_file)
    print(f"Saving updated AnnData → {output_file}")
    adata.write(output_file, compression="gzip")

    print(f"Writing run manifest to {manifest_file}")
    write_run_manifest(
        manifest_file=manifest_file,
        module="dbl",
        args=args,
        inputs={
            "input_file": str(input_file),
        },
        outputs={
            "adata": str(output_file),
            "plot(s)": [str(fig_path) for fig_path in fig_paths],
            "doublet qc plot": str(plot_file),
            "summary xlsx": str(xlsx_file)
        },
    )
