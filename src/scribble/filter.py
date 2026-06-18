#!/usr/bin/env python

from pathlib import Path


def run_filter(args):
    import scanpy as sc
    import pandas as pd
    import numpy as np
    import random
    from scribble.import_data import setup_environment

    PROJECT_DIR = Path(args.project_dir)
    PLOT_DIR = PROJECT_DIR / "sc_plots"
    setup_environment(sc, np, random, PLOT_DIR)

    input_file = Path(args.input)
    output_file = input_file.with_name(f"{input_file.stem}_filtered.h5ad")

    print(f"Loading {input_file}")
    adata = sc.read(input_file)

    # --------------------------------------------------
    # Load filters from Excel (optional)
    # --------------------------------------------------
    filter_dict = None

    if args.filter_xlsx:
        print(f"Loading filter definitions from {args.filter_xlsx}")

        filter_df = pd.read_excel(args.filter_xlsx, sheet_name="filters")

        required_cols = {"sample", "min_genes", "max_genes"}
        if not required_cols.issubset(filter_df.columns):
            raise ValueError(
                f"'filters' sheet must contain columns: {required_cols}"
            )

        filter_dict = (
            filter_df
            .set_index("sample")[["min_genes", "max_genes"]]
            .to_dict(orient="index")
        )

    # --------------------------------------------------
    # Build filtering masks
    # --------------------------------------------------
    print("Applying filters...")

    base_mask = (
        (~adata.obs["mt_outlier"]) &
        (~adata.obs["doublet_outlier"])
    )

    if filter_dict:
        gene_mask = np.zeros(adata.n_obs, dtype=bool)

        for sample, thresholds in filter_dict.items():
            idx = adata.obs["sample"] == sample

            gene_mask[idx] = (
                (adata.obs.loc[idx, "n_genes_by_counts"] > thresholds["min_genes"]) &
                (adata.obs.loc[idx, "n_genes_by_counts"] < thresholds["max_genes"])
            )

        # fallback for samples not listed
        unmatched = ~adata.obs["sample"].isin(filter_dict.keys())

        gene_mask[unmatched] = (
            (adata.obs.loc[unmatched, "n_genes_by_counts"] > args.mingenes) &
            (adata.obs.loc[unmatched, "n_genes_by_counts"] < args.maxgenes)
        )

    else:
        gene_mask = (
            (adata.obs["n_genes_by_counts"] > args.mingenes) &
            (adata.obs["n_genes_by_counts"] < args.maxgenes)
        )

    final_mask = base_mask & gene_mask

    # --------------------------------------------------
    # Apply filtering
    # --------------------------------------------------
    n_before = adata.n_obs
    adata = adata[final_mask].copy()
    n_after = adata.n_obs

    print(f"Retained {n_after} / {n_before} cells ({n_after/n_before:.1%})")

    print(f"Saving updated AnnData → {output_file}")
    adata.write(output_file)
