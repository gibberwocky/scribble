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
    counts_file = input_file.with_name(f"{input_file.stem}_filtered.tsv")

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
    # Capture BEFORE counts
    # --------------------------------------------------
    print("\nComputing pre-filter counts...")
    sample_counts_before = adata.obs["sample"].value_counts()

    # --------------------------------------------------
    # Build individual masks
    # --------------------------------------------------
    print("Building filter masks...")

    mt_mask = ~adata.obs["mt_outlier"]
    dbl_mask = ~adata.obs["doublet_outlier"]

    # gene mask (sample-aware)
    if filter_dict:
        gene_mask = np.zeros(adata.n_obs, dtype=bool)

        for sample, thresholds in filter_dict.items():
            idx = adata.obs["sample"] == sample

            gene_mask[idx] = (
                (adata.obs.loc[idx, "n_genes_by_counts"] > thresholds["min_genes"]) &
                (adata.obs.loc[idx, "n_genes_by_counts"] < thresholds["max_genes"])
            )

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

    # --------------------------------------------------
    # Combine masks
    # --------------------------------------------------
    final_mask = mt_mask & dbl_mask & gene_mask

    # --------------------------------------------------
    # Apply filtering
    # --------------------------------------------------
    n_before = adata.n_obs
    adata_filtered = adata[final_mask].copy()
    n_after = adata_filtered.n_obs

    print(f"\nRetained {n_after} / {n_before} cells ({n_after/n_before:.1%})")

    sample_counts_after = adata_filtered.obs["sample"].value_counts()

    # --------------------------------------------------
    # Per-sample logging + TSV export
    # --------------------------------------------------
    print("\nPer-sample QC summary:")
    print("-" * 50)

    rows = []

    for sample in sample_counts_before.index:
        idx = adata.obs["sample"] == sample

        before = sample_counts_before[sample]
        after = sample_counts_after.get(sample, 0)

        mt_fail = (~mt_mask & idx).sum()
        dbl_fail = (~dbl_mask & idx).sum()
        gene_fail = (~gene_mask & idx).sum()

        removed = before - after
        retained_frac = after / before if before > 0 else 0.0

        # Terminal output
        print(f"{sample}")
        print(f"  kept:    {after}/{before} ({retained_frac:.1%})")
        print(f"  removed: {removed} ({removed/before:.1%})")
        print(f"    - mt_outlier:      {mt_fail}")
        print(f"    - doublet_outlier: {dbl_fail}")
        print(f"    - gene_filter:     {gene_fail}")
        print()

        # Collect row for TSV
        rows.append({
            "sample": sample,
            "before": before,
            "after": after,
            "retained_frac": retained_frac,
            "removed": removed,
            "mt_outlier": mt_fail,
            "doublet_outlier": dbl_fail,
            "gene_filter": gene_fail,
        })

    rows.append({
        "sample": "GLOBAL",
        "before": n_before,
        "after": n_after,
        "retained_frac": n_after / n_before,
        "removed": n_before - n_after,
        "mt_outlier": (~mt_mask).sum(),
        "doublet_outlier": (~dbl_mask).sum(),
        "gene_filter": (~gene_mask).sum(),
    })

    # --------------------------------------------------
    # Write TSV file
    # --------------------------------------------------
    df = pd.DataFrame(rows)

    print(f"Writing QC summary → {counts_file}")
    df.to_csv(counts_file, sep="\t", index=False)

    # --------------------------------------------------
    # Save output
    # --------------------------------------------------
    print(f"Saving updated AnnData → {output_file}")
    adata_filtered.write(output_file)
