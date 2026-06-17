#!/usr/bin/env python

import argparse
from pathlib import Path
import scanpy as sc
import scvelo as scv
import numpy as np
import pandas as pd
import random
import scipy.sparse as sp
import matplotlib.pyplot as plt


# ---------- Environment setup ----------
def setup_environment(plot_dir):
    sc.settings.verbosity = 3
    scv.settings.verbosity = 3
    sc.logging.print_versions()
    sc.settings.figdir = plot_dir
    sc.set_figure_params(dpi=100, facecolor="white")

    np.random.seed(0)
    random.seed(0)


# ---------- Sample detection ----------
def get_samples(velo_dir):
    samples = sorted([d.name for d in velo_dir.iterdir() if d.is_dir()])
    if len(samples) == 0:
        raise RuntimeError("No sample folders found in velocyto directory")
    return samples

# ---------- Sample metadata -----------
def add_metadata(adata, metadata_file):
    print(f"Reading metadata: {metadata_file}")

    meta = pd.read_excel(metadata_file)

    # ---------- validation ----------
    if "sample" not in meta.columns:
        raise ValueError("Metadata file must contain a 'sample' column")

    meta = meta.set_index("sample")

    # check overlap
    adata_samples = set(adata.obs["sample"].unique())
    meta_samples = set(meta.index)

    missing = adata_samples - meta_samples
    if len(missing) > 0:
        raise ValueError(f"Metadata missing samples: {missing}")

    extra = meta_samples - adata_samples
    if extra:
        print(f"Warning: metadata contains unused samples: {extra}")

    # ---------- join ----------
    adata.obs = adata.obs.join(meta, on="sample")

    # ---------- type inference ----------
    for col in meta.columns:
        # skip sample column (already handled)
        if col == "sample":
            continue

        series = adata.obs[col]

        # try numeric conversion
        numeric = pd.to_numeric(series, errors="coerce")

        # heuristic: if most values convert → treat as numeric
        frac_numeric = numeric.notna().mean()

        if frac_numeric > 0.9:
            adata.obs[col] = numeric
            print(f"{col}: numeric")
        else:
            adata.obs[col] = series.astype("string").astype("category")
            print(f"{col}: categorical")

    # ensure sample is categorical
    adata.obs["sample"] = adata.obs["sample"].astype("category")

    return adata


# ---------- Barcode cleaning ----------
def clean_barcodes(idx, sample):
    return (
        idx
        .str.replace(r"^.*:", "", regex=True)
        .str.replace(r"x$", "", regex=True)
        .str.replace(r"-1$", "", regex=True)
        + "_" + sample
    )


# ---------- Process one sample ----------
def process_sample(sample, cellranger_dir, velo_dir):
    print(f"\nProcessing {sample}")

    sample_dir = velo_dir / sample

    loom_files = list(sample_dir.glob("*.loom"))
    if len(loom_files) != 1:
        raise RuntimeError(f"{sample}: Expected 1 loom file, found {len(loom_files)}")

    loom_file = loom_files[0]
    print(f"Using loom file: {loom_file}")

    # Load loom
    adata = sc.read(loom_file, cache=True)
    adata.obs_names = clean_barcodes(adata.obs_names, sample)

    # Load 10X
    adata_10X = sc.read_10x_mtx(
        cellranger_dir / sample / "outs/filtered_feature_bc_matrix"
    )
    adata_10X.obs_names = clean_barcodes(adata_10X.obs_names, sample)

    # QC metrics
    adata_10X.var["mt"] = adata_10X.var_names.str.contains("^MT-|^mt-")
    sc.pp.calculate_qc_metrics(adata_10X, qc_vars=["mt"], inplace=True)

    # Consistency check
    missing = set(adata.obs_names) - set(adata_10X.obs_names)
    if len(missing) > 0:
        raise ValueError(f"{sample}: {len(missing)} barcodes missing from 10X data")

    # Fix var names
    adata.var_names = adata.var.index
    adata.var_names_make_unique()

    # Extract layers
    spliced = sp.csr_matrix(adata.layers["spliced"], dtype=np.float32)
    unspliced = sp.csr_matrix(adata.layers["unspliced"], dtype=np.float32)

    counts = spliced + unspliced

    # Scrublet
    adata_counts = sc.AnnData(X=counts.copy())
    adata_counts.obs_names = adata.obs_names.copy()
    adata_counts.var_names = adata.var_names.copy()

    print(f"Running Scrublet for {sample}")
    sc.external.pp.scrublet(adata_counts, threshold=0.25)

    # Build clean object
    adata_clean = sc.AnnData(
        X=counts,
        obs=adata.obs.copy(),
        var=adata.var.copy()
    )

    adata_clean.layers["spliced"] = spliced
    adata_clean.layers["unspliced"] = unspliced

    # Transfer QC metrics
    adata_clean.obs["doublet_score"] = adata_counts.obs["doublet_score"].values
    adata_clean.obs["predicted_doublet"] = adata_counts.obs["predicted_doublet"].values

    meta_10X = adata_10X.obs.loc[adata_clean.obs_names]

    adata_clean.obs["pct_counts_mt"] = meta_10X["pct_counts_mt"].values
    adata_clean.obs["total_counts"] = meta_10X["total_counts"].values
    adata_clean.obs["n_genes_by_counts"] = meta_10X["n_genes_by_counts"].values

    return adata_clean

# ---------- Main ----------
def run_import(args):
    from pathlib import Path

    PROJECT_DIR = Path(args.project_dir)
    metadata_file = args.metadata_file

    CELLRANGER_DIR = PROJECT_DIR / "cellranger"
    VELO_DIR = PROJECT_DIR / "velocyto"
    PLOT_DIR = PROJECT_DIR / "sc_plots"
    ADATA_DIR = PROJECT_DIR / "adata"

    PLOT_DIR.mkdir(exist_ok=True)
    ADATA_DIR.mkdir(exist_ok=True)

    setup_environment(PLOT_DIR)

    samples = get_samples(VELO_DIR)
    print(f"Detected samples: {samples}")

    adatas = []

    for sample in samples:
        adata_clean = process_sample(sample, CELLRANGER_DIR, VELO_DIR)
        adata_clean.obs["sample"] = sample

        adatas.append(adata_clean)

    print("\nConcatenating samples...")
    adata = sc.concat(adatas, label="sample", keys=samples)

    # ---------- ADD METADATA ----------
    print("\nAppending metadata...")
    adata = add_metadata(adata, args.metadata_file)

    out_file = ADATA_DIR / "raw.h5ad"
    print(f"Saving: {out_file}")

    adata.write(out_file)

