#!/usr/bin/env python
from pathlib import Path


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


# ---------- Metadata ----------
def add_metadata(adata, metadata_file):
    print(f"Reading metadata: {metadata_file}")

    meta = pd.read_excel(metadata_file)

    if "sample" not in meta.columns:
        raise ValueError("Metadata must contain a 'sample' column")

    meta = meta.set_index("sample")

    adata.obs = adata.obs.join(meta, on="sample")

    # type inference
    for col in meta.columns:
        series = adata.obs[col]
        numeric = pd.to_numeric(series, errors="coerce")
        if numeric.notna().mean() > 0.9:
            adata.obs[col] = numeric
        else:
            adata.obs[col] = series.astype("string").astype("category")

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
    print(f"Processing {sample}")

    loom_file = list((velo_dir / sample).glob("*.loom"))[0]
    adata = sc.read(loom_file)
    adata.obs_names = clean_barcodes(adata.obs_names, sample)

    adata_10X = sc.read_10x_mtx(
        cellranger_dir / sample / "outs/filtered_feature_bc_matrix"
    )
    adata_10X.obs_names = clean_barcodes(adata_10X.obs_names, sample)

    # MT + QC
    adata_10X.var["mt"] = adata_10X.var_names.str.contains("^MT-|^mt-")
    sc.pp.calculate_qc_metrics(adata_10X, qc_vars=["mt"], inplace=True)

    # Layers
    spliced = sp.csr_matrix(adata.layers["spliced"], dtype=np.float32)
    unspliced = sp.csr_matrix(adata.layers["unspliced"], dtype=np.float32)
    counts = spliced + unspliced

    # Scrublet
    adata_counts = sc.AnnData(X=counts.copy())
    adata_counts.obs_names = adata.obs_names
    adata_counts.var_names = adata.var_names

    sc.external.pp.scrublet(adata_counts, threshold=0.25)

    # Build final object
    adata_clean = sc.AnnData(
        X=counts,
        obs=adata.obs.copy(),
        var=adata.var.copy()
    )

    adata_clean.layers["spliced"] = spliced
    adata_clean.layers["unspliced"] = unspliced

    adata_clean.obs["doublet_score"] = adata_counts.obs["doublet_score"]
    adata_clean.obs["predicted_doublet"] = adata_counts.obs["predicted_doublet"]

    meta_10X = adata_10X.obs.loc[adata_clean.obs_names]

    adata_clean.obs["pct_counts_mt"] = meta_10X["pct_counts_mt"]
    adata_clean.obs["total_counts"] = meta_10X["total_counts"]
    adata_clean.obs["n_genes_by_counts"] = meta_10X["n_genes_by_counts"]

    return adata_clean


# ---------- MAIN ENTRY ----------
def run_import(args):
    import scanpy as sc
    import numpy as np
    import pandas as pd
    import random
    import scipy.sparse as sp

    PROJECT_DIR = Path(args.project_dir)

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
        ad = process_sample(sample, CELLRANGER_DIR, VELO_DIR)
        ad.obs["sample"] = sample
        adatas.append(ad)

    print("Concatenating...")
    adata = sc.concat(adatas, label="sample", keys=samples)

    adata = add_metadata(adata, args.metadata_file)

    out_file = ADATA_DIR / "combined.h5ad"
    print(f"Saving to {out_file}")

    adata.write(out_file)
