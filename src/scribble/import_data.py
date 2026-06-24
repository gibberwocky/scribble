#!/usr/bin/env python
from pathlib import Path


# ---------- Environment setup ----------
def setup_environment(sc, np, random, plot_dir):
    sc.settings.verbosity = 3
    #scv.settings.verbosity = 3
    sc.logging.print_header()
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
def add_metadata(pd, adata, metadata_file):
    print(f"Reading metadata: {metadata_file}")

    meta = pd.read_excel(metadata_file, sheet_name="meta")

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


# ----------- Inflection point -----------
def get_inflection(np, df, lower=100):
    df_fit = df[df["total"] >= lower].copy()
    # log-transform
    log_total = np.log10(df_fit["total"].values)
    log_rank = np.log10(df_fit["rank"].values)
    # first derivative
    d1n = np.diff(log_total) / np.diff(log_rank)
    # locate minimum slope
    right_edge = np.argmin(d1n)
    inflection = 10 ** log_total[right_edge]
    return inflection


# ---------- Process one sample ----------
def process_sample(sc, sp, np, pd, sample, cellranger_dir, velo_dir, plot_dir):
    from scribble.plots import qc_hexbin_ax, knee_plot_ax
    import matplotlib.pyplot as plt

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

    # Use gene symbol for name and ensure unique
    adata.var_names = adata.var.index
    adata.var_names_make_unique()

    # Layers
    spliced = sp.csr_matrix(adata.layers["spliced"], dtype=np.float32)
    unspliced = sp.csr_matrix(adata.layers["unspliced"], dtype=np.float32)
    counts = spliced + unspliced

    # Build final object
    adata_clean = sc.AnnData(
        X=counts.copy(),
        obs=adata.obs.copy(),
        var=adata.var.copy()
    )
    adata_clean.layers["counts"] = counts.copy()
    adata_clean.layers["spliced"] = spliced
    adata_clean.layers["unspliced"] = unspliced

    # Transfer MT metrics from adta_10X to adata
    meta_10X = adata_10X.obs.loc[adata_clean.obs_names]
    adata_clean.obs["pct_counts_mt"] = meta_10X["pct_counts_mt"].values
    adata_clean.obs["total_counts"] = meta_10X["total_counts"].values
    adata_clean.obs["n_genes_by_counts"] = meta_10X["n_genes_by_counts"].values

    # Data for knee plot and inflection point
    knee_df = pd.DataFrame({
        "total": adata_clean.obs["total_counts"]
    }).sort_values("total", ascending=False).reset_index(drop=True)
    knee_df["rank"] = np.arange(1, len(knee_df) + 1)
    inflection = get_inflection(np, knee_df, lower=inflection_lower)

    # Multi-panel QC figure
    fig, axes = plt.subplots(1, 3, figsize=(15,5), constrained_layout=True)
    axes = axes.flatten()

    # ---------- 1. counts vs genes ----------
    hb = qc_hexbin_ax(
        axes[0],
        x=adata_clean.obs["total_counts"].values,
        y=adata_clean.obs["n_genes_by_counts"].values,
        xlabel="Total counts (log scale)",
        ylabel="Number of genes (log scale)",
        log_x=True,
        log_y=True
    )
    axes[0].set_title("Library complexity")

    # ---------- 2. knee plot ----------
    knee_plot_ax(axes[1], knee_df, inflection)

    # ---------- 3. counts vs MT ----------
    qc_hexbin_ax(
        axes[2],
        x=adata_clean.obs["total_counts"].values,
        y=adata_clean.obs["pct_counts_mt"].values,
        xlabel="Total counts (log scale)",
        ylabel="% mitochondrial counts",
        log_x=True,
        log_y=False
    )
    axes[2].set_title("MT content")

    # ---------- shared colorbar ----------
    cbar = fig.colorbar(
        hb,
        ax=axes,
        location="right",
        fraction=0.05,
        pad=0.04
    )
    cbar.set_label("log10(cell density)")

    # ---------- layout fix ----------
    plt.savefig(plot_dir / f"{sample}_qc_panel.png", dpi=300, bbox_inches="tight")
    plt.close()

    # Filter cells above inflection point
    adata_clean = adata_clean[adata_clean.obs["total_counts"] > inflection].copy()

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
    PLOT_DIR = PROJECT_DIR / "scribble/plots"
    ADATA_DIR = PROJECT_DIR / "scribble/adata"
    TABLE_DIR = PROJECT_DIR / "scribble/tables"

    PLOT_DIR.mkdir(exist_ok=True, parents=True)
    ADATA_DIR.mkdir(exist_ok=True, parents=True)
    TABLE_DIR.mkdir(exist_ok=True, parents=True)

    setup_environment(sc, np, random, PLOT_DIR)

    samples = get_samples(VELO_DIR)
    print(f"Detected samples: {samples}")

    adatas = []

    for sample in samples:
        ad = process_sample(sc, sp, np, pd, sample, CELLRANGER_DIR, VELO_DIR, PLOT_DIR)
        ad.obs["sample"] = sample
        adatas.append(ad)

    print("Concatenating...")
    adata = sc.concat(adatas, label="sample", keys=samples)

    adata = add_metadata(pd, adata, args.metadata_file)

    out_file = ADATA_DIR / "combined.h5ad"
    print(f"Saving to {out_file}")

    adata.write(out_file)
