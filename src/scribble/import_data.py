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
def process_sample(
    sc,
    sp,
    np,
    pd,
    sample,
    cellranger_dir,
    cellbender_dir,
    velo_dir,
    plot_dir,
    inflection_lower
):

    from scribble.plots import qc_hexbin_ax, knee_plot_ax
    import matplotlib.pyplot as plt

    print(f"Processing {sample}")

    # Read velocyto data
    loom_file = list((velo_dir / sample).glob("*.loom"))[0]
    adata_velo = sc.read(loom_file)
    adata_velo.obs_names = clean_barcodes(adata_velo.obs_names, sample)
    adata_velo.var_names = adata_velo.var.index
    adata_velo.var_names_make_unique()


    # Read CellBender or CellRanger data
    adata_primary = None
    # CellBender preferred
    if (cellbender_dir is not None and (cellbender_dir / sample / f"{sample}_filtered.h5").exists()):
        print(f"{sample}: using CellBender counts for X")
        adata_primary = sc.read_10x_h5(cellbender_dir / sample / f"{sample}_filtered.h5")
    # Fallback to Cell Ranger
    elif (cellranger_dir is not None and (cellranger_dir / sample / "outs/filtered_feature_bc_matrix").exists()):
        print(f"{sample}: using Cell Ranger counts for X")
        adata_primary = sc.read_10x_mtx(cellranger_dir / sample / "outs/filtered_feature_bc_matrix")
    # Else raise error
    else:
        raise ValueError(f"{sample}: no CellBender or CellRanger count matrix available")


    # Clean barcodes
    adata_primary.obs_names = clean_barcodes(adata_primary.obs_names, sample)

    # Check for shared cells with velocyto data
    shared_cells = (adata_primary.obs_names.intersection(adata_velo.obs_names))
    if len(shared_cells) == 0:
        raise ValueError(
            f"{sample}: no overlapping cells between "
            "primary counts and velocyto data"
        )
    adata_primary = adata_primary[shared_cells].copy()
    adata_velo = adata_velo[shared_cells].copy()

    # MT + QC
    adata_primary.var["mt"] = adata_primary.var_names.str.contains("^MT-|^mt-")
    sc.pp.calculate_qc_metrics(adata_primary, qc_vars=["mt"], inplace=True)

    # Use gene symbol for name and ensure unique
    adata_primary.var_names = adata_primary.var.index
    adata_primary.var_names_make_unique()

    # Velocyto layers
    spliced = sp.csr_matrix(adata_velo.layers["spliced"], dtype=np.float32)
    unspliced = sp.csr_matrix(adata_velo.layers["unspliced"], dtype=np.float32)

    # Build final object
    counts = sp.csr_matrix(adata_primary.X, dtype=np.float32)
    adata_clean = sc.AnnData(
        X=counts.copy(),
        obs=adata_primary.obs.copy(),
        var=adata_primary.var.copy()
    )
    adata_clean.layers["counts"] = counts.copy()
    adata_clean.layers["spliced"] = spliced
    adata_clean.layers["unspliced"] = unspliced

    # Transfer MT metrics from adata_10X to adata
    meta_primary = adata_primary.obs.loc[adata_clean.obs_names]
    adata_clean.obs["pct_counts_mt"] = meta_primary["pct_counts_mt"].values
    adata_clean.obs["total_counts"] = meta_primary["total_counts"].values
    adata_clean.obs["n_genes_by_counts"] = meta_primary["n_genes_by_counts"].values

    # Add provenance
    adata_clean.uns["count_source"] = (
        "cellbender"
        if cellbender_dir is not None
        and (
            cellbender_dir /
            sample /
            f"{sample}_filtered.h5"
        ).exists()
        else "cellranger"
    )

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

# Check sampels present in cellranger and velocyto directories
def get_samples(cellranger_dir, cellbender_dir, velo_dir, requested_samples=None):

    available_cr = {
        d.name for d in cellranger_dir.iterdir()
        if d.is_dir()
    }

    available_velo = {
        d.name for d in velo_dir.iterdir()
        if d.is_dir()
    }

    available_counts = set()
    if cellranger_dir.exists():
        available_counts |= available_cr
    if (
        cellbender_dir is not None
        and cellbender_dir.exists()
    ):
        available_cb = {
            d.name
            for d in cellbender_dir.iterdir()
            if d.is_dir()
        }
        available_counts |= available_cb

    available = (available_counts & available_velo)

    if requested_samples is None:
        samples = sorted(available)
    else:
        missing = sorted(set(requested_samples) - available)

        if missing:
            raise ValueError(
                f"Samples missing from CellBender|CellRanger and/or Velocyto:\n"
                + "\n".join(missing)
            )

        samples = sorted(requested_samples)

    return samples

# ---------- MAIN ENTRY ----------
def run_import(args):
    import scanpy as sc
    import numpy as np
    import pandas as pd
    import random
    import scipy.sparse as sp
    from scribble.utils import write_run_manifest


    PROJECT_DIR = Path(args.project_dir)

    CELLRANGER_DIR = Path(args.cellranger_dir)
    CELLBENDER_DIR = (Path(args.cellbender_dir) if args.cellbender_dir else None)
    VELO_DIR = Path(args.velocyto_dir)
    PLOT_DIR = PROJECT_DIR / "scribble/plots"
    ADATA_DIR = PROJECT_DIR / "scribble/adata"
    TABLE_DIR = PROJECT_DIR / "scribble/tables"
    MANIFEST_DIR = PROJECT_DIR / "scribble/manifests"

    PLOT_DIR.mkdir(exist_ok=True, parents=True)
    ADATA_DIR.mkdir(exist_ok=True, parents=True)
    TABLE_DIR.mkdir(exist_ok=True, parents=True)
    MANIFEST_DIR.mkdir(exist_ok=True, parents=True)

    setup_environment(sc, np, random, PLOT_DIR)


    requested_samples = None

    if args.samples:
        requested_samples = args.samples

    elif args.samples_file:
        with open(args.samples_file) as handle:
            requested_samples = [
                line.strip()
                for line in handle
                if line.strip()
            ]

    samples = get_samples(
        CELLRANGER_DIR,
        CELLBENDER_DIR,
        VELO_DIR,
        requested_samples
    )

    print(f"Detected samples: {samples}")

    adatas = []

    for sample in samples:
        ad = process_sample(sc, sp, np, pd, sample, CELLRANGER_DIR, CELLBENDER_DIR, VELO_DIR, PLOT_DIR, args.inflection_lower)
        ad.obs["sample"] = sample
        adatas.append(ad)

    print("Concatenating...")
    adata = sc.concat(adatas, label="sample", keys=samples)

    adata = add_metadata(pd, adata, args.metadata_file)

    out_file = ADATA_DIR / "combined.h5ad"
    manifest_file = (MANIFEST_DIR / "import_manifest.json")
    adata.uns["scribble_manifest"] = str(manifest_file)

    print(f"Saving to {out_file}")
    adata.write(out_file)

    print(f"Writing run manifest to {manifest_file}")
    write_run_manifest(
        manifest_file=manifest_file,
        module="import",
        args=args,
        inputs={
            "cellranger_dir": str(CELLRANGER_DIR),
            "cellbender_dir": str(CELLBENDER_DIR),
            "velocyto_dir": str(VELO_DIR),
            "metadata_file": str(args.metadata_file),
        },
        outputs={
            "adata": str(out_file),
        },
    )
