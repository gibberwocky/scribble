#!/usr/bin/env python

import anndata as ad
import numpy as np
import scanpy as sc
import scvelo as scv
import scipy.sparse as sp
from pathlib import Path


def load_velocyto_sample(sample, velo_dir):

    sample_dir = Path(velo_dir) / sample

    loom_files = list(sample_dir.glob("*.loom"))

    if len(loom_files) == 0:
        raise FileNotFoundError(
            f"No loom file found in {sample_dir}"
        )

    if len(loom_files) > 1:
        raise RuntimeError(
            f"Multiple loom files found in {sample_dir}: {loom_files}"
        )

    loom_file = loom_files[0]

    print(f"Using loom file for {sample}: {loom_file}")

    adata = sc.read(loom_file, cache=True)

    adata.var_names = adata.var.index
    adata.var_names_make_unique()

    adata.obs_names = (
        adata.obs_names
        .str.replace(r"^.*:", "", regex=True)
        .str.replace(r"x$", f"_{sample}", regex=True)
    )

    spliced = sp.csr_matrix(
        adata.layers["spliced"],
        dtype=np.float32,
    )

    unspliced = sp.csr_matrix(
        adata.layers["unspliced"],
        dtype=np.float32,
    )

    out = ad.AnnData(
        X=spliced + unspliced,
        obs=adata.obs.copy(),
        var=adata.var.copy(),
    )

    out.layers["spliced"] = spliced
    out.layers["unspliced"] = unspliced

    return out


def load_velocyto(samples, velo_dir):

    adatas = []
    for sample in samples:
        adatas.append(load_velocyto_sample(sample, velo_dir))

    return sc.concat(adatas, label="sample", keys=samples)


def build_monod(
    adata_loom,
    adata_traj,
    min_shared_counts=20
):

    shared_cells = (
        adata_traj.obs_names
        .intersection(adata_loom.obs_names)
    )

    print(
        f"Found {len(shared_cells)} overlapping cells"
    )

    if len(shared_cells) == 0:
        raise ValueError(
            "No shared cells between trajectory and loom data"
        )

    adata_monod = adata_loom[
        shared_cells,
        :
    ].copy()

    adata_monod.obs = (
        adata_traj.obs
        .loc[shared_cells]
        .copy()
    )

    scv.pp.filter_genes(
        adata_monod,
        min_shared_counts=min_shared_counts
    )

    return adata_monod


def report_balance(adata):

    spliced_total = adata.layers["spliced"].sum()
    unspliced_total = adata.layers["unspliced"].sum()

    print(
        f"Global ratio: "
        f"{unspliced_total / spliced_total:.4f}"
    )

    spliced_per_cell = np.asarray(
        adata.layers["spliced"].sum(axis=1)
    ).ravel()

    unspliced_per_cell = np.asarray(
        adata.layers["unspliced"].sum(axis=1)
    ).ravel()

    ratio = (
        unspliced_per_cell /
        np.maximum(spliced_per_cell, 1)
    )

    print(f"Median ratio: {np.median(ratio):.4f}")
    print(f"Mean ratio: {np.mean(ratio):.4f}")


def run_premonod(args):

    from pathlib import Path
    from scribble.import_data import setup_environment

    PROJECT_DIR = Path(args.project_dir)
    VELO_DIR = Path(args.velocyto_dir)
    LOOM_DIR = Path(args.loom_dir)
    PLOT_DIR = PROJECT_DIR / "scribble/plots"
    ADATA_DIR = PROJECT_DIR / "scribble/adata"
    TABLE_DIR = PROJECT_DIR / "scribble/tables"

    setup_environment(sc, np, random, PLOT_DIR)
    scv.settings.verbosity = 3
    sc.logging.print_versions()

    #  Import velocyto loom files
    adata_loom = load_velocyto(args.samples, VELO_DIR)

    # Import h5ad AnnData (traj, markers, pbs)
    input_file = Path(args.input)
    adata_traj = sc.read(args.input)

    adata_monod = build_monod(
        adata_loom,
        adata_traj,
        args.min_shared_counts
    )

    report_balance(adata_monod)

    ad.settings.allow_write_nullable_strings = True

    outfile = (
        Path(args.output)
        if args.output
        else datdir /
        f"{args.trajectory}_monod.h5ad"
    )

    adata_monod.write(outfile)

    print(f"Saved {outfile}")
