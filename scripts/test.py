#!/usr/bin/env python
from pathlib import Path
import scanpy as sc
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.image as mpimg
import random
from scribble.import_data import setup_environment
from scribble.refine import restore_counts
import os
import argparse
from scipy import sparse

def summarize_expression(adata, label="adata"):
    X = adata.X

    if sparse.issparse(X):
        values = X.data
        nnz = X.nnz
    else:
        values = np.asarray(X).ravel()
        nnz = np.count_nonzero(values)

    print(f"\n{label}.X summary")
    print(f"  shape: {adata.shape}")
    print(f"  nonzero entries: {nnz:,}")
    print(f"  min / max: {values.min():.4g} / {values.max():.4g}")
    print(f"  mean (nonzero): {values.mean():.4g}")
    print(f"  first values: {values[:20]}")
    print(f"  integer-like: {np.allclose(values, np.rint(values))}")
    print(f"  log1p metadata: {adata.uns.get('log1p', 'absent')}")


parser = argparse.ArgumentParser()

parser.add_argument("--input_file", type=str, required=True)
parser.add_argument("--plot_file", type=str, required=True)
parser.add_argument("--markers", nargs="+", required=True)
parser.add_argument("--samples", nargs="+", default=None)
parser.add_argument("--dotplot", action="store_true")
parser.add_argument("--umap", action="store_true")


args = parser.parse_args()

# Import adata
adata = sc.read(args.input_file)
summarize_expression(adata, "Loaded data")

try:
    adata_plot = restore_counts(adata)
    summarize_expression(adata_plot, "After restore_counts")

    # preserve annotations/embeddings
    adata_plot.obsm = adata.obsm.copy()
    adata_plot.obs = adata.obs.copy()
    adata_plot.uns = adata.uns.copy()

    adata = adata_plot
    adata.uns.pop("log1p", None)

    sc.pp.normalize_total(adata)
    sc.pp.log1p(adata)
    adata.raw = adata

    print(
        f"Restored full feature space: "
        f"{adata.n_vars:,} genes"
    )

except Exception as e:

    print(
        f"Could not restore counts: {e}"
    )


# Optional sample filtering
if args.samples is not None:
    if "sample" not in adata.obs.columns:
        raise ValueError(
            "Column 'sample' not found in adata.obs. "
            f"Available columns: {list(adata.obs.columns)}"
        )

    before_n = adata.n_obs

    adata = adata[adata.obs["sample"].isin(args.samples)].copy()

    after_n = adata.n_obs

    print(
        f"Subsetted to samples: {', '.join(args.samples)} "
        f"({after_n}/{before_n} cells retained)"
    )

    if after_n == 0:
        raise ValueError(
            f"No cells found for requested samples: {args.samples}"
        )


# Keep only markers present in the dataset
requested_markers = args.markers
available_markers = [g for g in requested_markers if g in adata.var_names]
missing_markers = [g for g in requested_markers if g not in adata.var_names]

if missing_markers:
    print(
        f"Warning: {len(missing_markers)} marker(s) not found and will be skipped: "
        + ", ".join(missing_markers)
    )

if not available_markers:
    raise ValueError("None of the requested markers were found in adata.var_names.")

if args.dotplot:
    sc.pl.dotplot(
        adata,
        available_markers,
        standard_scale="var",
        groupby="refine_label",
        dendrogram=False
    )
    plt.savefig(args.plot_file, dpi=300, bbox_inches="tight")
    plt.close()

if args.umap:
    for gene in available_markers:
        sc.pl.umap(
            adata,
            color=gene,
            cmap="Reds"
        )

        out_file = Path(args.plot_file).with_name(
            f"{Path(args.plot_file).stem}_{gene}.png"
        )
        plt.savefig(out_file, dpi=300, bbox_inches="tight")
        plt.close()
