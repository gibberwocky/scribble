#!/usr/bin/env python

from pathlib import Path


def pseudobulk_samples(
    adata,
    sample_col,
    batch_col=None,
    min_cells_per_sample=20,
):
    """
    Aggregate raw counts into sample-level pseudobulk profiles.

    Parameters
    ----------
    adata : AnnData
        AnnData object containing raw counts in
        adata.layers["counts"].

    sample_col : str
        Column defining biological samples.

    batch_col : str
        Column defining dataset, cohort or batch.
        Used to generate unique sample identifiers.

    min_cells_per_sample : int
        Minimum cells required to retain a sample.

    Returns
    -------
    pb : pd.DataFrame

        Rows:
            dataset:sample

        Columns:
            genes

        Values:
            integer pseudobulk counts

    meta_pb : pd.DataFrame

        Sample-level metadata aligned to pb.
    """

    import pandas as pd
    import numpy as np

    # --------------------------------------------------
    # Determine where raw counts are stored
    # --------------------------------------------------

    if "counts" in adata.layers:
        count_source = "layer"
    elif np.issubdtype(adata.X.dtype, np.integer):
        count_source = "X"
    else:
        raise ValueError(
            "Raw counts not found. Expected either "
            "adata.layers['counts'] or integer counts in adata.X"
        )


    if sample_col not in adata.obs:
        raise ValueError(
            f"{sample_col} not found in adata.obs"
        )

    if batch_col:
        if batch_col not in adata.obs:
            raise ValueError(
                f"{batch_col} not found in adata.obs"
            )


    obs = adata.obs.copy()

    obs["_sample"] = (
        obs[sample_col]
        .astype(str)
    )

    if batch_col:
        obs["_batch"] = (
            obs[batch_col]
            .astype(str)
        )

        obs["_sample_id"] = (
            obs["_batch"]
            + ":"
            + obs["_sample"]
        )
    else:
        obs["_sample_id"] = (
            obs["_sample"]
        )

    # --------------------------------------------------
    # Filter samples with insufficient cells
    # --------------------------------------------------

    cell_counts = (
        obs["_sample_id"]
        .value_counts()
    )

    keep_samples = cell_counts[
        cell_counts >= min_cells_per_sample
    ].index

    keep_mask = (
        obs["_sample_id"]
        .isin(keep_samples)
    )

    subset = adata[keep_mask].copy()

    if subset.n_obs == 0:
        return None, None

    subset.obs["_sample_id"] = (
        obs.loc[
            keep_mask,
            "_sample_id"
        ].values
    )

    pseudobulk = []
    sample_meta = []

    keep_samples = sorted(keep_samples)

    for sample_id in keep_samples:

        sample_mask = (
            subset.obs["_sample_id"]
            == sample_id
        )

        sample_adata = subset[sample_mask]

        if count_source == "layer":
            X = sample_adata.layers["counts"]
        else:
            X = sample_adata.X

        summed = np.asarray(X.sum(axis=0)).ravel()
        pseudobulk.append(summed)

        meta_row = {
            "sample": sample_id,
            "dataset": (
                sample_id.split(":", 1)[0]
                if batch_col
                else "dataset"
            ),
            "n_cells": int(
                sample_adata.n_obs
            ),
            "library_size": int(
                summed.sum()
            ),
        }

        sample_meta.append(meta_row)

    pb = pd.DataFrame(
        pseudobulk,
        columns=subset.var_names,
        index=keep_samples,
    )

    pb = (
        pb
        .round()
        .astype(np.int64)
    )

    meta_pb = (
        pd.DataFrame(sample_meta)
        .set_index("sample")
        .loc[pb.index]
    )

    return pb, meta_pb


def run_similarity(args):
    import random
    from datetime import datetime
    import scanpy as sc
    import pandas as pd
    import numpy as np
    import seaborn as sns
    import matplotlib.pyplot as plt
    from scribble.refine import restore_counts
    from scribble.import_data import setup_environment
    from scribble.utils import (
            serializable_args,
            write_run_manifest,
        )

    # ----------------------------
    # Setup
    # ----------------------------
    PROJECT_DIR = Path(args.project_dir)
    MANIFEST_DIR = PROJECT_DIR / "scribble/manifests"
    PLOT_DIR = PROJECT_DIR / "scribble/plots"
    TABLE_DIR = PROJECT_DIR / "scribble/tables"

    setup_environment(sc, np, random, PLOT_DIR)
    tag = datetime.now().strftime("%Y%m%d_%H%M%S")

    pearson_png = PLOT_DIR / f"similarity_pearson_{tag}.png"
    pearson_tsv = TABLE_DIR / f"similarity_pearson_{tag}.tsv"
    spearman_png = PLOT_DIR / f"similarity_spearman_{tag}.png"
    spearman_tsv = TABLE_DIR / f"similarity_spearman_{tag}.tsv"
    manifest_file = (MANIFEST_DIR / f"similarity_manifest_{tag}.json")
    metadata_tsv = (TABLE_DIR / f"similarity_{tag}_metadata.tsv")

    pseudobulks = []
    metadata_frames = []

    for infile in args.inputs:

        print(f"Loading {infile}")

        adata = sc.read(infile)

        adata = restore_counts(adata)
        print("X dtype:", adata.X.dtype)
        print("layers:", list(adata.layers.keys()))

        pb_df, meta_df = pseudobulk_samples(
            adata=adata,
            sample_col=args.sample_col,
            batch_col=args.batch_col,
            min_cells_per_sample=args.min_cells_per_sample,
        )

        if pb_df is None:

            print(
                f"Skipping {infile}: "
                "no samples passed filtering"
            )

            continue

        pseudobulks.append(pb_df)

        meta_df["source_file"] = str(infile)
        meta_df.index.name = "sample"

        metadata_frames.append(meta_df)

    if len(pseudobulks) == 0:

        raise ValueError(
            "No pseudobulk samples generated"
        )

    pseudobulk_df = pd.concat(
        pseudobulks,
        axis=0,
    )

    print(
        f"Pseudobulk samples: "
        f"{pseudobulk_df.shape[0]}"
    )

    print(
        f"Genes: "
        f"{pseudobulk_df.shape[1]}"
    )

    metadata_df = pd.concat(
        metadata_frames,
        axis=0,
    )

    metadata_df = metadata_df.loc[
        ~metadata_df.index.duplicated()
    ]

    metadata_df.to_csv(
        metadata_tsv,
        sep="\t",
        index=False,
    )

    pb = sc.AnnData(
        X=pseudobulk_df.values
    )

    pb.obs_names = (
        pseudobulk_df.index
    )

    pb.var_names = (
        pseudobulk_df.columns
    )

    pb.obs = metadata_df.loc[
        pb.obs_names
    ].copy()

    sc.pp.normalize_total(pb)

    sc.pp.log1p(pb)

    n_hvgs = min(
        args.hvgs,
        pb.n_vars,
    )

    print(
        f"HVGs retained: "
        f"{n_hvgs}"
    )

    sc.pp.highly_variable_genes(
        pb,
        n_top_genes=n_hvgs,
        flavor="seurat",
    )

    pb = pb[
        :,
        pb.var.highly_variable
    ].copy()

    expr = pd.DataFrame(
        pb.X.T,
        index=pb.var_names,
        columns=pb.obs_names,
    )

    pearson = expr.corr(
        method="pearson"
    )

    spearman = expr.corr(
        method="spearman"
    )

    pearson.to_csv(
        pearson_tsv,
        sep="\t",
    )

    spearman.to_csv(
        spearman_tsv,
        sep="\t",
    )


    sample_colors = None

    if (
        args.batch_col is not None
        and "dataset" in pb.obs.columns
    ):

        palette = dict(
            zip(
                pb.obs["dataset"].unique(),
                sns.color_palette(
                    "tab10",
                    pb.obs["dataset"].nunique(),
                ),
            )
        )

        sample_colors = (
            pb.obs["dataset"]
            .map(palette)
        )

    g = sns.clustermap(
        pearson,
        cmap="vlag",
        center=0,
        figsize=(10, 10),
        row_colors=sample_colors,
        col_colors=sample_colors,
    )

    g.fig.savefig(
        pearson_png,
        dpi=300,
        bbox_inches="tight",
    )

    plt.close(g.fig)

    g = sns.clustermap(
        spearman,
        cmap="vlag",
        center=0,
        figsize=(10, 10),
        row_colors=sample_colors,
        col_colors=sample_colors,
    )

    g.fig.savefig(
        spearman_png,
        dpi=300,
        bbox_inches="tight",
    )

    plt.close(g.fig)

    write_run_manifest(
        manifest_file=manifest_file,
        module="similarity",
        args=args,
        inputs={
            "inputs": [
                str(x)
                for x in args.inputs
            ]
        },
        outputs={
            "pearson_tsv":
                str(pearson_tsv),
            "spearman_tsv":
                str(spearman_tsv),
            "pearson_png":
                str(pearson_png),
            "spearman_png":
                str(spearman_png),
            "metadata_tsv":
                str(metadata_tsv),
        },
    )
