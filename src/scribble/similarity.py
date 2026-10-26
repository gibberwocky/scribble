#!/usr/bin/env python

from pathlib import Path


def pseudobulk_samples(
    adata,
    sample_col,
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

    obs = adata.obs.copy()

    obs["_sample"] = (
        obs[sample_col]
        .astype(str)
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

    spearman_png = PLOT_DIR / f"distance_spearman_{tag}.png"
    spearman_tsv = TABLE_DIR / f"distance_spearman_{tag}.tsv"
    nearest_tsv = TABLE_DIR / f"nearest_spearman_{tag}.tsv"
    manifest_file = (MANIFEST_DIR / f"similarity_manifest_{tag}.json")
    metadata_tsv = (TABLE_DIR / f"similarity_{tag}_metadata.tsv")

    pseudobulks = []
    metadata_frames = []
    dataset_samples = {}

    if len(args.inputs) != len(args.datasets):
        raise ValueError(
            "--inputs and --datasets "
            "must contain the same number of entries"
        )


    for infile, dataset_name in zip(
        args.inputs,
        args.datasets,
    ):

        print(f"Loading {infile}")

        adata = sc.read(infile)

        adata = restore_counts(adata)
        adata.layers["counts"] = adata.X.copy()

        pb_df, meta_df = pseudobulk_samples(
            adata=adata,
            sample_col=args.sample_col,
            min_cells_per_sample=args.min_cells_per_sample,
        )

        if pb_df is None:
            print(f"Skipping {infile}: no samples passed filtering")
            continue

        pb_df.index = [
            f"{dataset_name}:{x}"
            for x in pb_df.index
        ]

        meta_df.index = [
            f"{dataset_name}:{x}"
            for x in meta_df.index
        ]

        meta_df["dataset"] = dataset_name

        dataset_samples.setdefault(
            dataset_name,
            []
        ).extend(
            pb_df.index.tolist()
        )

        print(f"{dataset_name}: {len(pb_df.index)} samples")
        pseudobulks.append(pb_df)

        meta_df["dataset"] = dataset_name
        meta_df["source_file"] = str(infile)
        meta_df.index.name = "sample"

        metadata_frames.append(meta_df)

    print("\n[DATASETS]")
    print(dataset_samples.keys())

    if len(pseudobulks) == 0:
        raise ValueError("No pseudobulk samples generated")

    shared_genes = set(pseudobulks[0].columns)

    for pb in pseudobulks[1:]:
        shared_genes &= set(pb.columns)

    shared_genes = sorted(shared_genes)

    pseudobulks = [
        pb.loc[:, shared_genes]
        for pb in pseudobulks
    ]

    pseudobulk_df = pd.concat(
        pseudobulks,
        axis=0,
    )
    pseudobulk_df = pseudobulk_df.fillna(0)

    print(f"Pseudobulk samples: {pseudobulk_df.shape[0]}")
    print(f"Genes: {pseudobulk_df.shape[1]}")

    metadata_df = pd.concat(metadata_frames, axis=0)
    metadata_df = metadata_df.loc[~metadata_df.index.duplicated()]
    metadata_df.to_csv(metadata_tsv, sep="\t", index=False)

    pb = sc.AnnData(X=pseudobulk_df.values)
    pb.obs_names = (pseudobulk_df.index)
    pb.var_names = (pseudobulk_df.columns)
    pb.obs = metadata_df.loc[pb.obs_names].copy()

    sc.pp.normalize_total(pb)
    sc.pp.log1p(pb)

    if args.hvgs:
        n_hvgs = min(args.hvgs, pb.n_vars)
        print(f"HVGs retained: {n_hvgs}")
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

    spearman = expr.corr(method="spearman")
    distance = (1.0 - spearman)
    distance.to_csv(spearman_tsv, sep="\t")

    dataset_names = list(
        dataset_samples.keys()
    )

    if len(dataset_names) < 2:
        raise ValueError(
            "Need at least two datasets"
        )

    reference_dataset = dataset_names[0]

    other_samples = []

    for ds in dataset_names[1:]:
        other_samples.extend(
            dataset_samples[ds]
        )

    nearest_rows = []

    for sample in dataset_samples[
        reference_dataset
    ]:

        d = distance.loc[
            sample,
            other_samples
        ]

        nearest_sample = d.idxmin()

        nearest_distance = float(
            d.min()
        )

        nearest_rows.append({
            "sample": sample,
            "nearest_sample": nearest_sample,
            "distance": nearest_distance,
        })

    nearest_df = pd.DataFrame(
        nearest_rows
    )

    print(
        "\nNearest cross-dataset matches"
    )

    for _, row in (
        nearest_df
        .sort_values("distance")
        .iterrows()
    ):

        print(
            f"{row['sample']} "
            f"-> "
            f"{row['nearest_sample']} "
            f"(distance={row['distance']:.3f})"
        )

    nearest_df.to_csv(
        nearest_tsv,
        sep="\t",
        index=False,
    )



    sample_colors = None
    if "dataset" in pb.obs.columns:
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
        1-spearman,
        cmap="viridis",
        vmin = np.quantile(1-spearman.values, 0.05),
        vmax = np.quantile(1-spearman.values, 0.95),
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
            "spearman distance tsv": str(spearman_tsv),
            "spearman nearest tsv": str(nearest_tsv),
            "spearman distance png": str(spearman_png),
            "metadata tsv": str(metadata_tsv),
        },
    )
