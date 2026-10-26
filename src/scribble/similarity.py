#!/usr/bin/env python

from pathlib import Path


def run_similarity(args):
    import random
    import scanpy as sc
    import pandas as pd
    import numpy as np
    import seaborn as sns
    import matplotlib.pyplot as plt
    from scipy import sparse
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

    input_file = Path(args.input)
    pearson_png = PLOT_DIR / f"{input_file.stem}_similarity_pearson.png"
    pearson_tsv = TABLE_DIR / f"{input_file.stem}_similarity_pearson.tsv"
    spearman_png = PLOT_DIR / f"{input_file.stem}_similarity_spearman.png"
    spearman_tsv = TABLE_DIR / f"{input_file.stem}_similarity_spearman.tsv"
    manifest_file = (MANIFEST_DIR / f"{input_file.stem}_similarity_manifest.json")

    pseudobulks = []
    metadata_rows = []

    for infile in args.inputs:

        print(f"Loading {infile}")

        adata = sc.read(infile)

        if args.sample_col not in adata.obs:
            raise ValueError(
                f"{args.sample_col} not in obs"
            )

        if args.batch_col:
            if args.batch not in adata.obs:
                raise ValueError(
                    f"{args.batch_col} not present in obs"
                )


        counts = adata.X

        if sparse.issparse(counts):
            counts = counts.tocsr()

        samples = (
            adata.obs[args.sample_col]
            .astype(str)
        )
        combined_ids = samples

        if args.batch:
            batches = (
                adata.obs[args.batch_col]
                .astype(str)
            )

            combined_ids = (
                batches + ":" + samples
            )



        for sample in sorted(combined_ids.unique()):

            mask = (
                combined_ids == sample
            ).values

            sample_counts = np.asarray(
                counts[mask].sum(axis=0)
            ).ravel()

            pseudobulks.append(
                pd.Series(
                    sample_counts,
                    index=adata.var_names,
                    name=sample,
                )
            )

            dataset, sample_id = sample.split(
                ":",
                1,
            )

            metadata_rows.append({
                "sample": sample,
                "sample_id": sample_id,
                "dataset": dataset,
                "source_file": str(infile),
                "n_cells": int(mask.sum()),
            })

    pseudobulk_df = pd.concat(
        pseudobulks,
        axis=1
    )

    metadata_df = pd.DataFrame(
        metadata_rows
    )

    metadata_df.to_csv(
        metadata_tsv,
        sep="\t",
        index=False,
    )

    pb = sc.AnnData(
        X=pseudobulk_df.T.values
    )

    pb.obs_names = (
        pseudobulk_df.columns
    )

    pb.var_names = (
        pseudobulk_df.index
    )

    sc.pp.normalize_total(pb)

    sc.pp.log1p(pb)

    sc.pp.highly_variable_genes(
        pb,
        n_top_genes=args.hvgs,
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

    g = sns.clustermap(
        pearson,
        cmap="vlag",
        center=0,
        figsize=(10, 10),
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
