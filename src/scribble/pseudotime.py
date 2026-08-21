#!/usr/bin/env python

from pathlib import Path
import re

import scanpy as sc
import scFates as scf
import numpy as np
import pandas as pd
import decoupler as dc
import omnipath as op
import harmonypy as hm
import matplotlib.pyplot as plt

from scribble.refine import restore_counts


def safe_sheet_name(name):
    name = re.sub(r"[/\\?*\\[\\]:]", "_", name)
    return name[:31]


def load_dorothea():

    net = op.interactions.Dorothea.get(
        organism="human",
        dorothea_level=["A", "B", "C"],
        genesymbols=True
    )

    net["source"] = net["source_genesymbol"]
    net["target"] = net["target_genesymbol"]
    net["weight"] = np.where(
        net["is_inhibition"],
        -1,
        1
    )

    net = net[
        ["source", "target", "weight"]
    ].copy()

    net = net.dropna()
    net = net.drop_duplicates(
        subset=["source", "target"]
    )

    return net


def build_lineage(adata, clusters):

    return adata[
        adata.obs["refine_label"]
        .astype(str)
        .isin(clusters)
    ].copy()


def build_marker_dataset(
    adata_lineage,
    celltype,
    tf_list
):

    adata_markers = adata_lineage.copy()

    sc.pp.normalize_total(adata_markers)
    sc.pp.log1p(adata_markers)

    sc.tl.rank_genes_groups(
        adata_markers,
        groupby=celltype,
        method="wilcoxon"
    )

    markers = sc.get.rank_genes_groups_df(
        adata_markers,
        group=None
    )

    marker_tfs = markers[
        markers["names"].isin(tf_list)
    ].copy()

    adata_markers.uns["marker_tfs"] = (
        marker_tfs
    )

    return adata_markers


def build_pseudobulk(
    adata_lineage,
    groupby
):

    pbs = dc.pp.pseudobulk(
        adata_lineage,
        sample_col="sample",
        groups_col=groupby
    )

    pbs = pbs[
        np.asarray(
            pbs.X.sum(axis=1)
        ).ravel() > 0
    ].copy()

    sc.pp.normalize_total(
        pbs,
        target_sum=1e4
    )

    sc.pp.log1p(pbs)

    return pbs


def build_trajectory_dataset(
    adata_lineage,
    hvgs,
    npcs
):

    adata_traj = adata_lineage.copy()

    sc.pp.normalize_total(adata_traj)
    sc.pp.log1p(adata_traj)

    adata_traj.raw = adata_traj.copy()

    sc.pp.highly_variable_genes(
        adata_traj,
        n_top_genes=hvgs,
        flavor="seurat_v3",
        layer="counts"
    )

    adata_traj = adata_traj[
        :,
        adata_traj.var.highly_variable
    ].copy()

    sc.pp.scale(
        adata_traj,
        max_value=10
    )

    sc.tl.pca(
        adata_traj,
        n_comps=npcs
    )

    return adata_traj


def run_harmony_if_required(
    adata,
    batch
):

    nbatches = (
        adata.obs[batch]
        .astype(str)
        .nunique()
    )

    if nbatches < 2:

        print(
            "Single batch detected. "
            "Using PCA embedding."
        )

        return "X_pca"

    print(
        f"Running Harmony on "
        f"{nbatches} batches "
        f"using '{batch}'."
    )

    ho = hm.run_harmony(
        adata.obsm["X_pca"],
        adata.obs,
        batch
    )

    adata.obsm["X_pca_harmony"] = (
        ho.Z_corr
    )

    return "X_pca_harmony"


def infer_trajectory(
    adata_traj,
    representation,
    root,
    neighbors,
    nodes,
    celltype,
    plotdir,
    trajectory_name
):

    sc.pp.neighbors(
        adata_traj,
        n_neighbors=neighbors,
        use_rep=representation
    )

    sc.tl.umap(
        adata_traj,
        min_dist=0.3
    )

    scf.tl.tree(
        adata_traj,
        method="epg",
        Nodes=nodes,
        use_rep=representation
    )

    sc.pl.umap(
        adata_traj,
        color=celltype,
        legend_loc="right margin",
        frameon=False,
        show=False
    )

    plt.savefig(
        plotdir /
        f"umap_{trajectory_name}.png",
        dpi=300,
        bbox_inches="tight"
    )

    plt.close()

    scf.pl.graph(
        adata_traj,
        basis="umap"
    )

    plt.savefig(
        plotdir /
        f"umap_{trajectory_name}_tree.png",
        dpi=300,
        bbox_inches="tight"
    )

    plt.close()

    adata_traj.obs["root_cells"] = (
        adata_traj.obs["refine_label"]
        .astype(str)
        .isin([root])
    )

    scf.tl.root(
        adata_traj,
        "root_cells"
    )

    # This patches a bug related to epigraph format mismatch
    adata_traj.uns["epg"]["Edges"] = [
        adata_traj.uns["epg"]["Edges"]
    ]

    scf.tl.pseudotime(adata_traj)

    sc.pl.umap(
        adata_traj,
        color="t",
        cmap="viridis",
        show=False
    )

    plt.savefig(
        plotdir /
        f"umap_{trajectory_name}_pseudo.png",
        dpi=300,
        bbox_inches="tight"
    )

    plt.close()

    return adata_traj


def synchronise_objects(
    adata_markers,
    adata_traj
):

    adata_markers.obsm["X_umap"] = (
        adata_traj.obsm["X_umap"]
    )

    adata_markers.obs["t"] = (
        adata_traj.obs["t"]
    )

    adata_markers.obs["seg"] = (
        adata_traj.obs["seg"]
    )


def run_pseudotime(args):

    import scanpy as sc
    import numpy as np
    import pandas as pd
    import random

    from scribble.import_data import setup_environment

    PROJECT_DIR = Path(args.project_dir)
    ADATA_DIR = PROJECT_DIR / "scribble/adata"

    PLOT_DIR = PROJECT_DIR / "scribble/plots/pseudotime"
    PLOT_DIR.mkdir(parents=True, exist_ok=True)

    setup_environment(sc, np, random, PLOT_DIR)

    clusters = [
        str(x)
        for x in args.clusters
    ]

    trajectory_name = "_".join(clusters)

    root = (
        str(args.root)
        if args.root is not None
        else clusters[0]
    )

    print("Loading AnnData")
    input_file = Path(args.input)
    adata = sc.read(input_file)

    print(f"Building lineage: {trajectory_name}")

    adata_lineage = build_lineage(adata, clusters)
    adata_lineage = restore_counts(adata_lineage)
    adata_lineage.layers["counts"] = (adata_lineage.X.copy())

    tf_list = []
    if args.tf_list:
        tf_list = pd.read_csv(args.tf_list, header=None)[0].tolist()

    print("Building marker dataset")
    adata_markers = build_marker_dataset(adata_lineage, args.celltype, tf_list)

    print("Building pseudobulk dataset")
    pbs = build_pseudobulk(adata_lineage, args.celltype)

    print("Inferring TF activity")
    net = load_dorothea()
    dc.mt.ulm(data=pbs, net=net, tmin=5, verbose=True)

    adata_traj = build_trajectory_dataset(adata_lineage, args.hvgs, args.npcs)

    representation = (
        run_harmony_if_required(
            adata_traj,
            args.batch
        )
    )

    print("Inferring trajectory")
    adata_traj = infer_trajectory(
        adata_traj=adata_traj,
        representation=representation,
        root=root,
        neighbors=args.neighbors,
        nodes=args.nodes,
        celltype=args.celltype,
        plotdir=PLOT_DIR,
        trajectory_name=trajectory_name
    )

    synchronise_objects(adata_markers, adata_traj)

    if args.association:
        print("Running trajectory association analysis")
        scf.tl.test_association(adata_traj, n_jobs=args.n_jobs)

    print("Writing outputs")
    adata_markers.write(ADATA_DIR / f"{trajectory_name}_markers.h5ad")
    adata_traj.write(ADATA_DIR / f"{trajectory_name}_trajectory.h5ad")
    pbs.write(ADATA_DIR / f"{trajectory_name}_pseudobulk.h5ad")
