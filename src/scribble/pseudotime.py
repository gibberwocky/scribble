#!/usr/bin/env python

import os, sys
os.environ['R_HOME'] = sys.exec_prefix+"/lib/R/"

from pathlib import Path
import re
import random
import scanpy as sc
import scFates as scf
import numpy as np
import pandas as pd
import decoupler as dc
import seaborn as sns
import omnipath as op
import harmonypy as hm
import gseapy as gp
import matplotlib.pyplot as plt
from statsmodels.stats.multitest import multipletests
from pygam import LinearGAM, s
from scribble.refine import restore_counts
from scribble.import_data import setup_environment

#def safe_sheet_name(name):
#    name = re.sub(r"[/\\?*\\[\\]:]", "_", name)
#    return name[:31]

def safe_sheet_name(name):
    name = str(name)
    name = re.sub(r'[\\/*?:\[\]]', '_', name)
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


def run_gam(
    adata_traj,
    adata_markers
):

    marker_tfs = adata_markers.uns["marker_tfs"]
    trajectory_tfs = (
        marker_tfs
        .query(
            "logfoldchanges > 0"
        )["names"]
        .drop_duplicates()
        .tolist()
    )

    results = []
    x = (adata_markers.obs["t"].values.reshape(-1, 1))
    for tf in trajectory_tfs:
        y = adata_markers[:, tf].X
        if hasattr(y, "toarray"):
            y = y.toarray().ravel()
        # skip near-constant genes
        if np.var(y) < 1e-6:
            continue
        try:
            gam = LinearGAM(
                s(0),
                max_iter=200
            ).fit(x, y)
            results.append({
                "TF": tf,
                "explained_deviance": gam.statistics_["pseudo_r2"]["explained_deviance"],
                "p": gam.statistics_["p_values"][1]
            })
        except Exception as e:
            print(f"{tf}: {e}")

    tf_gam = pd.DataFrame(results)
    tf_gam["fdr"] = multipletests(
        tf_gam["p"],
        method="fdr_bh"
    )[1]
    tf_gam = tf_gam.sort_values(
        ["fdr", "explained_deviance"],
        ascending=[True, False]
    )

    return(tf_gam)


def run_pseudotime(args):

    PROJECT_DIR = Path(args.project_dir)
    ADATA_DIR = PROJECT_DIR / "scribble/adata"

    TABLE_DIR = PROJECT_DIR / "scribble/tables/pseudotime"
    TABLE_DIR.mkdir(parents=True, exist_ok=True)

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

    print("Writing outputs")
    adata_markers.write(ADATA_DIR / f"{trajectory_name}_markers.h5ad")
    adata_traj.write(ADATA_DIR / f"{trajectory_name}_trajectory.h5ad")
    pbs.write(ADATA_DIR / f"{trajectory_name}_pseudobulk.h5ad")


    if args.association:

        print("Running trajectory association analysis")
        scf.tl.test_association(adata_traj, n_jobs=args.n_jobs)

        print("Running general additive model analysis: fitting TF expression ~ spline(t)")
        tf_gam = run_gam(adata_traj, adata_markers)

        # Update adata_markers with GAM results
        adata_markers.uns["tf_gam"] = tf_gam
        adata_markers.write(ADATA_DIR / f"{trajectory_name}_markers.h5ad")

        # Merge results back to markers
        marker_tfs = adata_markers.uns["marker_tfs"]
        tf_summary = (
            marker_tfs[
                ["group", "names", "scores", "logfoldchanges"]
            ]
            .merge(
                adata_markers.uns["tf_gam"],
                left_on="names",
                right_on="TF",
                how="left"
            )
        )

        # 1. Which TFs are associated with progression along pseudotime? (lineage program)
        tf_summary.query("fdr < 0.05") \
              .sort_values(
                  "explained_deviance",
                  ascending=False
              ).to_excel(
                  TABLE_DIR / f"{trajectory_name}_trajectory_tfs.xlsx",
                  index=False
              )

        # 2. Which TFs distinguish this state from the others? (branch TFs)
        branch_tfs = (
            tf_summary
            .query("logfoldchanges > 0")
            .sort_values(
                ["group", "scores"],
                ascending=[True, False]
            )
        )
        with pd.ExcelWriter(
            TABLE_DIR / f"{trajectory_name}_branch_tfs.xlsx",
            engine="openpyxl"
        ) as writer:
            for group, df in (
                branch_tfs
                .groupby("group", observed=True)
            ):
                (
                    df.sort_values(
                        "scores",
                        ascending=False
                    )
                    .to_excel(
                        writer,
                        sheet_name=safe_sheet_name(group),
                        index=False
                    )
                )

        # 3. TFs that are both branch markers and significantly associated with the trajectory
        candidate_regulators = (
            tf_summary
            .query(
                "logfoldchanges > 0 and fdr < 0.05"
            )
        )
        outfile = (TABLE_DIR /f"{trajectory_name}_candidate_regulators.xlsx")
        with pd.ExcelWriter(outfile, engine="openpyxl") as writer:
            for group, df in candidate_regulators.groupby(
                "group",
                observed=True
            ):
                (
                    df.sort_values(
                        ["scores", "explained_deviance"],
                        ascending=False
                    )
                    .to_excel(
                        writer,
                        sheet_name=safe_sheet_name(group),
                        index=False
                    )
                )



        # ============================================================================
        # Visualisation
        # ============================================================================

        ordered_celltypes = (
            adata_markers.obs
            .groupby(args.celltype)["t"]
            .median()
            .sort_values()
            .index
            .tolist()
        )

        # Plot the top genes per cell type
        markers = sc.get.rank_genes_groups_df(adata_markers, group=None)

        top_markers = (
            markers
            .query("logfoldchanges > 0")
            .sort_values(
                ["group", "scores"],
                ascending=[True, False]
            )
            .groupby("group", observed=True)
            .head(3)
        )

        genes = top_markers["names"].tolist()

        sc.pl.dotplot(
            adata_markers,
            var_names=genes,
            groupby=args.celltype,
            standard_scale="var"
        )

        plt.savefig(
            PLOT_DIR / f"dotp_{trajectory_name}_genes.png",
            dpi=300,
            bbox_inches="tight"
        )

        plt.close()


        # Plot the top genes per cell type
        markers = sc.get.rank_genes_groups_df(adata_markers, group=None)
        marker_tfs = markers[
            markers["names"].isin(tf_list)
        ].copy()

        top_tfs = (
            marker_tfs
            .query("logfoldchanges > 0")
            .sort_values(
                ["group", "scores"],
                ascending=[True, False]
            )
            .groupby("group", observed=True)
            .head(3)
        )

        tf_genes = top_tfs["names"].tolist()

        sc.pl.dotplot(
            adata_markers,
            var_names=tf_genes,
            groupby=args.celltype,
            standard_scale="var"
        )

        plt.savefig(
            PLOT_DIR / f"dotp_{trajectory_name}_TF_genes.png",
            dpi=300,
            bbox_inches="tight"
        )

        plt.close()


        # Plot top active TFs specific to cell types, aggregated across samples
        ulm = pbs.obsm["score_ulm"].copy()

        ulm.index = pbs.obs_names
        ulm["cell_type"] = pbs.obs[args.celltype].values

        ct_means = (
            ulm
            .groupby("cell_type", observed=True)
            .mean()
        )

        specificity = pd.DataFrame(
            index=ct_means.index,
            columns=ct_means.columns
        )

        for ct in ct_means.index:
            specificity.loc[ct] = (
                ct_means.loc[ct]
                -
                ct_means.drop(ct).mean()
            )

        top_activity_tfs = (
            specificity
            .apply(
                lambda row:
                row.sort_values(
                    ascending=False
                ).head(3).index.tolist(),
                axis=1
            )
        )

        top_genes = sorted(
            set(
                tf
                for tfs in top_activity_tfs
                for tf in tfs
            )
        )

        ct_means = ct_means.loc[ordered_celltypes]

        sns.clustermap(
            ct_means[top_genes],
            cmap="RdBu_r",
            center=0,
            z_score=1,
            row_cluster=False,
        )

        plt.savefig(
            PLOT_DIR / f"ulm_{trajectory_name}_TF_ct_means.png",
            dpi=300,
            bbox_inches="tight"
        )

        plt.close()

        # Heatmap of top TFs ordered by pseudotime (cells and TFs)
        top_tfs = (
            candidate_regulators
            .sort_values(
                "scores",
                ascending=False
            )
            .groupby("group", observed=True)
            .head(3)
        )

        genes = (
            top_tfs["names"]
            .dropna()
            .unique()
            .tolist()
        )

        # Estimate peak pseudotime for each TF
        tf_peak_time = {}
        t = adata_markers.obs["t"].values
        for tf in genes:
            expr = adata_markers[:, tf].X
            if hasattr(expr, "toarray"):
                expr = expr.toarray().ravel()
            expr = np.maximum(expr, 0)
            if expr.sum() > 0:
                tf_peak_time[tf] = np.average(
                    t,
                    weights=expr
                )
            else:
                tf_peak_time[tf] = np.nan

        # Order TFs along this trajectory
        ordered_genes = [
            k for k, v in
            sorted(
                tf_peak_time.items(),
                key=lambda x: x[1]
            )
        ]

        # Plot
        adata_markers.obs[args.celltype] = pd.Categorical(
            adata_markers.obs[args.celltype],
            categories=ordered_celltypes,
            ordered=True
        )

        sc.pl.matrixplot(
            adata_markers,
            var_names=ordered_genes,
            groupby=args.celltype,
            standard_scale="var",
            cmap="viridis"
        )

        plt.savefig(
            PLOT_DIR / f"hmap_{trajectory_name}_TFs.png",
            dpi=300,
            bbox_inches="tight"
        )

        plt.close()


        # Dot plot the same TFs
        sc.pl.dotplot(
            adata_markers,
            var_names=genes,
            groupby=args.celltype,
            standard_scale="var"
        )

        plt.savefig(
            PLOT_DIR / f"dotp_{trajectory_name}_TFs.png",
            dpi=300,
            bbox_inches="tight"
        )

        plt.close()

        # UMAP of these markers
        for group, df in top_tfs.groupby("group", observed=True):

            sc.pl.umap(
                adata_markers,
                color=df["names"].tolist(),
                cmap="Reds",
                colorbar_loc=None,
                ncols=3,
                title=df["names"].tolist()
            )

            plt.savefig(
                PLOT_DIR / f"umap_{safe_sheet_name(group)}_TFs.png",
                dpi=300,
                bbox_inches="tight"
            )

            plt.close()

        # ============================================================================
        # Marker results
        # ============================================================================

        markers = sc.get.rank_genes_groups_df(
            adata_markers,
            group=None
        )

        markers = markers.rename(
            columns={
                "names": "gene",
                "scores": "marker_score",
                "logfoldchanges": "marker_logFC",
                "pvals": "marker_p",
                "pvals_adj": "marker_fdr"
            }
        )

        markers["expr_sig"] = markers["marker_fdr"] < 0.05


        # ============================================================================
        # GAM results
        # ============================================================================

        gam = (
            adata_markers.uns["tf_gam"]
            .rename(columns={"TF": "gene"})
            .copy()
        )

        gam["gam_sig"] = gam["fdr"] < 0.05

        gam = gam.rename(
            columns={
                "explained_deviance": "gam_explained_deviance",
                "p": "gam_p",
                "fdr": "gam_fdr"
            }
        )


        # ============================================================================
        # scFates association results
        # ============================================================================

        assoc = (
            adata_traj.var[
                ["A", "p_val", "fdr", "signi"]
            ]
            .copy()
        )

        assoc["gene"] = assoc.index

        assoc = assoc.rename(
            columns={
                "A": "scfates_A",
                "p_val": "scfates_p",
                "fdr": "scfates_fdr",
                "signi": "scfates_sig"
            }
        )


        # ============================================================================
        # ULM activity
        # ============================================================================

        ulm_score = pbs.obsm["score_ulm"].copy()
        ulm_padj = pbs.obsm["padj_ulm"].copy()

        ulm_long = []

        for cell_type in ulm_score.index:
            tmp = pd.DataFrame({
                "group": cell_type,
                "gene": ulm_score.columns,
                "ulm_score": ulm_score.loc[cell_type].values,
                "ulm_fdr": ulm_padj.loc[cell_type].values
            })
            ulm_long.append(tmp)


        ulm_long = pd.concat(ulm_long, ignore_index=True)
        ulm_long["ulm_sig"] = ulm_long["ulm_fdr"] < 0.05
        ulm_long["group"] = (
            ulm_long["group"]
            .str.replace(r"^[^_]+_", "", regex=True)
        )


        # ============================================================================
        # Merge everything
        # ============================================================================

        summary = (
            markers
            .merge(gam, on="gene", how="left")
            .merge(assoc, on="gene", how="left")
            .merge(
                ulm_long,
                on=["group", "gene"],
                how="left"
            )
        )


        # ============================================================================
        # DoRothEA coverage
        # ============================================================================

        regulon_size = (
            net.groupby("source")
               .size()
               .rename("regulon_size")
               .reset_index()
               .rename(columns={"source": "gene"})
        )

        summary = summary.merge(
            regulon_size,
            on="gene",
            how="left"
        )

        summary["regulon_size"] = (
            summary["regulon_size"]
            .fillna(0)
            .astype(int)
        )


        # ============================================================================
        # Testing availability
        # ============================================================================

        summary["scfates_tested"] = (
            summary["gene"]
            .isin(adata_traj.var_names)
        )

        summary["ulm_tested"] = (
            summary["regulon_size"] > 0
        )


        # ============================================================================
        # Significance flags
        # ============================================================================

        summary["expr_sig"] = summary["expr_sig"].fillna(False)
        summary["gam_sig"] = summary["gam_sig"].fillna(False)

        summary["scfates_sig"] = (
            summary["scfates_sig"]
            .fillna(False)
        )

        summary["ulm_sig"] = (
            summary["ulm_sig"]
            .fillna(False)
        )

        for col in [
            "gam_sig",
            "scfates_sig",
            "ulm_sig"
        ]:
            summary[col] = (
                summary[col]
                .fillna(False)
                .astype(bool)
            )

        # ============================================================================
        # Available tests
        # ============================================================================

        # expr + gam are always available
        summary["tests_available"] = (
            2
            + summary["scfates_tested"].astype(int)
            + summary["ulm_tested"].astype(int)
        )


        # ============================================================================
        # Significant tests
        # ============================================================================

        summary["tests_significant"] = (
            summary["expr_sig"].astype(int)
            + summary["gam_sig"].astype(int)
            + (
                summary["scfates_sig"]
                & summary["scfates_tested"]
            ).astype(int)
            + (
                summary["ulm_sig"]
                & summary["ulm_tested"]
            ).astype(int)
        )

        # ============================================================================
        # Support fraction
        # ============================================================================

        summary["support_fraction"] = (
            summary["tests_significant"]
            / summary["tests_available"]
        )

        summary["evidence_level"] = np.select(
            [
                summary["support_fraction"] >= 1.0,
                summary["support_fraction"] >= 0.75,
                summary["support_fraction"] >= 0.5
            ],
            [
                "Very High",
                "High",
                "Moderate"
            ],
            default="Low"
        )

        summary["evidence_pattern"] = (
            summary["expr_sig"].map({True:"E", False:"-"})
            + summary["gam_sig"].map({True:"G", False:"-"})
            + summary["scfates_sig"].fillna(False).map({True:"S", False:"-"})
            + summary["ulm_sig"].fillna(False).map({True:"U", False:"-"})
        )

        # EGSU = supported by all four analyses
        # EGS- = expression + GAM + scFates
        # EG-- = expression + GAM only
        # E--U = marker + ULM activity

        summary.loc[
            summary["gene"] == "RFX3",
            ["group", "gene", "marker_score", "marker_logFC", "marker_fdr"]
        ]


        # Plot the top scFates-associated genes per cell type
        top_scfates_tfs = (
            summary.loc[
                summary["scfates_tested"]
                & summary["scfates_sig"]
                & summary["gene"].isin(tf_list)
            ]
            .sort_values(
                "scfates_A",
                ascending=False
            )
            .drop_duplicates("gene")
            .head(12)
        )

        genes = top_scfates_tfs["gene"].tolist()

        sc.pl.dotplot(
            adata_markers,
            var_names=genes,
            groupby=args.celltype,
            standard_scale="var"
        )

        plt.savefig(
            PLOT_DIR / f"dotp_{trajectory_name}_scFates-genes.png",
            dpi=300,
            bbox_inches="tight"
        )

        plt.close()


        # ============================================================================
        # Save
        # ============================================================================

        outfile = (TABLE_DIR / f"{trajectory_name}_integrated_candidate_regulators.xlsx")
        with pd.ExcelWriter(outfile, engine="openpyxl") as writer:
            for group, df in summary.groupby(
                "group",
                observed=True
            ):
                (
                    df.query(
                        "marker_logFC > 0 and marker_fdr < 0.05"
                    )
                    .sort_values(
                        [
                            "support_fraction",
                            "gam_explained_deviance",
                            "marker_score"
                        ],
                        ascending=[False, False, False]
                    )
                    .to_excel(
                        writer,
                        sheet_name=safe_sheet_name(group),
                        index=False
                    )
                )



        # --------------------------------------------------------------------------
        # Activated regulators
        # --------------------------------------------------------------------------

        # Requires the cell-type specificty ULM TFs generated earlier (repeated here)
        ulm = pbs.obsm["score_ulm"].copy()
        ulm.index = pbs.obs_names
        ulm["cell_type"] = pbs.obs[args.celltype].values
        ct_means = (
            ulm
            .groupby("cell_type", observed=True)
            .mean()
        )
        specificity = pd.DataFrame(
            index=ct_means.index,
            columns=ct_means.columns
        )
        for ct in ct_means.index:
            specificity.loc[ct] = (
                ct_means.loc[ct]
                -
                ct_means.drop(ct).mean()
            )


        specificity_long = (
            specificity
            .stack()
            .reset_index()
        )

        specificity_long.columns = [
            "group",
            "gene",
            "ulm_specificity"
        ]

        summary = summary.merge(
            specificity_long,
            on=["group", "gene"],
            how="left"
        )

        outfile = (TABLE_DIR / f"{trajectory_name}_integrated_ulm_activated_regulators.xlsx")
        with pd.ExcelWriter(outfile, engine="openpyxl") as writer:
            for group, df in summary.groupby(
                "group",
                observed=True
            ):
                (
                    df.query(
                        "ulm_tested and ulm_sig and (gam_sig or expr_sig) and ulm_specificity > 0"
                    )
                    .sort_values(
                        [
                            "support_fraction",
                            "ulm_fdr",
                            "gam_fdr"
                        ],
                        ascending=[False, True, True]
                    )
                    .to_excel(
                        writer,
                        sheet_name=safe_sheet_name(group),
                        index=False
                    )
                )

        # --------------------------------------------------------------------------
        # Activated downstream targets per group
        # --------------------------------------------------------------------------

        # Define activated TFs
        activated_tfs = summary.loc[
            summary["ulm_tested"]
            & summary["ulm_sig"]
            & (summary["ulm_specificity"] > 0)
        ].copy()

        # Marker genes by group
        marker_genes = {
            group: set(
                df.loc[
                    (df["marker_logFC"] > 0)
                    & (df["marker_fdr"] < 0.05),
                    "gene"
                ]
            )
            for group, df in summary.groupby(
                "group",
                observed=True
            )
        }

        # Build activated target summary
        rows = []
        group_active_targets = {}
        for _, row in activated_tfs.iterrows():
            group = row["group"]
            tf = row["gene"]
            # DoRothEA targets for TF
            tf_targets = net.loc[
                net["source"] == tf,
                "target"
            ].unique()
            n_dorothea_targets = len(tf_targets)
            # Overlap with significant positive markers
            active_targets = sorted(
                set(tf_targets)
                & marker_genes[group]
            )
            # Save for later enrichment
            if group not in group_active_targets:
                group_active_targets[group] = set()
            group_active_targets[group].update(active_targets)
            rows.append({
                **row.to_dict(),
                "n_dorothea_targets":
                    n_dorothea_targets,
                "n_active_targets":
                    len(active_targets),
                "active_targets":
                    "; ".join(active_targets),
                "target_fraction":
                    (
                        len(active_targets)
                        / n_dorothea_targets
                    )
                    if n_dorothea_targets > 0
                    else np.nan,
                "target_score":
                    (
                        len(active_targets)
                        *
                        (
                            len(active_targets)
                            / n_dorothea_targets
                        )
                    )
                    if n_dorothea_targets > 0
                    else 0
            })

        activated_tf_summary = pd.DataFrame(rows)

        outfile = (TABLE_DIR / f"{trajectory_name}_integrated_ulm_activated_targets.xlsx")
        with pd.ExcelWriter(outfile, engine="openpyxl") as writer:
            for group, target_df in activated_tf_summary.groupby(
                "group",
                observed=True
            ):
                (
                    target_df
                    .sort_values(
                        [
                            "target_score",
                            "n_active_targets",
                            "ulm_specificity"
                        ],
                        ascending=False
                    )
                    .to_excel(
                        writer,
                        sheet_name=safe_sheet_name(group),
                        index=False
                    )
                )

        # GSEA for each cell type
        gene_sets = [
            "GO_Biological_Process_2026",
            "Reactome_Pathways_2024"
        ]

        for group, active_targets in group_active_targets.items():
            active_targets = sorted(active_targets)
            if len(active_targets) < 10:
                print(
                    f"Skipping {group}: "
                    f"only {len(active_targets)} targets"
                )
                continue
            print(
                f"Running GSEA for {group}: "
                f"{len(active_targets)} targets"
            )
            try:
                enr = gp.enrichr(
                    gene_list=active_targets,
                    gene_sets=gene_sets,
                    organism="human",
                    outdir=None
                )
                results = (
                    enr.results
                    .sort_values(
                        "Adjusted P-value"
                    )
                )
                # Save enrichment table
                results.to_excel(
                    TABLE_DIR / f"{safe_sheet_name(group)}_gsea.xlsx",
                    index=False
                )
                # Dotplot
                ax = gp.dotplot(
                    results,
                    column="Adjusted P-value",
                    top_term=15,
                    title=f"Activated programmes in {group}"
                )
                plt.savefig(
                    PLOT_DIR / f"{safe_sheet_name(group)}_gsea_dotplot.png",
                    dpi=300,
                    bbox_inches="tight"
                )
                plt.close()
            except Exception as e:
                print(
                    f"GSEA failed for {group}: "
                    f"{e}"
                )
