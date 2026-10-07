#!/usr/bin/env python
"""
Perform DE and GSEA analyses.

This CLI supports arbitrary straightforward model adjustment formulas but simple three-part categorical contrasts only.
Interaction-specific hypotheses require explicit additional implementation.

Returns:
  <celltype>_DE_GSEA_summary.csv                      Summary of DE and GSEA analyses
                                                      - DE: n_genes_tested, n_FDR_005, n_FDR_005_up, n_FDR_005_down, n_FDR_005_log2FC2, n_nominal
                                                      - GSEA: gsea_duplicate_rank_pct, Hallmark_sig_pathways, Hallmark_sig_up, Hallmark_sig_down,
                                                              KEGG_sig_pathways, KEGG_sig_up, KEGG_sig_down,
                                                              GO_BP_sig_pathways, GO_BP_sig_up, GO_BP_sig_down,
                                                              Reactome_sig_pathways, Reactome_sig_up, Reactome_sig_down,
                                                              Wiki_sig_pathways,Wiki_sig_up,Wiki_sig_down
                                                      - QC:   qc_status, qc_flags, n_samples_qc, n_genes_qc,
                                                              library_size_ratio, size_factor_ratio, size_factor_library_corr, cell_count_library_corr,
                                                              PC1_variance_pct, PC2_variance_pct,
                                                      - Fit:  n_dispersion_outliers, pct_dispersion_outliers,
                                                              n_cooks_outlier_genes, pct_cooks_outlier_genes,
                                                              n_lfc_not_converged, pct_lfc_not_converged,
                                                              n_MAP_not_converged, pct_MAP_not_converged,
                                                              n_genewise_not_converged, pct_genewise_not_converged,
                                                      - Cook: worst_cooks_sample, worst_sample_n_cooks_gt1, worst_cooks_q99_sample,
                                                              max_cooks_q99, max_cooks_sample, max_cooks_value
  DE/<annotation>/
      *_all_genes.csv                                 All genes that pass prefilter and tested by PyDESeq2: (pb >= 10).sum(axis=0) >= 3
      *_nominal.csv                                   Genes with pval < 0.05 and abs(log2FoldChange) > 0.5
      *_significant.csv                               Genes with FDR < 0.05
      *_significant_log2FC1.csv                       Genes with FDR < 0.05 and abs(log2FoldChange) > 1
  GSEA/<annotation>/
      *cell_ranking.csv                               PyDESeq2 Wald statistics for TEST-versus-REFERENCE contrast unfiltered for significance, rankings are input for gseapy
                                                      - large +ve is strong evidence gene increased in TEST relative to REFERENCE
                                                      - Near zero indicates little evidence for treatment-associated change
                                                      - large -ve is strong evidence gene decreased in TEST relative to REFRENCE
      <celltype>/
          <geneset>/
              all_pathways.csv                        GSEA result table listing every pathway that passed GSEA size filters (min_size=10, max_size=500)
              gene_sets.gmt                           GSEApy gene set Gene Matrix Transposed definition file (list of genes per pathway)
              gseapy.gene_set.prerank.report.csv      GSEApy-generated equiavalent to all_pathways.csv
              significant...down.csv                  Pathways that are significantly depleted (FDR < 0.05 and NES < 0)
              significant...up.csv                    Pathways that are significantly enriched (FDR < 0.05 and NES > 0)
              top25...down.csv                        Top 25 pathways that are depleted (NES < 0, sorted by FDR)
              top25...up.csv                          Top 25 pathways that are enriched (NES > 0, sorted by FDR)
              <preank>/
                  <pathway>.pdf                       GSEA enrichment plot
                                                          Bottom panel:   left/red/positive indicates genes increased in TEST level
                                                                          right/blue/negative indicates genes decreased in TEST level (or realtively higher in REFERENCE level)
                                                          Middle panel:   small vertical black lines represent genes belonging to the pathway
                                                          Top panel:      running enrichment score
  inspection/<annotation>/
      <celltype>/
          cooks_by_sample.png                         Indicates how influential a sample's count data is on the fitted DE model, samples should be similar
          diagnostics.txt                             Descriptive QC summary for PyDESeq2 fit
          dispersion.png                              Shows count variability in PyDESeq2 estimates for each gene relative to its average expression
                                                      - gene-wise is the initial dispersion estimate calcualted independently for each gene
                                                      - final is the dispersion estimate used by the model after PyDESeq2's dispersion trend/shrinkage procedure
          pca_donor.png                               Pseudobulk PCA coloured by donor
          pca_treatment.png                           Pseudobulk PCA coloured by treatment
          pca.csv                                     Sample PC1, PC2, and metadata for generating PCA plots
          sample_correlation.csv                      Pairwise sample Pearson correlation matrix acorss 1000 log-normalised HVGs used for PCA (above)
          sample_correlation.png                      Heatmap of sample correlation matrix
          sample_qc.csv                               Sample-level QC table
          size_factors.png                            Top-panel pseudobulk library size, bottom-panel PyDESeq2 size factors where dashed line indicates size factor 1
"""

import argparse
import traceback
from pathlib import Path
import scanpy as sc
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import re
import seaborn as sns
from scribble.refine import restore_counts
from scipy import sparse
import gseapy as gp
from pydeseq2.dds import DeseqDataSet
from pydeseq2.ds import DeseqStats
from pydeseq2.default_inference import DefaultInference
import json

def safe_filename(value):
    """Convert a string into a filesystem-safe filename component."""
    value = str(value).strip()
    value = re.sub(r"[^\w.-]+", "_", value)
    value = re.sub(r"_+", "_", value)
    return value.strip("_.")

def get_formula_variables(formula):
    """
    Extract simple variable names from a Formulaic-style
    design formula.

    Intended for validation of formulas such as:
      ~ donor + treatment
      ~ treatment
      ~ donor + treatment + batchInfo

    Formulaic itself remains responsible for actual parsing.
    """

    tokens = re.findall(
        r"\b[A-Za-z_]\w*\b",
        formula
    )

    reserved = {
        "C",
        "I",
    }

    return {
        token
        for token in tokens
        if token not in reserved
    }

# --------------------------------------------------
# CLI parameters
# --------------------------------------------------
parser = argparse.ArgumentParser()
parser.add_argument("--project_dir", type=str, required=True)
parser.add_argument("--input_file", type=str, required=True)
parser.add_argument("--annotation", type=str, required=True)
parser.add_argument("--analysis_name", type=str, required=True,
    help=("Name identifying this analysis, e.g. 'combined', 'organoid', 'KOLF2', 'ex_vivo'."),
)
parser.add_argument("--samples", nargs="+", default=None,
    help=("Optional list of sample IDs to include. "
        "If omitted, all samples containing the contrast "
        "levels are considered.")
)
parser.add_argument("--formula", type=str, required=True,
    help=("PyDESeq2 design formula, e.g. '~ donor + treatment' or '~ treatment'.")
)
parser.add_argument("--contrast", nargs=3, required=True,
    metavar=("VARIABLE", "TEST_LEVEL", "REFERENCE_LEVEL",),
    help=("PyDESeq2 contrast. Example: --contrast treatment IL17A veh")
)
parser.add_argument("--min_gene_count", type=int, default=10)
parser.add_argument("--min_gene_samples", type=int, default=3)
parser.add_argument("--min_cells_per_sample", type=int, default=20)
parser.add_argument("--min_samples_per_group", type=int, default=2)
parser.add_argument("--min_informative_donors", type=int, default=3)
parser.add_argument("--require_informative_donors", action=argparse.BooleanOptionalAction, default=True,
    help=("Require min_informative_donors donors "
        "represented in both contrast levels. "
        "Use --no-require_informative_donors for "
        "an explicitly unblocked analysis.")
)
parser.add_argument("--gsea_permutations", type=int, default=10000)
parser.add_argument("--n_cpus", type=int, default=8)
args = parser.parse_args()

CONTRAST_VARIABLE = args.contrast[0]
CONTRAST_TEST = args.contrast[1]
CONTRAST_REFERENCE = args.contrast[2]
FORMULA_VARIABLES = get_formula_variables(args.formula)
DONOR_IN_MODEL = ("donor" in FORMULA_VARIABLES)
REQUIRED_MODEL_COLUMNS = (FORMULA_VARIABLES | {CONTRAST_VARIABLE})
if (
    CONTRAST_VARIABLE
    not in FORMULA_VARIABLES
):
    raise ValueError(
        f"Contrast variable "
        f"'{CONTRAST_VARIABLE}' "
        "is not present in design formula "
        f"'{args.formula}'."
    )

RUN_INFO = {
    "analysis_name": args.analysis_name,
    "formula": args.formula,
    "contrast_variable": (
        CONTRAST_VARIABLE
    ),
    "contrast_test": (
        CONTRAST_TEST
    ),
    "contrast_reference": (
        CONTRAST_REFERENCE
    ),
}

# Create directory structure
ANALYSIS_NAME = safe_filename(args.analysis_name)
OUTDIR = (Path(args.project_dir) / "scribble" / "DGE" / ANALYSIS_NAME)
DE_DIR = (OUTDIR / "DE" / args.annotation)
GSEA_DIR = (OUTDIR / "GSEA" / args.annotation)
PB_DIR = (OUTDIR / "pseudobulk" / args.annotation)
INSPECTION_DIR = (OUTDIR / "inspection" / args.annotation)
DE_DIR.mkdir(exist_ok=True, parents=True)
GSEA_DIR.mkdir(exist_ok=True, parents=True)
PB_DIR.mkdir(exist_ok=True, parents=True)
INSPECTION_DIR.mkdir(parents=True, exist_ok=True)


gene_sets = {
    "Hallmark": "MSigDB_Hallmark_2020",
    "KEGG": "KEGG_2026",
    "GO_BP": "GO_Biological_Process_2026",
    "Reactome": "Reactome_Pathways_2024",
    "Wiki": "WikiPathways_2024_Human"
}


# Run manifest
run_manifest = {
    "analysis_name": args.analysis_name,
    "input_file": str(
        Path(args.input_file)
        .resolve()
    ),
    "annotation": args.annotation,
    "samples": args.samples,
    "formula": args.formula,
    "formula_variables": sorted(FORMULA_VARIABLES),
    "contrast_variable": (CONTRAST_VARIABLE),
    "contrast_test": (CONTRAST_TEST),
    "contrast_reference": (CONTRAST_REFERENCE),
    "require_informative_donors": (args.require_informative_donors),
    "min_gene_count": (args.min_gene_count),
    "min_gene_samples": (args.min_gene_samples),
    "min_cells_per_sample": (args.min_cells_per_sample),
    "min_samples_per_group": (args.min_samples_per_group),
    "min_informative_donors": (args.min_informative_donors),
    "gsea_permutations": (args.gsea_permutations),
    "n_cpus": args.n_cpus,
    "gene_sets": gene_sets,
}

with open(
    OUTDIR /
    "run_manifest.json",
    "w",
) as handle:

    json.dump(
        run_manifest,
        handle,
        indent=2,
    )

# --------------------------------------------------
# helpers
# --------------------------------------------------

def pseudobulk_celltype(
    adata,
    cell_type,
):
    """
    Aggregate raw counts within each sample for one cell type.

    Returns
    -------
    pb : DataFrame
        samples x genes integer pseudobulk counts

    meta_pb : DataFrame
        sample-level metadata aligned exactly with pb
    """
    subset = adata[
        adata.obs[args.annotation] == cell_type
    ]
    # ----------------------------------------------
    # nuclei per sample
    # ----------------------------------------------
    cell_counts = (
        subset.obs["sample"]
        .value_counts()
    )
    keep_samples = cell_counts[
        cell_counts >= args.min_cells_per_sample
    ].index
    subset = subset[
        subset.obs["sample"].isin(keep_samples)
    ]
    if subset.n_obs == 0:
        return None
    pseudobulk = []
    sample_meta = []
    # Sort to make outputs deterministic
    keep_samples = sorted(keep_samples)
    for sample in keep_samples:
        sample_adata = subset[
            subset.obs["sample"] == sample
        ]
        X = sample_adata.layers["counts"]
        if sparse.issparse(X):
            summed = np.asarray(
                X.sum(axis=0)
            ).ravel()
        else:
            summed = np.asarray(
                X.sum(axis=0)
            ).ravel()
        pseudobulk.append(summed)

        meta_row = {
            "sample": sample,
            "n_cells": sample_adata.n_obs,
            "library_size": summed.sum(),
        }

        metadata_to_retain = (
            REQUIRED_MODEL_COLUMNS
            |
            {
                "donor",
                "cell_line",
                "weeks",
                "sex",
                "age",
                "batchInfo",
                "treatment",
            }
        )

        for col in sorted(metadata_to_retain):

            if col in sample_adata.obs.columns:

                values = (
                    sample_adata.obs[col]
                    .dropna()
                    .unique()
                )

                if len(values) > 1:
                    raise ValueError(
                        f"{sample}: inconsistent {col} "
                        "within sample"
                    )

                if len(values) == 1:
                    meta_row[col] = values[0]

        sample_meta.append(meta_row)

    pb = pd.DataFrame(
        pseudobulk,
        columns=subset.var_names,
        index=keep_samples,
    )
    # Counts should be integer-like and DESeq2 expects integers.
    pb = pb.round().astype(np.int64)
    meta_pb = (
        pd.DataFrame(sample_meta)
        .set_index("sample")
        .loc[pb.index]
    )
    return pb, meta_pb


def assess_design(meta_pb):
    """
    Assess whether sufficient replication remains after
    cell-count filtering for the requested contrast.

    If donor metadata is present, also report donors represented
    under both contrast levels.
    """

    variable = CONTRAST_VARIABLE
    test_level = CONTRAST_TEST
    reference_level = CONTRAST_REFERENCE

    if variable not in meta_pb.columns:
        raise ValueError(
            f"Contrast variable '{variable}' "
            "is not present in pseudobulk metadata."
        )

    # ----------------------------------------------
    # Samples per contrast level
    # ----------------------------------------------

    groups = (
        meta_pb[variable]
        .astype(str)
    )

    n_test = int(
        groups.eq(
            test_level
        ).sum()
    )

    n_reference = int(
        groups.eq(
            reference_level
        ).sum()
    )

    result = {
        "n_samples": len(meta_pb),
        f"n_{reference_level}": (
            n_reference
        ),
        f"n_{test_level}": (
            n_test
        ),
    }

    # ----------------------------------------------
    # Donor-aware design assessment
    # ----------------------------------------------

    if "donor" in meta_pb.columns:

        donor_group = pd.crosstab(
            meta_pb["donor"],
            groups,
        )

        reference_counts = (
            donor_group[
                reference_level
            ]
            if reference_level
            in donor_group.columns
            else pd.Series(
                0,
                index=donor_group.index,
            )
        )

        test_counts = (
            donor_group[
                test_level
            ]
            if test_level
            in donor_group.columns
            else pd.Series(
                0,
                index=donor_group.index,
            )
        )

        informative = (
            donor_group.index[
                (reference_counts > 0)
                &
                (test_counts > 0)
            ]
            .tolist()
        )

        n_informative = len(
            informative
        )

        paired_1to1 = int(
            (
                (reference_counts == 1)
                &
                (test_counts == 1)
            )
            .sum()
        )

        replicated = (
            donor_group.index[
                (
                    (
                        reference_counts > 1
                    )
                    |
                    (
                        test_counts > 1
                    )
                )
                &
                (reference_counts > 0)
                &
                (test_counts > 0)
            ]
            .tolist()
        )

        result.update(
            {
                "n_donors": int(
                    meta_pb[
                        "donor"
                    ]
                    .nunique()
                ),
                "n_informative_donors": (
                    n_informative
                ),
                "informative_donors": (
                    "|".join(
                        map(
                            str,
                            informative,
                        )
                    )
                ),
                "n_1to1_paired_donors": (
                    paired_1to1
                ),
                "replicated_donors": (
                    "|".join(
                        map(
                            str,
                            replicated,
                        )
                    )
                ),
            }
        )

        if args.require_informative_donors:

            if not DONOR_IN_MODEL:
                raise ValueError(
                    "--require_informative_donors was enabled, "
                    "but 'donor' is not present in the model formula. "
                    "Either include donor in the formula or use "
                    "--no-require_informative_donors."
                )

            enough_design = (
                n_informative
                >= args.min_informative_donors
            )

        else:

            enough_design = True

    else:

        if args.require_informative_donors:
            raise ValueError(
                "--require_informative_donors was enabled, "
                "but donor metadata is unavailable."
            )

        enough_design = True

    result["eligible"] = (
        (
            n_reference
            >=
            args.min_samples_per_group
        )
        and
        (
            n_test
            >=
            args.min_samples_per_group
        )
        and
        enough_design
    )

    return result


def filter_genes(pb):
    """
    Mild prefilter removing essentially uninformative genes.

    Retain genes with at least args.min_gene_count counts
    in at least args.min_gene_samples retained pseudobulk
    samples.
    """
    keep = (
        (pb >= args.min_gene_count)
        .sum(axis=0)
        >= args.min_gene_samples
    )
    return pb.loc[:, keep]

def run_pydeseq2(
    pb,
    meta_pb,
):
    """
    Run PyDESeq2 using the user-supplied design formula
    and contrast.
    """

    metadata = meta_pb.copy()

    # ----------------------------------------------
    # Defensive alignment check
    # ----------------------------------------------

    if not pb.index.equals(
        metadata.index
    ):
        raise ValueError(
            "Pseudobulk counts and metadata "
            "are not aligned."
        )

    # ----------------------------------------------
    # Ensure contrast variable exists
    # ----------------------------------------------

    if (
        CONTRAST_VARIABLE
        not in metadata.columns
    ):
        raise ValueError(
            f"Contrast variable "
            f"'{CONTRAST_VARIABLE}' "
            "not present in metadata."
        )


    # Ensure contrast variable is categorical
    metadata[CONTRAST_VARIABLE] = (
        metadata[CONTRAST_VARIABLE]
        .astype(str)
    )

    # Ensure contrast levels exist
    contrast_levels = set(
        metadata[
            CONTRAST_VARIABLE
        ]
        .astype(str)
    )

    required_levels = {
        CONTRAST_TEST,
        CONTRAST_REFERENCE,
    }

    missing_levels = (
        required_levels
        -
        contrast_levels
    )

    if missing_levels:
        raise ValueError(
            "Missing contrast levels: "
            f"{sorted(missing_levels)}"
        )

    missing_model_columns = (
        REQUIRED_MODEL_COLUMNS
        -
        set(metadata.columns)
    )

    if missing_model_columns:
        raise ValueError(
            "Required model metadata missing from "
            "pseudobulk metadata: "
            f"{sorted(missing_model_columns)}"
        )

    missing_model_values = (
        metadata[
            sorted(REQUIRED_MODEL_COLUMNS)
        ]
        .isna()
        .any()
    )

    bad_columns = (
        missing_model_values[
            missing_model_values
        ]
        .index
        .tolist()
    )

    if bad_columns:
        raise ValueError(
            "Missing values present in required model "
            "metadata columns: "
            f"{bad_columns}"
        )


    # ----------------------------------------------
    # PyDESeq2
    # ----------------------------------------------

    inference = DefaultInference(
        n_cpus=args.n_cpus
    )

    dds = DeseqDataSet(
        counts=pb,
        metadata=metadata,
        design=args.formula,
        refit_cooks=True,
        inference=inference,
        quiet=True,
    )

    design_matrix = np.asarray(
        dds.obsm[
            "design_matrix"
        ]
    )

    rank = np.linalg.matrix_rank(
        design_matrix
    )

    if rank < design_matrix.shape[1]:
        raise ValueError(
            "Design matrix is rank-deficient: "
            f"rank={rank}, "
            f"columns={design_matrix.shape[1]}. "
            f"Formula: {args.formula}"
        )

    dds.deseq2()

    # ----------------------------------------------
    # Explicit contrast
    # ----------------------------------------------

    ds = DeseqStats(
        dds,
        contrast=[
            CONTRAST_VARIABLE,
            CONTRAST_TEST,
            CONTRAST_REFERENCE,
        ],
        inference=inference,
        quiet=True,
    )

    ds.summary()

    de = ds.results_df.copy()

    de.index.name = "gene"

    de = (
        de
        .reset_index()
        .rename(
            columns={
                "pvalue": "pval",
                "padj": "FDR",
                "stat": "Wald_stat",
            }
        )
    )

    return de, dds, ds

def make_gsea_ranking(de):
    ranking = (
        de[
            [
                "gene",
                "Wald_stat",
            ]
        ]
        .rename(
            columns={
                "Wald_stat": "rank"
            }
        )
        .replace(
            [np.inf, -np.inf],
            np.nan
        )
        .dropna()
        .drop_duplicates(
            subset="gene"
        )
        .sort_values(
            ["rank", "gene"],
            ascending=[False, True],
        )
    )

    duplicate_rank_pct = (
        ranking["rank"]
        .duplicated(
            keep=False
        )
        .mean()
        * 100
    )

    return ranking, duplicate_rank_pct




# ------------------------------------------------------------------
# Data preparation
# ------------------------------------------------------------------

def prepare_inspection_data(dds, meta_pb):
    """Extract matrices and sample metadata in the same order as dds."""

    sample_names = list(dds.obs_names)
    gene_names = list(dds.var_names)

    normed = pd.DataFrame(
        np.asarray(dds.layers["normed_counts"]),
        index=sample_names,
        columns=gene_names,
    )

    cooks = pd.DataFrame(
        np.asarray(dds.layers["cooks"]),
        index=sample_names,
        columns=gene_names,
    )

    sample_info = meta_pb.loc[sample_names].copy()

    sample_info["size_factor"] = np.asarray(
        dds.obs["size_factors"],
        dtype=float,
    )

    if "library_size" not in sample_info.columns:
        sample_info["library_size"] = np.nan

    return sample_names, normed, cooks, sample_info


# ------------------------------------------------------------------
# Sample-level QC
# ------------------------------------------------------------------

def calculate_sample_qc(sample_info, normed, cooks):
    """Add descriptive sample-level QC metrics."""

    sample_info["n_detected_genes"] = (
        (normed > 0).sum(axis=1)
    )

    sample_info["max_cooks"] = cooks.max(axis=1)
    sample_info["median_cooks"] = cooks.median(axis=1)
    sample_info["cooks_q95"] = cooks.quantile(0.95, axis=1)
    sample_info["cooks_q99"] = cooks.quantile(0.99, axis=1)
    sample_info["n_cooks_gt_1"] = (cooks > 1).sum(axis=1)
    sample_info["n_cooks_gt_5"] = (cooks > 5).sum(axis=1)

    return sample_info.sort_values(
        "cooks_q99",
        ascending=False,
    )


# ------------------------------------------------------------------
# PCA
# ------------------------------------------------------------------

def calculate_pca(normed, sample_names, sample_info, n_genes=1000):
    """
    Calculate PCA using log2-normalised counts.

    Returns
    -------
    pca_df
        PC scores and selected sample metadata.
    explained
        Fraction of variance explained by each PC.
    pca_genes
        Genes used for PCA.
    """

    log_normed = np.log2(normed + 1)

    gene_variance = (
        log_normed
        .var(axis=0)
        .sort_values(ascending=False)
    )

    n_pca_genes = min(n_genes, len(gene_variance))

    pca_genes = gene_variance.head(n_pca_genes).index

    X = log_normed[pca_genes].values

    # Centre each gene
    X = X - X.mean(axis=0)

    # SVD instead of sklearn
    U, S, _ = np.linalg.svd(
        X,
        full_matrices=False,
    )

    scores = U * S

    variance = S ** 2
    explained = variance / variance.sum()

    pca_df = pd.DataFrame(
        {
            "PC1": scores[:, 0],
            "PC2": scores[:, 1],
        },
        index=sample_names,
    )

    pca_metadata_cols = {
        CONTRAST_VARIABLE,
        "donor",
        "treatment",
        "cell_line",
        "n_cells",
        "library_size",
    }

    for col in sorted(
        pca_metadata_cols
    ):
        if col in sample_info.columns:
            pca_df[col] = sample_info[col]

    return pca_df, explained, pca_genes, log_normed


def save_pca_data(pca_df, cell_dir):
    """Save PCA coordinates."""

    pca_df.to_csv(
        cell_dir / "pca.csv"
    )


def plot_pca(
    pca_df,
    explained,
    cell_type,
    cell_dir,
    group_by,
    filename,
    legend_kwargs=None,
):
    """Create a PCA scatter plot grouped by a sample variable."""

    fig, ax = plt.subplots(figsize=(8, 6))

    for group in sorted(pca_df[group_by].unique()):

        idx = pca_df[group_by] == group

        scatter_kwargs = {
            "label": group,
            "s": 70,
        }

        if group_by == "donor":
            cmap = plt.get_cmap("tab20")
            group_values = sorted(
                pca_df[group_by].unique()
            )
            scatter_kwargs["color"] = cmap(
                group_values.index(group) % 20
            )

        ax.scatter(
            pca_df.loc[idx, "PC1"],
            pca_df.loc[idx, "PC2"],
            **scatter_kwargs,
        )

    annotate_samples(ax, pca_df)

    ax.set_xlabel(
        f"PC1 ({explained[0] * 100:.1f}%)"
    )
    ax.set_ylabel(
        f"PC2 ({explained[1] * 100:.1f}%)"
    )

    ax.set_title(
        f"{cell_type}\n"
        f"Pseudobulk PCA by {group_by}"
    )

    if legend_kwargs is None:
        legend_kwargs = {}

    ax.legend(
        title=group_by.capitalize(),
        **legend_kwargs,
    )

    fig.tight_layout()

    fig.savefig(
        cell_dir / filename,
        dpi=200,
        bbox_inches="tight",
    )

    plt.close(fig)


def annotate_samples(ax, pca_df):
    """Annotate PCA points with sample names."""

    for sample, row in pca_df.iterrows():
        ax.annotate(
            sample,
            (
                row["PC1"],
                row["PC2"],
            ),
            fontsize=7,
            xytext=(3, 3),
            textcoords="offset points",
        )


# ------------------------------------------------------------------
# Sample correlation
# ------------------------------------------------------------------

def plot_sample_correlation(
    log_normed,
    pca_genes,
    cell_type,
    cell_dir,
):
    """Calculate, save and plot sample-sample correlation."""

    corr = (
        log_normed[pca_genes]
        .T
        .corr()
    )

    corr.to_csv(
        cell_dir / "sample_correlation.csv"
    )

    fig, ax = plt.subplots(figsize=(9, 8))

    sns.heatmap(
        corr,
        cmap="vlag",
        center=0,
        vmin=-1,
        vmax=1,
        square=True,
        ax=ax,
    )

    ax.set_title(
        f"{cell_type}\n"
        "Pseudobulk sample correlation"
    )

    fig.tight_layout()

    fig.savefig(
        cell_dir / "sample_correlation.png",
        dpi=200,
        bbox_inches="tight",
    )

    plt.close(fig)

    return corr


# ------------------------------------------------------------------
# Library sizes / size factors
# ------------------------------------------------------------------

def plot_size_factors(
    sample_info,
    cell_type,
    cell_dir,
):
    """Plot raw library sizes and DESeq2 size factors."""

    plot_df = (
        sample_info
        .reset_index()
        .rename(columns={"index": "sample"})
    )

    fig, axes = plt.subplots(
        2,
        1,
        figsize=(10, 8),
        sharex=True,
    )

    axes[0].bar(
        plot_df["sample"],
        plot_df["library_size"],
    )

    axes[0].set_ylabel("Raw library size")
    axes[0].set_title(
        f"{cell_type}\n"
        "Pseudobulk library sizes"
    )

    axes[1].bar(
        plot_df["sample"],
        plot_df["size_factor"],
    )

    axes[1].axhline(
        1,
        color="black",
        linestyle="--",
        linewidth=1,
    )

    axes[1].set_ylabel("DESeq2 size factor")

    axes[1].tick_params(
        axis="x",
        rotation=90,
    )

    fig.tight_layout()

    fig.savefig(
        cell_dir / "size_factors.png",
        dpi=200,
        bbox_inches="tight",
    )

    plt.close(fig)


def calculate_library_metrics(sample_info):
    """
    Calculate descriptive library-size, size-factor,
    and nuclei-count relationships.
    """

    # --------------------------------------------------
    # Size factors
    # --------------------------------------------------

    size_factors = np.asarray(
        sample_info["size_factor"],
        dtype=float,
    )

    sf_min = np.nanmin(size_factors)
    sf_max = np.nanmax(size_factors)

    if sf_min > 0:
        sf_ratio = sf_max / sf_min
    else:
        sf_ratio = np.nan

    # --------------------------------------------------
    # Library sizes
    # --------------------------------------------------

    library_sizes = np.asarray(
        sample_info["library_size"],
        dtype=float,
    )

    lib_min = np.nanmin(library_sizes)
    lib_max = np.nanmax(library_sizes)

    if lib_min > 0:
        library_ratio = lib_max / lib_min
    else:
        library_ratio = np.nan

    # --------------------------------------------------
    # Library size vs size factor
    # --------------------------------------------------

    valid = (
        np.isfinite(library_sizes)
        &
        np.isfinite(size_factors)
        &
        (library_sizes > 0)
        &
        (size_factors > 0)
    )

    if valid.sum() >= 3:
        size_factor_library_corr = np.corrcoef(
            np.log10(
                library_sizes[valid]
            ),
            np.log10(
                size_factors[valid]
            ),
        )[0, 1]
    else:
        size_factor_library_corr = np.nan

    # --------------------------------------------------
    # Nuclei count vs library size
    # --------------------------------------------------

    n_cells = np.asarray(
        sample_info["n_cells"],
        dtype=float,
    )

    valid_cells = (
        np.isfinite(n_cells)
        &
        np.isfinite(library_sizes)
        &
        (n_cells > 0)
        &
        (library_sizes > 0)
    )

    if valid_cells.sum() >= 3:
        cell_count_library_corr = np.corrcoef(
            np.log10(
                n_cells[valid_cells]
            ),
            np.log10(
                library_sizes[valid_cells]
            ),
        )[0, 1]
    else:
        cell_count_library_corr = np.nan

    # --------------------------------------------------
    # Return all metrics
    # --------------------------------------------------

    return {
        "library_size_ratio": library_ratio,
        "size_factor_ratio": sf_ratio,
        "size_factor_library_corr": (
            size_factor_library_corr
        ),
        "cell_count_library_corr": (
            cell_count_library_corr
        ),
    }


# ------------------------------------------------------------------
# Dispersion
# ------------------------------------------------------------------

def plot_dispersion(
    dds,
    cell_type,
    cell_dir,
):
    """Plot gene-wise and final dispersion estimates."""

    var = dds.var

    x = np.asarray(
        var["_normed_means"],
        dtype=float,
    )

    y_genewise = np.asarray(
        var["genewise_dispersions"],
        dtype=float,
    )

    y_final = np.asarray(
        var["dispersions"],
        dtype=float,
    )

    valid_gene = (
        (x > 0)
        & np.isfinite(x)
        & np.isfinite(y_genewise)
        & (y_genewise > 0)
    )

    valid_final = (
        (x > 0)
        & np.isfinite(x)
        & np.isfinite(y_final)
        & (y_final > 0)
    )

    fig, ax = plt.subplots(figsize=(7, 6))

    ax.scatter(
        x[valid_gene],
        y_genewise[valid_gene],
        s=5,
        alpha=0.2,
        label="Gene-wise",
    )

    ax.scatter(
        x[valid_final],
        y_final[valid_final],
        s=5,
        alpha=0.3,
        label="Final",
    )

    ax.set_xscale("log")
    ax.set_yscale("log")

    ax.set_xlabel("Mean normalized count")
    ax.set_ylabel("Dispersion")

    ax.set_title(
        f"{cell_type}\n"
        "Dispersion estimates"
    )

    ax.legend()

    fig.tight_layout()

    fig.savefig(
        cell_dir / "dispersion.png",
        dpi=200,
        bbox_inches="tight",
    )

    plt.close(fig)


# ------------------------------------------------------------------
# Cook's distance
# ------------------------------------------------------------------

def plot_cooks(
    cooks,
    sample_names,
    cell_type,
    cell_dir,
):
    """Plot Cook's distance distribution for each sample."""

    cooks_plot = [
        (
            cooks.loc[sample]
            .replace(
                [np.inf, -np.inf],
                np.nan,
            )
            .dropna()
            .values
        )
        for sample in sample_names
    ]

    fig, ax = plt.subplots(figsize=(10, 6))

    ax.boxplot(
        cooks_plot,
        tick_labels=sample_names,
        showfliers=False,
    )

    ax.set_yscale("log")
    ax.set_ylabel("Cook's distance")

    ax.set_title(
        f"{cell_type}\n"
        "Cook's distance by sample"
    )

    ax.tick_params(
        axis="x",
        rotation=90,
    )

    fig.tight_layout()

    fig.savefig(
        cell_dir / "cooks_by_sample.png",
        dpi=200,
        bbox_inches="tight",
    )

    plt.close(fig)


def calculate_cooks_metrics(sample_info):
    """Calculate descriptive sample-level Cook's metrics."""

    return {
        "worst_cooks_sample": (
            sample_info["n_cooks_gt_1"].idxmax()
        ),
        "worst_sample_n_cooks_gt1": int(
            sample_info["n_cooks_gt_1"].max()
        ),
        "worst_cooks_q99_sample": (
            sample_info["cooks_q99"].idxmax()
        ),
        "max_cooks_q99": float(
            sample_info["cooks_q99"].max()
        ),
        "max_cooks_sample": (
            sample_info["max_cooks"].idxmax()
        ),
        "max_cooks_value": float(
            sample_info["max_cooks"].max()
        ),
    }


# ------------------------------------------------------------------
# Gene-level model diagnostics
# ------------------------------------------------------------------

def count_true(dds, column):
    """Count True values, treating missing values as False."""

    if column not in dds.var.columns:
        return np.nan

    values = pd.Series(
        np.asarray(dds.var[column])
    ).fillna(False)

    return int(
        values.astype(bool).sum()
    )


def count_false(dds, column):
    """Count explicit False values, ignoring missing values."""

    if column not in dds.var.columns:
        return np.nan

    values = pd.Series(
        np.asarray(dds.var[column])
    )

    values = values[values.notna()]

    if len(values) == 0:
        return np.nan

    return int(
        (~values.astype(bool)).sum()
    )


def calculate_model_diagnostics(dds):
    """Collect PyDESeq2 model-fitting diagnostics."""

    n_genes = dds.n_vars

    diagnostics = {
        "n_genes": n_genes,

        "n_dispersion_outliers": count_true(
            dds,
            "_outlier_genes",
        ),

        "n_cooks_outlier_genes": count_true(
            dds,
            "_pvalue_cooks_outlier",
        ),

        "n_lfc_not_converged": count_false(
            dds,
            "_LFC_converged",
        ),

        "n_MAP_not_converged": count_false(
            dds,
            "_MAP_converged",
        ),

        "n_genewise_not_converged": count_false(
            dds,
            "_genewise_converged",
        ),
    }

    return diagnostics


# ------------------------------------------------------------------
# Percentages
# ------------------------------------------------------------------

def percentage(value, denominator):
    """Calculate a percentage safely."""

    if (
        pd.isna(value)
        or denominator == 0
    ):
        return np.nan

    return 100 * value / denominator


def add_model_percentages(diagnostics):
    """Add percentages to gene-level diagnostic counts."""

    n_genes = diagnostics["n_genes"]

    diagnostics["pct_dispersion_outliers"] = percentage(
        diagnostics["n_dispersion_outliers"],
        n_genes,
    )

    diagnostics["pct_cooks_outliers"] = percentage(
        diagnostics["n_cooks_outlier_genes"],
        n_genes,
    )

    diagnostics["pct_lfc_not_converged"] = percentage(
        diagnostics["n_lfc_not_converged"],
        n_genes,
    )

    diagnostics["pct_map_not_converged"] = percentage(
        diagnostics["n_MAP_not_converged"],
        n_genes,
    )

    diagnostics["pct_genewise_not_converged"] = percentage(
        diagnostics["n_genewise_not_converged"],
        n_genes,
    )

    return diagnostics


# ------------------------------------------------------------------
# Automated triage
# ------------------------------------------------------------------

def assess_qc(diagnostics):
    """
    Determine automated review status from model-fitting diagnostics.

    Descriptive library-size, size-factor and Cook's thresholds
    intentionally do not affect this classification.
    """

    flags = []

    if (
        np.isfinite(diagnostics["pct_dispersion_outliers"])
        and diagnostics["pct_dispersion_outliers"] > 5
    ):
        flags.append("DISPERSION_OUTLIERS")

    if (
        np.isfinite(diagnostics["pct_cooks_outliers"])
        and diagnostics["pct_cooks_outliers"] > 1
    ):
        flags.append("COOKS_OUTLIERS")

    if (
        np.isfinite(diagnostics["pct_lfc_not_converged"])
        and diagnostics["pct_lfc_not_converged"] > 1
    ):
        flags.append("LFC_CONVERGENCE")

    if (
        np.isfinite(diagnostics["pct_map_not_converged"])
        and diagnostics["pct_map_not_converged"] > 1
    ):
        flags.append("MAP_CONVERGENCE")

    if (
        np.isfinite(diagnostics["pct_genewise_not_converged"])
        and diagnostics["pct_genewise_not_converged"] > 5
    ):
        flags.append("GENEWISE_CONVERGENCE")

    if not flags:
        review_status = "OK"
    elif len(flags) == 1:
        review_status = "CHECK"
    else:
        review_status = "REVIEW"

    return review_status, flags


# ------------------------------------------------------------------
# Diagnostics report
# ------------------------------------------------------------------

def write_diagnostics_report(
    cell_dir,
    cell_type,
    dds,
    library_metrics,
    model_diagnostics,
    cooks_metrics,
    review_status,
    flags,
):
    """Write human-readable diagnostics.txt."""

    diagnostics = [
        f"Cell type: {cell_type}",
        f"Samples: {dds.n_obs}",
        f"Genes fitted: {model_diagnostics['n_genes']}",
        "",

        "DESCRIPTIVE SAMPLE METRICS",

        (
            "Library size max/min ratio: "
            f"{library_metrics['library_size_ratio']:.2f}"
        ),

        (
            "Size factor max/min ratio: "
            f"{library_metrics['size_factor_ratio']:.2f}"
        ),

        (
            "log10(library size) vs "
            "log10(size factor) correlation: "
            f"{library_metrics['size_factor_library_corr']:.3f}"
        ),

        (
            "log10(nuclei count) vs "
            "log10(library size) correlation: "
            f"{library_metrics['cell_count_library_corr']:.3f}"
        ),

        "",

        "MODEL DIAGNOSTICS",

        (
            "Dispersion outlier genes: "
            f"{model_diagnostics['n_dispersion_outliers']} "
            f"({model_diagnostics['pct_dispersion_outliers']:.3f}%)"
        ),

        (
            "Cook's-filtered genes: "
            f"{model_diagnostics['n_cooks_outlier_genes']} "
            f"({model_diagnostics['pct_cooks_outliers']:.3f}%)"
        ),

        (
            "LFC non-converged genes: "
            f"{model_diagnostics['n_lfc_not_converged']} "
            f"({model_diagnostics['pct_lfc_not_converged']:.3f}%)"
        ),

        (
            "MAP dispersion non-converged genes: "
            f"{model_diagnostics['n_MAP_not_converged']} "
            f"({model_diagnostics['pct_map_not_converged']:.3f}%)"
        ),

        (
            "Gene-wise dispersion non-converged genes: "
            f"{model_diagnostics['n_genewise_not_converged']} "
            f"({model_diagnostics['pct_genewise_not_converged']:.3f}%)"
        ),

        "",

        "COOK'S DISTANCE CONTEXT",

        (
            "Sample with most Cook's >1 values: "
            f"{cooks_metrics['worst_cooks_sample']}"
        ),

        (
            "Genes with Cook's >1 in that sample: "
            f"{cooks_metrics['worst_sample_n_cooks_gt1']}"
        ),

        (
            "Sample containing maximum Cook's value: "
            f"{cooks_metrics['max_cooks_sample']}"
        ),

        (
            "Maximum Cook's distance: "
            f"{cooks_metrics['max_cooks_value']:.3f}"
        ),

        "",

        "AUTOMATED TRIAGE",

        (
            "Review status: "
            f"{review_status}"
        ),

        (
            "Flags: "
            f"{', '.join(flags) if flags else 'None'}"
        ),

        "",

        (
            "NOTE: library-size ratio, size-factor ratio, "
            "and arbitrary Cook's >1 counts are descriptive "
            "and do not trigger QC failure."
        ),
    ]

    with open(
        cell_dir / "diagnostics.txt",
        "w",
    ) as handle:
        handle.write(
            "\n".join(diagnostics)
        )


# ------------------------------------------------------------------
# Global summary
# ------------------------------------------------------------------

def build_inspection_summary(
    cell_type,
    dds,
    library_metrics,
    pca_explained,
    model_diagnostics,
    cooks_metrics,
    review_status,
    flags,
):
    """Build the row returned to the global inspection summary."""

    return {
        "cell_type": cell_type,

        # Automated triage
        "qc_status": review_status,
        "qc_flags": "|".join(flags) if flags else "",

        # General
        "n_samples_qc": dds.n_obs,
        "n_genes_qc": model_diagnostics["n_genes"],

        # Descriptive library metrics
        "library_size_ratio": (
            library_metrics["library_size_ratio"]
        ),
        "size_factor_ratio": (
            library_metrics["size_factor_ratio"]
        ),
        "size_factor_library_corr": (
            library_metrics["size_factor_library_corr"]
        ),
        "cell_count_library_corr": (
            library_metrics["cell_count_library_corr"]
        ),

        # PCA
        "PC1_variance_pct": (
            pca_explained[0] * 100
        ),
        "PC2_variance_pct": (
            pca_explained[1] * 100
        ),

        # Dispersion
        "n_dispersion_outliers": (
            model_diagnostics["n_dispersion_outliers"]
        ),
        "pct_dispersion_outliers": (
            model_diagnostics["pct_dispersion_outliers"]
        ),

        # Cook's model-level result
        "n_cooks_outlier_genes": (
            model_diagnostics["n_cooks_outlier_genes"]
        ),
        "pct_cooks_outlier_genes": (
            model_diagnostics["pct_cooks_outliers"]
        ),

        # Convergence
        "n_lfc_not_converged": (
            model_diagnostics["n_lfc_not_converged"]
        ),
        "pct_lfc_not_converged": (
            model_diagnostics["pct_lfc_not_converged"]
        ),
        "n_MAP_not_converged": (
            model_diagnostics["n_MAP_not_converged"]
        ),
        "pct_MAP_not_converged": (
            model_diagnostics["pct_map_not_converged"]
        ),
        "n_genewise_not_converged": (
            model_diagnostics["n_genewise_not_converged"]
        ),
        "pct_genewise_not_converged": (
            model_diagnostics["pct_genewise_not_converged"]
        ),

        # Cook's descriptive context
        "worst_cooks_sample": (
            cooks_metrics["worst_cooks_sample"]
        ),
        "worst_sample_n_cooks_gt1": (
            cooks_metrics["worst_sample_n_cooks_gt1"]
        ),
        "worst_cooks_q99_sample": (
            cooks_metrics["worst_cooks_q99_sample"]
        ),
        "max_cooks_q99": (
            cooks_metrics["max_cooks_q99"]
        ),
        "max_cooks_sample": (
            cooks_metrics["max_cooks_sample"]
        ),
        "max_cooks_value": (
            cooks_metrics["max_cooks_value"]
        ),
    }




# ------------------------------------------------------------------
# Main function
# ------------------------------------------------------------------

def generate_de_inspection(
    dds,
    meta_pb,
    cell_type,
    output_dir=INSPECTION_DIR,
):
    """
    Generate pseudobulk / PyDESeq2 inspection outputs for one cell type.
    """

    cell_dir = (
        Path(output_dir)
        / safe_filename(cell_type)
    )

    cell_dir.mkdir(
        parents=True,
        exist_ok=True,
    )

    # --------------------------------------------------------------
    # Data
    # --------------------------------------------------------------

    sample_names, normed, cooks, sample_info = (
        prepare_inspection_data(
            dds,
            meta_pb,
        )
    )

    # --------------------------------------------------------------
    # Sample QC
    # --------------------------------------------------------------

    sample_info = calculate_sample_qc(
        sample_info,
        normed,
        cooks,
    )

    sample_info.to_csv(
        cell_dir / "sample_qc.csv"
    )

    # --------------------------------------------------------------
    # PCA
    # --------------------------------------------------------------

    pca_df, explained, pca_genes, log_normed = calculate_pca(
        normed,
        sample_names,
        sample_info,
    )

    save_pca_data(
        pca_df,
        cell_dir,
    )

    plot_pca(
        pca_df,
        explained,
        cell_type,
        cell_dir,
        group_by=CONTRAST_VARIABLE,
        filename=(f"pca_{safe_filename(CONTRAST_VARIABLE)}.png")
    )

    if "donor" in pca_df.columns:

        plot_pca(
            pca_df,
            explained,
            cell_type,
            cell_dir,
            group_by="donor",
            filename="pca_donor.png",
            legend_kwargs={
                "bbox_to_anchor": (1.02, 1),
                "loc": "upper left",
                "fontsize": 7,
            },
        )

    # --------------------------------------------------------------
    # Correlation
    # --------------------------------------------------------------

    plot_sample_correlation(
        log_normed,
        pca_genes,
        cell_type,
        cell_dir,
    )

    # --------------------------------------------------------------
    # Library / size factors
    # --------------------------------------------------------------

    plot_size_factors(
        sample_info,
        cell_type,
        cell_dir,
    )

    library_metrics = calculate_library_metrics(
        sample_info,
    )

    # --------------------------------------------------------------
    # Model diagnostics
    # --------------------------------------------------------------

    plot_dispersion(
        dds,
        cell_type,
        cell_dir,
    )

    plot_cooks(
        cooks,
        sample_names,
        cell_type,
        cell_dir,
    )

    model_diagnostics = calculate_model_diagnostics(
        dds,
    )

    model_diagnostics = add_model_percentages(
        model_diagnostics,
    )

    cooks_metrics = calculate_cooks_metrics(
        sample_info,
    )

    # --------------------------------------------------------------
    # Automated triage
    # --------------------------------------------------------------

    review_status, flags = assess_qc(
        model_diagnostics,
    )

    # --------------------------------------------------------------
    # Human-readable report
    # --------------------------------------------------------------

    write_diagnostics_report(
        cell_dir=cell_dir,
        cell_type=cell_type,
        dds=dds,
        library_metrics=library_metrics,
        model_diagnostics=model_diagnostics,
        cooks_metrics=cooks_metrics,
        review_status=review_status,
        flags=flags,
    )

    # --------------------------------------------------------------
    # Global summary row
    # --------------------------------------------------------------

    return build_inspection_summary(
        cell_type=cell_type,
        dds=dds,
        library_metrics=library_metrics,
        pca_explained=explained,
        model_diagnostics=model_diagnostics,
        cooks_metrics=cooks_metrics,
        review_status=review_status,
        flags=flags,
    )





# Import adata
adata = sc.read_h5ad(args.input_file)

# Check samples requested are present
if "sample" not in adata.obs.columns:
    raise ValueError(
        "Required column 'sample' "
        "is not present in adata.obs."
    )

# Restore raw counts
adata_counts = restore_counts(
    adata
)

# Defensive check
if not adata_counts.obs_names.equals(
    adata.obs_names
):
    raise ValueError(
        "restore_counts() changed observation order "
        "or observation identities."
    )

# Ensure corrected metadata is retained
adata_counts.obs = adata.obs.copy()

# --------------------------------------------------
# Select requested samples
# --------------------------------------------------

if args.samples is not None:

    requested_samples = set(
        args.samples
    )

    available_samples = set(
        adata_counts.obs[
            "sample"
        ]
        .astype(str)
        .unique()
    )

    missing_samples = (
        requested_samples
        -
        available_samples
    )

    if missing_samples:
        raise ValueError(
            "Requested samples not present "
            "in AnnData: "
            f"{sorted(missing_samples)}"
        )

    adata_de = adata_counts[
        adata_counts.obs[
            "sample"
        ]
        .astype(str)
        .isin(
            requested_samples
        )
    ].copy()

else:

    adata_de = (
        adata_counts.copy()
    )

# --------------------------------------------------
# Keep only requested contrast levels
# --------------------------------------------------
if (
    CONTRAST_VARIABLE
    not in adata_de.obs.columns
):
    raise ValueError(
        f"Contrast variable "
        f"'{CONTRAST_VARIABLE}' "
        "not present in adata.obs."
    )

adata_de = adata_de[
    adata_de.obs[
        CONTRAST_VARIABLE
    ]
    .astype(str)
    .isin(
        [
            CONTRAST_TEST,
            CONTRAST_REFERENCE,
        ]
    )
].copy()

# Store counts explicitly
adata_de.layers["counts"] = (
    adata_de.X.copy()
)

if adata_de.n_obs == 0:
    raise ValueError(
        "No observations remain after sample "
        "and contrast filtering."
    )

if (
    args.annotation
    not in adata_de.obs.columns
):
    raise ValueError(
        f"Annotation column "
        f"'{args.annotation}' "
        "not found in adata.obs."
    )

# Get sample cell counts
sample_counts = (
    adata_de.obs
    .groupby([args.annotation, "sample"], observed=True)
    .size()
    .unstack(fill_value=0)
)
file = (OUTDIR / f"{Path(args.input_file).stem}_{args.annotation}_counts.csv")
sample_counts.to_csv(file)

# --------------------------------------------------
# use raw counts
# --------------------------------------------------

counts = adata_de.layers["counts"]

if sparse.issparse(counts):
    count_values = counts.data
else:
    count_values = np.asarray(counts).ravel()

if np.any(count_values < 0):
    raise ValueError(
        "Negative values detected in counts layer."
    )

if not np.allclose(
    count_values,
    np.round(count_values)
):
    raise ValueError(
        "The counts layer does not contain integer-like raw counts. "
        "Do not use normalized/log-transformed expression for pseudobulk DE."
    )

print(
    "Counts validated:",
    f"min={count_values.min()},",
    f"max={count_values.max()}"
)

# --------------------------------------------------
# Verify metadata present for all samples
# --------------------------------------------------

sample_design_cols = [
    "sample"
]

for col in sorted(
    REQUIRED_MODEL_COLUMNS
    |
    {
        "donor",
        "cell_line",
    }
):

    if (
        col in adata_de.obs.columns
        and
        col not in sample_design_cols
    ):
        sample_design_cols.append(
            col
        )

sample_design = (
    adata_de.obs[
        sample_design_cols
    ]
    .drop_duplicates()
    .sort_values(
        "sample"
    )
)

print(
    sample_design.to_string(
        index=False
    )
)

sample_design.to_csv(
    OUTDIR /
    "sample_design.csv",
    index=False,
)


# --------------------------------------------------
# Validate sample-level metadata
# --------------------------------------------------

metadata_cols = sorted(
    REQUIRED_MODEL_COLUMNS
    |
    {
        "sample",
        "donor",
        "cell_line",
    }
)

for col in metadata_cols:

    if col not in adata_de.obs.columns:
        continue

    n_values = (
        adata_de.obs
        .groupby(
            "sample",
            observed=True,
        )[col]
        .nunique(
            dropna=False
        )
    )

    bad = n_values[
        n_values != 1
    ]

    if len(bad) > 0:
        raise ValueError(
            f"Inconsistent {col} metadata for samples: "
            f"{bad.index.tolist()}"
        )

# --------------------------------------------------
# Check required model metadata is not missing
# --------------------------------------------------

for col in sorted(
    REQUIRED_MODEL_COLUMNS
):

    if col not in adata_de.obs.columns:
        raise ValueError(
            f"Required model metadata column "
            f"'{col}' is absent from adata.obs."
        )

    missing_samples = (
        adata_de.obs
        .groupby(
            "sample",
            observed=True,
        )[col]
        .apply(
            lambda x: x.isna().all()
        )
    )

    missing_samples = (
        missing_samples[
            missing_samples
        ]
        .index
        .tolist()
    )

    if missing_samples:
        raise ValueError(
            f"Missing required model metadata "
            f"'{col}' for samples: "
            f"{missing_samples}"
        )

# --------------------------------------------------
# Validate contrast levels
# --------------------------------------------------

observed_contrast_levels = set(
    adata_de.obs[
        CONTRAST_VARIABLE
    ]
    .astype(str)
)

expected_contrast_levels = {
    CONTRAST_TEST,
    CONTRAST_REFERENCE,
}

unexpected_levels = (
    observed_contrast_levels
    -
    expected_contrast_levels
)

if unexpected_levels:
    raise ValueError(
        f"Unexpected {CONTRAST_VARIABLE} levels "
        "after contrast subsetting: "
        f"{sorted(unexpected_levels)}"
    )

missing_levels = (
    expected_contrast_levels
    -
    observed_contrast_levels
)

if missing_levels:
    raise ValueError(
        f"Missing required {CONTRAST_VARIABLE} levels: "
        f"{sorted(missing_levels)}"
    )

# Update manifest
selected_samples = sorted(
    adata_de.obs["sample"]
    .astype(str)
    .unique()
    .tolist()
)

run_manifest["selected_samples"] = (selected_samples)
run_manifest["n_selected_samples"] = (len(selected_samples))

with open(
    OUTDIR / "run_manifest.json",
    "w",
) as handle:

    json.dump(
        run_manifest,
        handle,
        indent=2,
    )



# --------------------------------------------------
# Differential expression + GSEA
# --------------------------------------------------

all_summary = []
inspection_summary = []

cell_types = sorted(
    adata_de.obs[args.annotation]
    .dropna()
    .unique()
)

for cell_type in cell_types:
    print()
    print("=" * 80)
    print(cell_type)
    print("=" * 80)
    cell_safe = safe_filename(cell_type)
    # --------------------------------------------------
    # Pseudobulk
    # --------------------------------------------------
    out = pseudobulk_celltype(
        adata_de,
        cell_type,
    )
    if out is None:

        print(
            f"{cell_type}: "
            "no eligible pseudobulk samples"
        )

        all_summary.append(
            {
                **RUN_INFO,
                "cell_type": cell_type,
                "analysis_status": "NO_PSEUDOBULK",
                "inspection_status": "NOT_RUN",
                "n_samples": 0,
                "n_genes_tested": 0,
                "n_FDR_005": 0,
                f"n_FDR_005_{safe_filename(CONTRAST_TEST)}_up": 0,
                f"n_FDR_005_{safe_filename(CONTRAST_TEST)}_down": 0,
                "n_FDR_005_log2FC1": 0,
                "n_nominal": 0,
            }
        )

        continue

    pb, meta_pb = out
    # --------------------------------------------------
    # Save sample-level information
    # --------------------------------------------------
    meta_pb.to_csv(PB_DIR / f"{cell_safe}_samples.csv")
    # --------------------------------------------------
    # Assess experimental design
    # --------------------------------------------------
    design_info = assess_design(meta_pb)
    print(
        f"{cell_type}: "
        f"{design_info['n_samples']} samples; "
        f"{design_info[f'n_{CONTRAST_REFERENCE}']} "
        f"{CONTRAST_REFERENCE}; "
        f"{design_info[f'n_{CONTRAST_TEST}']} "
        f"{CONTRAST_TEST}; "
        f"{design_info.get('n_informative_donors', 'NA')} "
        "informative donors; "
        f"{design_info.get('n_1to1_paired_donors', 'NA')} "
        "1:1 pairs"
    )
    if not design_info["eligible"]:
        print(
            f"{cell_type}: skipping DE "
            f"due to insufficient replication "
            f"for {CONTRAST_TEST} vs "
            f"{CONTRAST_REFERENCE}"
        )
        all_summary.append(
            {
                **RUN_INFO,
                "cell_type": cell_type,
                **design_info,
                "analysis_status": "INELIGIBLE",
                "inspection_status": "NOT_RUN",
                "n_genes_tested": 0,
                "n_FDR_005": 0,
                f"n_FDR_005_{safe_filename(CONTRAST_TEST)}_up": 0,
                f"n_FDR_005_{safe_filename(CONTRAST_TEST)}_down": 0,
                "n_FDR_005_log2FC1": 0,
                "n_nominal": 0,
            }
        )
        continue
    # --------------------------------------------------
    # Gene filtering
    # --------------------------------------------------
    pb_filtered = filter_genes(pb)

    if pb_filtered.shape[1] == 0:
        print(
            f"{cell_type}: skipping DE; "
            "no genes passed expression filtering"
        )

        all_summary.append(
            {
                **RUN_INFO,
                "cell_type": cell_type,
                **design_info,
                "analysis_status": "NO_GENES",
                "inspection_status": "NOT_RUN",
                "n_genes_tested": 0,
                "n_FDR_005": 0,
                f"n_FDR_005_{safe_filename(CONTRAST_TEST)}_up": 0,
                f"n_FDR_005_{safe_filename(CONTRAST_TEST)}_down": 0,
                "n_FDR_005_log2FC1": 0,
                "n_nominal": 0,
            }
        )

        continue

    print(
        f"{cell_type}: "
        f"{pb_filtered.shape[1]} genes retained"
    )
    # Save pseudobulk counts if useful for auditing
    pb_filtered.to_csv(
        PB_DIR /
        f"{cell_safe}_counts.csv"
    )
    # --------------------------------------------------
    # Differential expression
    # --------------------------------------------------
    inspection_status = "NOT_RUN"
    try:
        de, dds, ds = run_pydeseq2(
            pb_filtered,
            meta_pb,
        )
        # --------------------------------------------------
        # Automated DE inspection/QC
        # --------------------------------------------------
        inspection_status = "SUCCESS"

        try:

            qc = generate_de_inspection(
                dds=dds,
                meta_pb=meta_pb,
                cell_type=cell_type,
            )

            inspection_summary.append(
                qc
            )

        except Exception:

            inspection_status = "FAILED"

            print(
                f"{cell_type}: "
                "inspection output failed"
            )

            traceback.print_exc()
    except Exception:

        print(
            f"{cell_type}: PyDESeq2 failed"
        )

        traceback.print_exc()

        all_summary.append(
            {
                **RUN_INFO,
                "cell_type": cell_type,
                **design_info,
                "inspection_status": (inspection_status),
                "analysis_status": "DE_FAILED",
                "n_genes_tested": 0,
                "n_FDR_005": np.nan,
                f"n_FDR_005_{safe_filename(CONTRAST_TEST)}_up": np.nan,
                f"n_FDR_005_{safe_filename(CONTRAST_TEST)}_down": np.nan,
                "n_FDR_005_log2FC1": np.nan,
                "n_nominal": np.nan,
            }
        )

        continue
    # --------------------------------------------------
    # Sort results
    # --------------------------------------------------
    de = de.sort_values(["FDR", "pval"], na_position="last")
    # --------------------------------------------------
    # Complete DE table
    # --------------------------------------------------
    de.to_csv(DE_DIR / f"{cell_safe}_all_genes.csv", index=False)
    # --------------------------------------------------
    # FDR-significant DE
    # --------------------------------------------------
    sig = de[de["FDR"] < 0.05].copy()
    sig.to_csv(DE_DIR / f"{cell_safe}_significant.csv", index=False)
    # --------------------------------------------------
    # Direction of FDR-significant DE
    # --------------------------------------------------
    n_sig_up = ( ( (de["FDR"] < 0.05) & (de["log2FoldChange"] > 0) ).sum())
    n_sig_down = ( ( (de["FDR"] < 0.05) & (de["log2FoldChange"] < 0) ).sum())

    # --------------------------------------------------
    # FDR-significant + large effect
    # --------------------------------------------------
    sig_large = de[ (de["FDR"] < 0.05) & (np.abs(de["log2FoldChange"]) > 1) ].copy()
    sig_large.to_csv(
        DE_DIR /
        (
            f"{cell_safe}_"
            "significant_log2FC1.csv"
        ),
        index=False,
    )
    # --------------------------------------------------
    # Nominal/exploratory DE
    # --------------------------------------------------
    nominal = de[ (de["pval"] < 0.05) & ( np.abs(de["log2FoldChange"]) > 0.5) ].copy()
    nominal.to_csv(
        DE_DIR /
        f"{cell_safe}_nominal.csv",
        index=False,
    )
    # --------------------------------------------------
    # GSEA ranking
    # --------------------------------------------------
    ranking, duplicate_rank_pct = (
        make_gsea_ranking(de)
    )

    print(
        f"{cell_type}: "
        f"{duplicate_rank_pct:.2f}% "
        "duplicated GSEA ranking values"
    )

    ranking.to_csv(
        GSEA_DIR /
        f"{cell_safe}_ranking.csv",
        index=False,
    )

    # --------------------------------------------------
    # GSEA
    # --------------------------------------------------

    gsea_summary = {}

    if len(ranking) < 10:

        print(
            f"{cell_type}: skipping GSEA; "
            f"only {len(ranking)} genes have "
            "usable Wald statistics"
        )

        for gs_name in gene_sets:

            gsea_summary[f"{gs_name}_sig_pathways"] = np.nan
            gsea_summary[f"{gs_name}_sig_up"] = np.nan
            gsea_summary[f"{gs_name}_sig_down"] = np.nan

    else:

        for gs_name, gs_db in gene_sets.items():

            print(
                f"{cell_type}: "
                f"running {gs_name}"
            )

            gs_dir = (GSEA_DIR / cell_safe / gs_name)
            gs_dir.mkdir(parents=True, exist_ok=True,)

            try:

                pre_res = gp.prerank(
                    rnk=ranking,
                    gene_sets=gs_db,
                    min_size=10,
                    max_size=500,
                    permutation_num=args.gsea_permutations,
                    seed=0,
                    threads=args.n_cpus,
                    outdir=str(gs_dir),
                    verbose=False,
                )

                gsea_res = (pre_res.res2d.copy())

                # --------------------------------------------------
                # Save all pathways
                # --------------------------------------------------
                gsea_res.to_csv(
                    gs_dir /
                    "all_pathways.csv",
                    index=False,
                )

                # --------------------------------------------------
                # FDR-significant pathways
                # --------------------------------------------------
                sig_pathways = gsea_res[gsea_res["FDR q-val"] < 0.05].copy()
                sig_pathways.to_csv(
                    gs_dir /
                    "significant_pathways.csv",
                    index=False,
                )

                # --------------------------------------------------
                # Significant pathways increased in TEST level
                #
                # Positive NES:
                # genes towards the positive end of the ranking,
                # i.e. TEST > REFERENCE
                # --------------------------------------------------
                sig_up = (
                    sig_pathways[sig_pathways["NES"] > 0]
                    .sort_values(["FDR q-val", "NES"], ascending=[True, False],)
                )
                sig_up.to_csv(
                    gs_dir /
                    f"significant_{safe_filename(CONTRAST_TEST)}_up.csv",
                    index=False,
                )

                # --------------------------------------------------
                # Significant pathways decreased in TEST level
                #
                # Negative NES:
                # genes towards the negative end of the ranking,
                # i.e. TEST < REFERENCE / REFERENCE enriched
                # --------------------------------------------------
                sig_down = (
                    sig_pathways[sig_pathways["NES"] < 0]
                    .sort_values(["FDR q-val", "NES"], ascending=[True, True],)
                )
                sig_down.to_csv(
                    gs_dir /
                    f"significant_{safe_filename(CONTRAST_TEST)}_down.csv",
                    index=False,
                )

                # --------------------------------------------------
                # Top 25 positively enriched pathways
                #
                # These are drawn from ALL pathways, not only those
                # reaching FDR < 0.05, so that an inspection table
                # is always produced.
                #
                # Primary ordering: FDR
                # Secondary ordering: strongest positive NES
                # --------------------------------------------------
                top_up = (
                    gsea_res[gsea_res["NES"] > 0]
                    .sort_values(["FDR q-val", "NES"], ascending=[True, False],)
                    .head(25)
                    .copy()
                )
                top_up.to_csv(
                    gs_dir /
                    f"top25_{safe_filename(CONTRAST_TEST)}_up.csv",
                    index=False,
                )

                # --------------------------------------------------
                # Top 25 negatively enriched pathways
                #
                # Primary ordering: FDR
                # Secondary ordering: strongest negative NES
                # --------------------------------------------------
                top_down = (
                    gsea_res[gsea_res["NES"] < 0]
                    .sort_values(["FDR q-val", "NES"],ascending=[True, True],)
                    .head(25)
                    .copy()
                )
                top_down.to_csv(
                    gs_dir /
                    f"top25_{safe_filename(CONTRAST_TEST)}_down.csv",
                    index=False,
                )

                # --------------------------------------------------
                # Summary statistics
                # --------------------------------------------------
                gsea_summary[f"{gs_name}_sig_pathways"] = len(sig_pathways)
                gsea_summary[f"{gs_name}_sig_up"] = len(sig_up)
                gsea_summary[f"{gs_name}_sig_down"] = len(sig_down)
                print(
                    f"{cell_type}: "
                    f"{len(sig_pathways)} significant "
                    f"{gs_name} pathways "
                    f"({len(sig_up)} "
                    f"{CONTRAST_TEST}-up, "
                    f"{len(sig_down)} "
                    f"{CONTRAST_TEST}-down)"
                )

            except Exception as e:

                print(
                    f"{cell_type}: "
                    f"{gs_name} failed"
                )

                traceback.print_exc()

                gsea_summary[
                    f"{gs_name}_sig_pathways"
                ] = np.nan

                gsea_summary[
                    f"{gs_name}_sig_up"
                ] = np.nan

                gsea_summary[
                    f"{gs_name}_sig_down"
                ] = np.nan

    # --------------------------------------------------
    # Summary
    # --------------------------------------------------
    all_summary.append(
        {
            **RUN_INFO,
            "cell_type": cell_type,
            **design_info,
            "analysis_status": "SUCCESS",
            "inspection_status": inspection_status,
            "n_genes_tested": len(de),
            "n_FDR_005": len(sig),
            f"n_FDR_005_{safe_filename(CONTRAST_TEST)}_up": (int(n_sig_up)),
            f"n_FDR_005_{safe_filename(CONTRAST_TEST)}_down": (int(n_sig_down)),
            "n_FDR_005_log2FC1": len(sig_large),
            "n_nominal": len(nominal),
            "gsea_duplicate_rank_pct": duplicate_rank_pct,
            **gsea_summary,
        }
    )

summary = pd.DataFrame(
    all_summary
)

inspection_df = pd.DataFrame(
    inspection_summary
)

if not inspection_df.empty:

    summary = summary.merge(
        inspection_df,
        on="cell_type",
        how="left",
    )

summary.to_csv(
    OUTDIR /
    f"{args.annotation}_DE_GSEA_summary.csv",
    index=False,
)
