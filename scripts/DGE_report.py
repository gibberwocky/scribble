#!/usr/bin/env python3
"""
Generate filtered interpretation tables, plots, and an HTML report from
compare_DGE.py outputs.

This is a reporting/triage layer only. It does not alter DE or GSEA results.
Filtering thresholds determine which results are surfaced for interpretation,
not statistical significance in the underlying analyses.
"""

import argparse
import html
import json
import re
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns


# -----------------------------------------------------------------------------
# Helpers
# -----------------------------------------------------------------------------

def safe_filename(value):
    value = str(value).strip()
    value = re.sub(r"[^\w.-]+", "_", value)
    value = re.sub(r"_+", "_", value)
    return value.strip("_.")


def numeric(series):
    return pd.to_numeric(series, errors="coerce")


def parse_args():
    parser = argparse.ArgumentParser()

    parser.add_argument("--comparison_dir", required=True, help="Output directory created by compare_DGE.py.")
    parser.add_argument("--reference_analysis", required=True, help="Prespecified primary/reference analysis.")
    parser.add_argument("--analyses", nargs="+", required=True, help="Analyses to include, in preferred display order.")
    parser.add_argument("--report_name", required=True, help="Name for this interpretation report.")

    # Statistical/effect thresholds used only for reporting/triage
    parser.add_argument("--fdr", type=float, default=0.05)
    parser.add_argument("--nominal_p", type=float, default=0.05)
    parser.add_argument("--min_abs_log2fc", type=float, default=0.5)
    parser.add_argument("--min_abs_wald", type=float, default=2.0)
    parser.add_argument("--min_delta_log2fc", type=float, default=0.5)
    parser.add_argument("--min_abs_nes", type=float, default=1.5)
    parser.add_argument("--min_delta_nes", type=float, default=0.5)
    parser.add_argument("--min_supporting_analyses", type=int, default=2)

    # Display limits
    parser.add_argument("--top_genes", type=int, default=20)
    parser.add_argument("--top_pathways", type=int, default=15)
    parser.add_argument("--heatmap_genes", type=int, default=25)
    parser.add_argument("--heatmap_pathways", type=int, default=25)

    # Plot control
    parser.add_argument("--make_plots", action=argparse.BooleanOptionalAction, default=True)

    # Leading-edge analysis parameters
    parser.add_argument("--leading_edge_min_analyses", type=int, default=2,
        help=(
            "Minimum number of analyses in which a gene must occur "
            "in the leading edge of concordant significant pathways "
            "to be included in recurrent leading-edge outputs."),
    )

    parser.add_argument("--leading_edge_min_pathways", type=int, default=2,
        help=(
            "Minimum number of concordant significant pathways in "
            "which a gene must occur to be included in recurrent "
            "leading-edge outputs."
        ),
    )

    parser.add_argument("--top_leading_edge_genes", type=int, default=25,
        help=(
            "Maximum number of recurrent leading-edge genes shown "
            "per cell type in the HTML report."
        ),
    )

    parser.add_argument("--leading_edge_nonribosomal_genes", type=int, default=15,
        help=(
            "Maximum number of non-ribosomal recurrent "
            "leading-edge genes shown per direction."
        ),
    )

    parser.add_argument("--leading_edge_ribosomal_genes", type=int, default=5,
        help=(
            "Maximum number of ribosomal recurrent "
            "leading-edge genes shown per direction."
        ),
    )

    return parser.parse_args()


def require_columns(df, columns, source):
    missing = [c for c in columns if c not in df.columns]
    if missing:
        raise ValueError(f"Missing required columns {missing} in {source}")


def html_table(df, columns, n):
    columns = [c for c in columns if c in df.columns]
    x = df.loc[:, columns].head(n).copy()
    if x.empty:
        return '<p class="empty">No results met the reporting criteria.</p>'

    for col in x.columns:
        if pd.api.types.is_numeric_dtype(x[col]):
            x[col] = x[col].map(
                lambda v: "" if pd.isna(v) else f"{v:.3g}"
            )

    return x.to_html(
        index=False,
        escape=True,
        border=0,
        classes="data",
    )


def median_abs_delta_from_reference(df, analyses, reference, metric):
    ref_col = f"{reference}_{metric}"
    if ref_col not in df.columns:
        return pd.Series(np.nan, index=df.index)

    differences = []
    for analysis in analyses:
        if analysis == reference:
            continue
        col = f"{analysis}_{metric}"
        if col in df.columns:
            differences.append(
                (numeric(df[col]) - numeric(df[ref_col])).abs()
            )

    if not differences:
        return pd.Series(np.nan, index=df.index)

    return pd.concat(differences, axis=1).median(
        axis=1,
        skipna=True,
    )


# -----------------------------------------------------------------------------
# Gene interpretation
# -----------------------------------------------------------------------------

def gene_pair_evidence(df, analysis, reference, args):
    """Return pairwise direction/evidence columns for one secondary analysis."""

    ref_lfc = numeric(df[f"{reference}_log2FoldChange"])
    sec_lfc = numeric(df[f"{analysis}_log2FoldChange"])

    ref_fdr = numeric(df[f"{reference}_FDR"])
    sec_fdr = numeric(df[f"{analysis}_FDR"])

    ref_p = numeric(df[f"{reference}_pval"])
    sec_p = numeric(df[f"{analysis}_pval"])

    ref_wald = numeric(df[f"{reference}_Wald_stat"])
    sec_wald = numeric(df[f"{analysis}_Wald_stat"])

    available = ref_lfc.notna() & sec_lfc.notna()
    same_direction = available & (np.sign(ref_lfc) == np.sign(sec_lfc))
    opposite_direction = available & (np.sign(ref_lfc) != np.sign(sec_lfc))

    sec_meaningful = (
        sec_lfc.abs().ge(args.min_abs_log2fc)
        & sec_wald.abs().ge(args.min_abs_wald)
        & (sec_fdr.lt(args.fdr) | sec_p.lt(args.nominal_p))
    )

    ref_meaningful = (
        ref_lfc.abs().ge(args.min_abs_log2fc)
        & ref_wald.abs().ge(args.min_abs_wald)
        & (ref_fdr.lt(args.fdr) | ref_p.lt(args.nominal_p))
    )

    delta = (sec_lfc - ref_lfc).abs()

    meaningful_discordance = (
        opposite_direction
        & ref_lfc.abs().ge(args.min_abs_log2fc)
        & sec_lfc.abs().ge(args.min_abs_log2fc)
        & delta.ge(args.min_delta_log2fc)
        & (ref_meaningful | sec_meaningful)
    )

    secondary_fdr_hit = (
        sec_fdr.lt(args.fdr)
        & sec_lfc.abs().ge(args.min_abs_log2fc)
    )

    supported_subset = secondary_fdr_hit & same_direction
    potential_subset_specific = (
        secondary_fdr_hit
        & (
            ref_lfc.abs().lt(args.min_abs_log2fc / 2)
            | opposite_direction
        )
    )

    return {
        "available": available,
        "same_direction": same_direction,
        "opposite_direction": opposite_direction,
        "delta": delta,
        "meaningful_discordance": meaningful_discordance,
        "secondary_fdr_hit": secondary_fdr_hit,
        "supported_subset": supported_subset,
        "potential_subset_specific": potential_subset_specific,
    }


def classify_gene_table(df, analyses, reference, args):
    out = df.copy()

    required = [
        f"{reference}_log2FoldChange",
        f"{reference}_lfcSE",
        f"{reference}_Wald_stat",
        f"{reference}_pval",
        f"{reference}_FDR",
    ]
    require_columns(out, required, "DE comparison table")

    lfc_cols = [
        f"{analysis}_log2FoldChange"
        for analysis in analyses
        if f"{analysis}_log2FoldChange" in out.columns
    ]

    fdr_cols = [
        f"{analysis}_FDR"
        for analysis in analyses
        if f"{analysis}_FDR" in out.columns
    ]

    out["n_available_effects"] = out[lfc_cols].notna().sum(axis=1)
    out["n_significant"] = sum(
        numeric(out[col]).lt(args.fdr)
        for col in fdr_cols
    )

    lfc_matrix = out[lfc_cols].apply(numeric)
    out["n_positive"] = np.sign(lfc_matrix).gt(0).sum(axis=1)
    out["n_negative"] = np.sign(lfc_matrix).lt(0).sum(axis=1)
    out["mean_abs_log2FC"] = lfc_matrix.abs().mean(axis=1)
    out["effect_range"] = (
        lfc_matrix.max(axis=1) - lfc_matrix.min(axis=1)
    )
    out["median_abs_delta_from_reference"] = (
        median_abs_delta_from_reference(
            out,
            analyses,
            reference,
            "log2FoldChange",
        )
    )

    ref_lfc = numeric(out[f"{reference}_log2FoldChange"])
    ref_fdr = numeric(out[f"{reference}_FDR"])
    ref_p = numeric(out[f"{reference}_pval"])
    ref_wald = numeric(out[f"{reference}_Wald_stat"])

    out["reference_significant"] = ref_fdr.lt(args.fdr)
    out["reference_reportable_effect"] = (
        ref_lfc.abs().ge(args.min_abs_log2fc)
    )
    out["reference_nominal_support"] = (
        ref_p.lt(args.nominal_p)
        & ref_wald.abs().ge(args.min_abs_wald)
    )

    supporting = pd.DataFrame(index=out.index)
    discordant = pd.DataFrame(index=out.index)
    supported_subset = pd.DataFrame(index=out.index)
    potential_subset = pd.DataFrame(index=out.index)

    for analysis in analyses:
        if analysis == reference:
            continue

        required_pair = [
            f"{analysis}_log2FoldChange",
            f"{analysis}_Wald_stat",
            f"{analysis}_pval",
            f"{analysis}_FDR",
        ]
        if not all(col in out.columns for col in required_pair):
            continue

        evidence = gene_pair_evidence(
            out,
            analysis,
            reference,
            args,
        )

        out[f"{analysis}_delta_log2FC"] = evidence["delta"]
        out[f"{analysis}_same_direction_as_reference"] = (
            evidence["same_direction"]
        )
        out[f"{analysis}_meaningful_discordance"] = (
            evidence["meaningful_discordance"]
        )

        supporting[analysis] = evidence["same_direction"]
        discordant[analysis] = evidence["meaningful_discordance"]
        supported_subset[analysis] = evidence["supported_subset"]
        potential_subset[analysis] = evidence["potential_subset_specific"]

    out["n_supporting_secondary"] = (
        supporting.sum(axis=1) if not supporting.empty else 0
    )
    out["n_meaningfully_discordant_secondary"] = (
        discordant.sum(axis=1) if not discordant.empty else 0
    )
    out["n_supported_subset_hits"] = (
        supported_subset.sum(axis=1) if not supported_subset.empty else 0
    )
    out["n_potential_subset_specific_hits"] = (
        potential_subset.sum(axis=1) if not potential_subset.empty else 0
    )

    out["report_category"] = "OTHER"

    ref_concordant = (
        out["reference_significant"]
        & out["reference_reportable_effect"]
        & out["n_supporting_secondary"].ge(args.min_supporting_analyses)
        & out["n_meaningfully_discordant_secondary"].eq(0)
    )

    meaningful_discordance = (
        out["n_meaningfully_discordant_secondary"].gt(0)
    )

    potential_subset_specific = (
        ~out["reference_significant"]
        & out["n_potential_subset_specific_hits"].gt(0)
    )

    supported_subset_hit = (
        ~out["reference_significant"]
        & ~potential_subset_specific
        & out["n_supported_subset_hits"].gt(0)
    )

    out.loc[ref_concordant, "report_category"] = (
        "REFERENCE_HIT_CONCORDANT"
    )
    out.loc[meaningful_discordance, "report_category"] = (
        "MEANINGFUL_DISCORDANCE"
    )
    out.loc[supported_subset_hit, "report_category"] = (
        "SUPPORTED_SUBSET_HIT"
    )
    out.loc[potential_subset_specific, "report_category"] = (
        "POTENTIAL_SUBSET_SPECIFIC"
    )

    # Ranking score is only for ordering report tables, not inference.
    out["report_priority"] = (
        out["reference_significant"].astype(int) * 100
        + out["n_supporting_secondary"] * 10
        + out["n_significant"] * 5
        + out["mean_abs_log2FC"].fillna(0) * 2
        + out["n_meaningfully_discordant_secondary"] * 20
        - out["median_abs_delta_from_reference"].fillna(0)
    )

    return out


# -----------------------------------------------------------------------------
# Pathway interpretation
# -----------------------------------------------------------------------------

def pathway_pair_evidence(df, analysis, reference, args):
    ref_nes = numeric(df[f"{reference}_NES"])
    sec_nes = numeric(df[f"{analysis}_NES"])
    ref_fdr = numeric(df[f"{reference}_FDR q-val"])
    sec_fdr = numeric(df[f"{analysis}_FDR q-val"])

    available = ref_nes.notna() & sec_nes.notna()
    same_direction = available & (np.sign(ref_nes) == np.sign(sec_nes))
    opposite_direction = available & (np.sign(ref_nes) != np.sign(sec_nes))
    delta = (sec_nes - ref_nes).abs()

    meaningful_discordance = (
        opposite_direction
        & ref_nes.abs().ge(args.min_abs_nes)
        & sec_nes.abs().ge(args.min_abs_nes)
        & delta.ge(args.min_delta_nes)
        & (ref_fdr.lt(args.fdr) | sec_fdr.lt(args.fdr))
    )

    secondary_fdr_hit = (
        sec_fdr.lt(args.fdr)
        & sec_nes.abs().ge(args.min_abs_nes)
    )

    supported_subset = secondary_fdr_hit & same_direction
    potential_subset_specific = (
        secondary_fdr_hit
        & (
            ref_nes.abs().lt(args.min_abs_nes / 2)
            | opposite_direction
        )
    )

    return {
        "available": available,
        "same_direction": same_direction,
        "opposite_direction": opposite_direction,
        "delta": delta,
        "meaningful_discordance": meaningful_discordance,
        "secondary_fdr_hit": secondary_fdr_hit,
        "supported_subset": supported_subset,
        "potential_subset_specific": potential_subset_specific,
    }


def classify_pathway_table(df, analyses, reference, args):
    out = df.copy()

    require_columns(
        out,
        [f"{reference}_NES", f"{reference}_FDR q-val"],
        "GSEA comparison table",
    )

    nes_cols = [
        f"{analysis}_NES"
        for analysis in analyses
        if f"{analysis}_NES" in out.columns
    ]
    fdr_cols = [
        f"{analysis}_FDR q-val"
        for analysis in analyses
        if f"{analysis}_FDR q-val" in out.columns
    ]

    nes_matrix = out[nes_cols].apply(numeric)
    out["n_available_NES"] = nes_matrix.notna().sum(axis=1)
    out["n_significant"] = sum(
        numeric(out[col]).lt(args.fdr)
        for col in fdr_cols
    )
    out["n_positive"] = np.sign(nes_matrix).gt(0).sum(axis=1)
    out["n_negative"] = np.sign(nes_matrix).lt(0).sum(axis=1)
    out["mean_abs_NES"] = nes_matrix.abs().mean(axis=1)
    out["NES_range"] = (
        nes_matrix.max(axis=1) - nes_matrix.min(axis=1)
    )
    out["median_abs_delta_from_reference"] = (
        median_abs_delta_from_reference(
            out,
            analyses,
            reference,
            "NES",
        )
    )

    ref_nes = numeric(out[f"{reference}_NES"])
    ref_fdr = numeric(out[f"{reference}_FDR q-val"])

    out["reference_significant"] = ref_fdr.lt(args.fdr)
    out["reference_reportable_NES"] = (
        ref_nes.abs().ge(args.min_abs_nes)
    )

    supporting = pd.DataFrame(index=out.index)
    discordant = pd.DataFrame(index=out.index)
    supported_subset = pd.DataFrame(index=out.index)
    potential_subset = pd.DataFrame(index=out.index)

    for analysis in analyses:
        if analysis == reference:
            continue
        if not all(
            col in out.columns
            for col in [
                f"{analysis}_NES",
                f"{analysis}_FDR q-val",
            ]
        ):
            continue

        evidence = pathway_pair_evidence(
            out,
            analysis,
            reference,
            args,
        )

        out[f"{analysis}_delta_NES"] = evidence["delta"]
        out[f"{analysis}_same_direction_as_reference"] = (
            evidence["same_direction"]
        )
        out[f"{analysis}_meaningful_discordance"] = (
            evidence["meaningful_discordance"]
        )

        supporting[analysis] = evidence["same_direction"]
        discordant[analysis] = evidence["meaningful_discordance"]
        supported_subset[analysis] = evidence["supported_subset"]
        potential_subset[analysis] = evidence["potential_subset_specific"]

    out["n_supporting_secondary"] = (
        supporting.sum(axis=1) if not supporting.empty else 0
    )
    out["n_meaningfully_discordant_secondary"] = (
        discordant.sum(axis=1) if not discordant.empty else 0
    )
    out["n_supported_subset_hits"] = (
        supported_subset.sum(axis=1) if not supported_subset.empty else 0
    )
    out["n_potential_subset_specific_hits"] = (
        potential_subset.sum(axis=1) if not potential_subset.empty else 0
    )

    out["report_category"] = "OTHER"

    ref_concordant = (
        out["reference_significant"]
        & out["reference_reportable_NES"]
        & out["n_supporting_secondary"].ge(args.min_supporting_analyses)
        & out["n_meaningfully_discordant_secondary"].eq(0)
    )

    meaningful_discordance = (
        out["n_meaningfully_discordant_secondary"].gt(0)
    )

    potential_subset_specific = (
        ~out["reference_significant"]
        & out["n_potential_subset_specific_hits"].gt(0)
    )

    supported_subset_hit = (
        ~out["reference_significant"]
        & ~potential_subset_specific
        & out["n_supported_subset_hits"].gt(0)
    )

    out.loc[ref_concordant, "report_category"] = (
        "REFERENCE_HIT_CONCORDANT"
    )
    out.loc[meaningful_discordance, "report_category"] = (
        "MEANINGFUL_DISCORDANCE"
    )
    out.loc[supported_subset_hit, "report_category"] = (
        "SUPPORTED_SUBSET_HIT"
    )
    out.loc[potential_subset_specific, "report_category"] = (
        "POTENTIAL_SUBSET_SPECIFIC"
    )

    out["report_priority"] = (
        out["reference_significant"].astype(int) * 100
        + out["n_supporting_secondary"] * 10
        + out["n_significant"] * 5
        + out["mean_abs_NES"].fillna(0) * 2
        + out["n_meaningfully_discordant_secondary"] * 20
        - out["median_abs_delta_from_reference"].fillna(0)
    )

    return out


def parse_leading_edge_genes(value):
    """
    Parse a GSEApy Lead_genes field into a set of genes.

    GSEApy commonly stores leading-edge genes separated by
    semicolons. Commas are accepted defensively as well.
    """
    if pd.isna(value):
        return set()

    value = str(value).strip()

    if not value:
        return set()

    genes = re.split(
        r"[;,]",
        value,
    )

    return {
        gene.strip()
        for gene in genes
        if gene.strip()
    }


def synthesize_leading_edges(
    pathways,
    analyses,
    reference_analysis,
    args,
):
    """
    Synthesize recurrent leading-edge genes across concordant
    significant pathways and analyses.

    Only pathways classified as REFERENCE_HIT_CONCORDANT are
    considered.

    For each pathway, an analysis contributes leading-edge genes
    only when:
      1. the pathway is FDR-significant in that analysis;
      2. |NES| >= args.min_abs_nes;
      3. the NES direction matches the reference analysis.

    Returns
    -------
    occurrences : DataFrame
        One row per gene x pathway x analysis occurrence.

    recurrent : DataFrame
        One row per cell type x gene summarising recurrence.

    pathway_summary : DataFrame
        One row per cell type x pathway summarising leading-edge
        overlap across analyses.
    """

    occurrence_rows = []
    pathway_rows = []

    if pathways.empty:
        return (
            pd.DataFrame(),
            pd.DataFrame(),
            pd.DataFrame(),
        )

    selected = pathways[
        pathways["report_category"]
        == "REFERENCE_HIT_CONCORDANT"
    ].copy()

    if selected.empty:
        return (
            pd.DataFrame(),
            pd.DataFrame(),
            pd.DataFrame(),
        )

    for _, row in selected.iterrows():

        cell_type = row["cell_type"]
        gene_set = row["gene_set"]
        term = row["Term"]

        ref_nes_col = (
            f"{reference_analysis}_NES"
        )

        if (
            ref_nes_col not in row.index
            or pd.isna(row[ref_nes_col])
        ):
            continue

        reference_nes = float(
            row[ref_nes_col]
        )

        reference_direction = int(
            np.sign(reference_nes)
        )

        if reference_direction == 0:
            continue

        genes_by_analysis = {}

        for analysis in analyses:

            nes_col = (
                f"{analysis}_NES"
            )

            fdr_col = (
                f"{analysis}_FDR q-val"
            )

            lead_col = (
                f"{analysis}_Lead_genes"
            )

            required_columns = [
                nes_col,
                fdr_col,
                lead_col,
            ]

            if not all(
                col in row.index
                for col in required_columns
            ):
                continue

            nes = pd.to_numeric(
                row[nes_col],
                errors="coerce",
            )

            fdr = pd.to_numeric(
                row[fdr_col],
                errors="coerce",
            )

            if (
                pd.isna(nes)
                or pd.isna(fdr)
            ):
                continue

            same_direction = (
                np.sign(nes)
                == reference_direction
            )

            significant = (
                fdr < args.fdr
            )

            strong_enrichment = (
                abs(nes)
                >= args.min_abs_nes
            )

            if not (
                same_direction
                and significant
                and strong_enrichment
            ):
                continue

            genes = parse_leading_edge_genes(
                row[lead_col]
            )

            if not genes:
                continue

            genes_by_analysis[
                analysis
            ] = genes

            for gene in genes:

                occurrence_rows.append(
                    {
                        "cell_type": (
                            cell_type
                        ),
                        "gene_set": gene_set,
                        "Term": term,
                        "direction": (
                            "POSITIVE"
                            if reference_direction > 0
                            else "NEGATIVE"
                        ),
                        "analysis": analysis,
                        "gene": gene,
                        "NES": float(nes),
                        "FDR": float(fdr),
                    }
                )

        if not genes_by_analysis:
            continue

        all_genes = set().union(
            *genes_by_analysis.values()
        )

        intersection = set.intersection(
            *genes_by_analysis.values()
        )

        n_analyses = len(
            genes_by_analysis
        )

        pathway_rows.append(
            {
                "cell_type": cell_type,
                "gene_set": gene_set,
                "Term": term,
                "direction": (
                    "POSITIVE"
                    if reference_direction > 0
                    else "NEGATIVE"
                ),
                "n_contributing_analyses": (
                    n_analyses
                ),
                "contributing_analyses": (
                    "|".join(
                        sorted(
                            genes_by_analysis
                        )
                    )
                ),
                "n_union_leading_edge_genes": (
                    len(all_genes)
                ),
                "n_shared_leading_edge_genes": (
                    len(intersection)
                ),
                "shared_leading_edge_genes": (
                    "|".join(
                        sorted(intersection)
                    )
                ),
            }
        )

    occurrences = pd.DataFrame(
        occurrence_rows
    )

    pathway_summary = pd.DataFrame(
        pathway_rows
    )

    if occurrences.empty:
        return (
            occurrences,
            pd.DataFrame(),
            pathway_summary,
        )

    # --------------------------------------------------
    # Per-gene recurrence
    # --------------------------------------------------

    grouped = (
        occurrences
        .groupby(
            [
                "cell_type",
                "gene",
                "direction",
            ],
            observed=True,
        )
    )

    recurrent = (
        grouped
        .agg(
            n_occurrences=(
                "gene",
                "size",
            ),
            n_analyses=(
                "analysis",
                "nunique",
            ),
            n_pathways=(
                "Term",
                "nunique",
            ),
            n_gene_set_collections=(
                "gene_set",
                "nunique",
            ),
            mean_abs_NES=(
                "NES",
                lambda x:
                np.mean(
                    np.abs(x)
                ),
            ),
            min_FDR=(
                "FDR",
                "min",
            ),
        )
        .reset_index()
    )

    recurrent[
        "ribosomal_gene"
    ] = (
        recurrent["gene"]
        .map(
            is_ribosomal_gene
        )
    )

    # --------------------------------------------------
    # Record exactly which analyses/pathways support
    # each recurrent gene.
    # --------------------------------------------------

    analyses_used = (
        grouped["analysis"]
        .agg(
            lambda x:
            "|".join(
                sorted(
                    set(x)
                )
            )
        )
        .reset_index(
            name="analyses"
        )
    )

    pathways_used = (
        grouped["Term"]
        .agg(
            lambda x:
            "|".join(
                sorted(
                    set(x)
                )
            )
        )
        .reset_index(
            name="pathways"
        )
    )

    collections_used = (
        grouped["gene_set"]
        .agg(
            lambda x:
            "|".join(
                sorted(
                    set(x)
                )
            )
        )
        .reset_index(
            name="gene_set_collections"
        )
    )

    recurrent = recurrent.merge(
        analyses_used,
        on=[
            "cell_type",
            "gene",
            "direction",
        ],
        how="left",
    )

    recurrent = recurrent.merge(
        pathways_used,
        on=[
            "cell_type",
            "gene",
            "direction",
        ],
        how="left",
    )

    recurrent = recurrent.merge(
        collections_used,
        on=[
            "cell_type",
            "gene",
            "direction",
        ],
        how="left",
    )

    # --------------------------------------------------
    # Reporting filter
    # --------------------------------------------------

    recurrent["reportable"] = (
        (
            recurrent["n_analyses"]
            >=
            args.leading_edge_min_analyses
        )
        &
        (
            recurrent["n_pathways"]
            >=
            args.leading_edge_min_pathways
        )
    )

    recurrent = recurrent.sort_values(
        [
            "cell_type",
            "reportable",
            "n_analyses",
            "n_pathways",
            "n_gene_set_collections",
            "n_occurrences",
            "mean_abs_NES",
        ],
        ascending=[
            True,
            False,
            False,
            False,
            False,
            False,
            False,
        ],
    )

    return (
        occurrences,
        recurrent,
        pathway_summary,
    )


def is_ribosomal_gene(gene):
    """
    Identify canonical ribosomal protein gene symbols.

    RPL* and RPS* genes are retained in the analysis but can
    be limited separately in interpretation plots so that a
    broad ribosomal programme does not dominate display space.
    """
    gene = str(
        gene
    ).upper()

    return bool(
        re.match(
            r"^RP[LS][0-9A-Z-]+$",
            gene,
        )
    )

# -----------------------------------------------------------------------------
# Plotting
# -----------------------------------------------------------------------------

def plot_de_concordance(df, cell_type, analyses, reference, outdir, args):
    ref_col = f"{reference}_log2FoldChange"
    ref_fdr = f"{reference}_FDR"

    for analysis in analyses:
        if analysis == reference:
            continue
        sec_col = f"{analysis}_log2FoldChange"
        if sec_col not in df.columns:
            continue

        plot_df = pd.DataFrame(
            {
                "ref": numeric(df[ref_col]),
                "sec": numeric(df[sec_col]),
                "ref_sig": numeric(df[ref_fdr]).lt(args.fdr),
                "discordant": df.get(
                    f"{analysis}_meaningful_discordance",
                    False,
                ),
            }
        ).dropna(subset=["ref", "sec"])

        if plot_df.empty:
            continue

        fig, ax = plt.subplots(figsize=(6.5, 6.2))
        ax.scatter(
            plot_df["ref"],
            plot_df["sec"],
            s=10,
            alpha=0.20,
            color="0.55",
            label="Other genes",
        )

        sig = plot_df[plot_df["ref_sig"]]
        if not sig.empty:
            ax.scatter(
                sig["ref"],
                sig["sec"],
                s=34,
                alpha=0.85,
                color="#2B6CB0",
                label="Reference FDR < threshold",
            )

        disc = plot_df[plot_df["discordant"]]
        if not disc.empty:
            ax.scatter(
                disc["ref"],
                disc["sec"],
                s=45,
                alpha=0.9,
                color="#C53030",
                marker="x",
                label="Meaningful discordance",
            )

        lim = np.nanmax(
            np.abs(
                pd.concat([plot_df["ref"], plot_df["sec"]])
            )
        )
        lim = max(float(lim), 0.5) * 1.05
        ax.plot([-lim, lim], [-lim, lim], "--", color="0.35", lw=1)
        ax.axhline(0, color="0.75", lw=0.8)
        ax.axvline(0, color="0.75", lw=0.8)
        ax.set_xlim(-lim, lim)
        ax.set_ylim(-lim, lim)
        ax.set_xlabel(f"{reference} log2FC")
        ax.set_ylabel(f"{analysis} log2FC")
        ax.set_title(f"{cell_type}\nDE effect-size concordance")
        ax.legend(fontsize=8, loc="best")
        fig.tight_layout()
        fig.savefig(
            outdir / f"{safe_filename(reference)}_vs_{safe_filename(analysis)}.png",
            dpi=180,
            bbox_inches="tight",
        )
        plt.close(fig)


def plot_nes_concordance(
    df,
    cell_type,
    gene_set,
    analyses,
    reference,
    outdir,
    args,
):
    ref_col = f"{reference}_NES"
    ref_fdr = f"{reference}_FDR q-val"

    for analysis in analyses:
        if analysis == reference:
            continue
        sec_col = f"{analysis}_NES"
        if sec_col not in df.columns:
            continue

        plot_df = pd.DataFrame(
            {
                "ref": numeric(df[ref_col]),
                "sec": numeric(df[sec_col]),
                "ref_sig": numeric(df[ref_fdr]).lt(args.fdr),
                "discordant": df.get(
                    f"{analysis}_meaningful_discordance",
                    False,
                ),
            }
        ).dropna(subset=["ref", "sec"])

        if plot_df.empty:
            continue

        fig, ax = plt.subplots(figsize=(6.5, 6.2))
        ax.scatter(
            plot_df["ref"],
            plot_df["sec"],
            s=14,
            alpha=0.25,
            color="0.55",
            label="Other pathways",
        )

        sig = plot_df[plot_df["ref_sig"]]
        if not sig.empty:
            ax.scatter(
                sig["ref"],
                sig["sec"],
                s=38,
                alpha=0.85,
                color="#2B6CB0",
                label="Reference FDR < threshold",
            )

        disc = plot_df[plot_df["discordant"]]
        if not disc.empty:
            ax.scatter(
                disc["ref"],
                disc["sec"],
                s=48,
                alpha=0.9,
                color="#C53030",
                marker="x",
                label="Meaningful discordance",
            )

        lim = np.nanmax(
            np.abs(pd.concat([plot_df["ref"], plot_df["sec"]]))
        )
        lim = max(float(lim), 1.0) * 1.05
        ax.plot([-lim, lim], [-lim, lim], "--", color="0.35", lw=1)
        ax.axhline(0, color="0.75", lw=0.8)
        ax.axvline(0, color="0.75", lw=0.8)
        ax.set_xlim(-lim, lim)
        ax.set_ylim(-lim, lim)
        ax.set_xlabel(f"{reference} NES")
        ax.set_ylabel(f"{analysis} NES")
        ax.set_title(f"{cell_type} - {gene_set}\nGSEA NES concordance")
        ax.legend(fontsize=8, loc="best")
        fig.tight_layout()
        fig.savefig(
            outdir / f"{safe_filename(reference)}_vs_{safe_filename(analysis)}.png",
            dpi=180,
            bbox_inches="tight",
        )
        plt.close(fig)


def plot_gene_heatmap(df, cell_type, analyses, reference, outpath, args):
    candidate = df[
        df["report_category"].isin(
            [
                "REFERENCE_HIT_CONCORDANT",
                "MEANINGFUL_DISCORDANCE",
                "SUPPORTED_SUBSET_HIT",
                "POTENTIAL_SUBSET_SPECIFIC",
            ]
        )
    ].copy()

    if candidate.empty:
        return

    candidate = candidate.sort_values(
        ["reference_significant", "report_priority"],
        ascending=[False, False],
    ).head(args.heatmap_genes)

    cols = [
        f"{analysis}_log2FoldChange"
        for analysis in analyses
        if f"{analysis}_log2FoldChange" in candidate.columns
    ]

    matrix = candidate.set_index("gene")[cols].apply(numeric)
    matrix.columns = [c.replace("_log2FoldChange", "") for c in matrix.columns]

    if matrix.empty:
        return

    height = max(4.5, 0.28 * len(matrix) + 2)
    fig, ax = plt.subplots(figsize=(7.2, height))
    sns.heatmap(
        matrix,
        cmap="vlag",
        center=0,
        annot=False,
        linewidths=0.25,
        linecolor="white",
        ax=ax,
    )
    ax.set_title(f"{cell_type}\nSelected gene log2FC values")
    ax.set_xlabel("Analysis")
    ax.set_ylabel("Gene")
    fig.tight_layout()
    fig.savefig(outpath, dpi=180, bbox_inches="tight")
    plt.close(fig)


def plot_pathway_heatmap(
    df,
    cell_type,
    gene_set,
    analyses,
    reference,
    outpath,
    args,
):
    candidate = df[
        df["report_category"].isin(
            [
                "REFERENCE_HIT_CONCORDANT",
                "MEANINGFUL_DISCORDANCE",
                "SUPPORTED_SUBSET_HIT",
                "POTENTIAL_SUBSET_SPECIFIC",
            ]
        )
    ].copy()

    if candidate.empty:
        return

    candidate = candidate.sort_values(
        ["reference_significant", "report_priority"],
        ascending=[False, False],
    ).head(args.heatmap_pathways)

    cols = [
        f"{analysis}_NES"
        for analysis in analyses
        if f"{analysis}_NES" in candidate.columns
    ]

    matrix = candidate.set_index("Term")[cols].apply(numeric)
    matrix.columns = [c.replace("_NES", "") for c in matrix.columns]

    if matrix.empty:
        return

    height = max(5, 0.26 * len(matrix) + 2)
    width = 8
    fig, ax = plt.subplots(figsize=(width, height))
    sns.heatmap(
        matrix,
        cmap="vlag",
        center=0,
        annot=False,
        linewidths=0.25,
        linecolor="white",
        ax=ax,
    )
    ax.set_title(f"{cell_type} - {gene_set}\nSelected pathway NES values")
    ax.set_xlabel("Analysis")
    ax.set_ylabel("Pathway")
    fig.tight_layout()
    fig.savefig(outpath, dpi=180, bbox_inches="tight")
    plt.close(fig)

def plot_leading_edge_heatmap(
    occurrences,
    recurrent,
    cell_type,
    direction,
    analyses,
    outpath,
    args,
):
    if (
        occurrences.empty
        or recurrent.empty
    ):
        return

    recurrent_sub = (
        recurrent[
            (
                recurrent["cell_type"]
                == cell_type
            )
            &
            (
                recurrent["direction"]
                == direction
            )
            &
            recurrent["reportable"]
        ]
        .copy()
    )

    if recurrent_sub.empty:
        return

    # --------------------------------------------------
    # Rank recurrent genes
    # --------------------------------------------------

    sort_columns = [
        "n_analyses",
        "n_pathways",
        "n_gene_set_collections",
        "n_occurrences",
        "mean_abs_NES",
    ]

    recurrent_sub = (
        recurrent_sub
        .sort_values(
            sort_columns,
            ascending=False,
        )
    )

    # --------------------------------------------------
    # Select non-ribosomal and ribosomal genes
    # separately so ribosomal programmes remain visible
    # without dominating the display.
    # --------------------------------------------------

    non_ribosomal = (
        recurrent_sub[
            ~recurrent_sub[
                "ribosomal_gene"
            ]
        ]
        .head(
            args.leading_edge_nonribosomal_genes
        )
    )

    ribosomal = (
        recurrent_sub[
            recurrent_sub[
                "ribosomal_gene"
            ]
        ]
        .head(
            args.leading_edge_ribosomal_genes
        )
    )

    selected = pd.concat(
        [
            non_ribosomal,
            ribosomal,
        ],
        ignore_index=True,
    )

    if selected.empty:
        return

    genes = set(
        selected["gene"]
    )

    occurrence_sub = (
        occurrences[
            (
                occurrences["cell_type"]
                == cell_type
            )
            &
            (
                occurrences["direction"]
                == direction
            )
            &
            occurrences["gene"].isin(
                genes
            )
        ]
    )

    if occurrence_sub.empty:
        return

    # --------------------------------------------------
    # Count unique significant concordant pathways
    # containing each gene in each analysis.
    # --------------------------------------------------

    matrix = (
        occurrence_sub
        .groupby(
            [
                "gene",
                "analysis",
            ],
            observed=True,
        )["Term"]
        .nunique()
        .unstack(
            fill_value=0
        )
    )

    matrix = matrix.reindex(
        columns=[
            analysis
            for analysis in analyses
            if analysis in matrix.columns
        ]
    )

    # Preserve selected ranking.
    gene_order = (
        selected["gene"]
        .tolist()
    )

    matrix = matrix.reindex(
        gene_order
    )

    height = max(
        4.5,
        0.30 * len(matrix) + 2,
    )

    fig, ax = plt.subplots(
        figsize=(7.5, height)
    )

    sns.heatmap(
        matrix,
        cmap="Blues",
        annot=True,
        fmt="g",
        linewidths=0.25,
        linecolor="white",
        ax=ax,
    )

    direction_label = (
        "Positive"
        if direction == "POSITIVE"
        else "Negative"
    )

    ax.set_title(
        f"{cell_type}\n"
        f"{direction_label} recurrent leading-edge "
        "pathway membership"
    )

    ax.set_xlabel(
        "Analysis"
    )

    ax.set_ylabel(
        "Gene"
    )

    fig.tight_layout()

    fig.savefig(
        outpath,
        dpi=180,
        bbox_inches="tight",
    )

    plt.close(fig)

# -----------------------------------------------------------------------------
# HTML documentation
# -----------------------------------------------------------------------------

def explanatory_html(args):
    return f"""
    <section class="guide">
      <h2>How to interpret this report</h2>
      <p>
        <b>{html.escape(args.reference_analysis)}</b> is treated as the
        prespecified reference analysis. Secondary analyses are used to assess
        robustness across sample subsets or alternative model specifications.
        Agreement is assessed primarily from the <b>direction and magnitude</b>
        of effects, not by requiring independent FDR significance in every
        analysis, because secondary analyses may have fewer samples and lower
        statistical power.
      </p>
      <p>
        For genes, positive log2FoldChange values indicate higher expression
        in the test condition and negative values indicate lower expression.
        lfcSE is the uncertainty of the log2 fold-change estimate. Wald_stat
        represents the model coefficient relative to its uncertainty. FDR is
        the multiple-testing-adjusted P value used for formal gene-level
        significance.
      </p>
      <p>
        For pathways, positive NES values indicate enrichment toward the
        positive/test-condition end of the ranked list and negative NES values
        indicate enrichment toward the reference-condition end. Pathway FDR is
        used for formal enrichment significance.
      </p>
      <p>
        <b>Recurrent leading-edge genes</b> are genes repeatedly
        contributing to the leading edges of concordant significant
        pathways across analyses. A gene must occur in at least
        {args.leading_edge_min_pathways} qualifying pathways and at
        least {args.leading_edge_min_analyses} analyses to be surfaced
        in the recurrent leading-edge tables. These counts identify
        recurring contributors to pathway-level responses and are not
        additional gene-level significance tests.
      </p>
      <p>
        Recurrent leading-edge genes are displayed separately for
        positive and negative pathway enrichment. Ribosomal protein
        genes are retained in all recurrence calculations and output
        tables, but the number displayed in each heatmap is capped
        separately from non-ribosomal genes. This prevents a broad
        ribosomal programme from using most available display rows
        while preserving the ribosomal signal for interpretation.
      </p>
      <p>
        <b>Concordant reference hits</b> are significant in the reference
        analysis, exceed the report effect-size threshold, and retain the same
        direction in at least {args.min_supporting_analyses} secondary analyses.
        Secondary analyses do not need to be independently FDR-significant.
      </p>
      <p>
        <b>Meaningful discordance</b> requires more than a trivial sign change:
        both compared effects must exceed the configured effect-size threshold,
        their difference must exceed the configured delta threshold, and at
        least one analysis must contain statistical support. This avoids
        highlighting thousands of near-zero direction changes.
      </p>
      <p>
        <b>Supported subset hits</b> are FDR-significant in a secondary analysis
        and point in the same direction as the reference estimate. These may
        reflect power differences. <b>Potential subset-specific hits</b> are
        significant in a secondary analysis while the reference estimate is
        small or points in the opposite direction. Neither category proves a
        statistically different treatment response between systems; formal
        evidence of such a difference requires an interaction test.
      </p>

      <h3>Fields used in interpretation</h3>
      <dl>
        <dt>n_significant</dt><dd>Number of analyses with FDR below the report threshold.</dd>
        <dt>n_available_effects / n_available_NES</dt><dd>Number of analyses contributing an estimated gene effect or pathway NES.</dd>
        <dt>n_supporting_secondary</dt><dd>Number of secondary analyses with the same direction as the reference estimate.</dd>
        <dt>mean_abs_log2FC / mean_abs_NES</dt><dd>Mean absolute effect magnitude across available analyses.</dd>
        <dt>effect_range / NES_range</dt><dd>Maximum minus minimum estimate across analyses.</dd>
        <dt>median_abs_delta_from_reference</dt><dd>Median absolute difference between secondary estimates and the reference estimate.</dd>
        <dt>report_priority</dt><dd>Ordering score used only to rank report tables; it is not a statistical test.</dd>
      </dl>

      <h3>Report thresholds</h3>
      <ul>
        <li>FDR: {args.fdr}</li>
        <li>Nominal P threshold: {args.nominal_p}</li>
        <li>Minimum |log2FC|: {args.min_abs_log2fc}</li>
        <li>Minimum |Wald statistic|: {args.min_abs_wald}</li>
        <li>Minimum |delta log2FC|: {args.min_delta_log2fc}</li>
        <li>Minimum |NES|: {args.min_abs_nes}</li>
        <li>Minimum |delta NES|: {args.min_delta_nes}</li>
      </ul>
    </section>
    """


# -----------------------------------------------------------------------------
# Main
# -----------------------------------------------------------------------------

def main():
    args = parse_args()

    if args.reference_analysis not in args.analyses:
        raise ValueError(
            "--reference_analysis must also be included in --analyses"
        )

    comparison_dir = Path(args.comparison_dir).resolve()
    if not comparison_dir.exists():
        raise FileNotFoundError(
            "Comparison directory does not exist: "
            f"{comparison_dir}"
        )

    # --------------------------------------------------
    # Validate comparison inputs
    # --------------------------------------------------

    de_source = (
        comparison_dir
        / "DE"
        / "by_cell_type"
    )

    gsea_source = (
        comparison_dir
        / "GSEA"
        / "by_cell_type"
    )

    if not de_source.exists():
        raise FileNotFoundError(
            "DE comparison directory not found: "
            f"{de_source}\n"
            "The --comparison_dir argument should point "
            "to the annotation-specific comparison directory "
            "containing DE/ and GSEA/."
        )

    if not gsea_source.exists():
        raise FileNotFoundError(
            "GSEA comparison directory not found: "
            f"{gsea_source}\n"
            "The --comparison_dir argument should point "
            "to the annotation-specific comparison directory "
            "containing DE/ and GSEA/."
        )

    de_files = sorted(
        de_source.glob("*.csv")
    )

    if not de_files:
        raise FileNotFoundError(
            "No DE comparison CSV files found in: "
            f"{de_source}"
        )

    gsea_cell_dirs = sorted(
        path
        for path in gsea_source.iterdir()
        if path.is_dir()
    )

    if not gsea_cell_dirs:
        raise FileNotFoundError(
            "No GSEA cell-type directories found in: "
            f"{gsea_source}"
        )

    report_dir = (
        comparison_dir
        / "reports"
    )

    gene_dir = report_dir / "genes"
    pathway_dir = report_dir / "pathways"
    plot_dir = report_dir / "plots"
    de_plot_dir = plot_dir / "DE_concordance"
    gsea_plot_dir = plot_dir / "GSEA_concordance"
    gene_heatmap_dir = plot_dir / "gene_heatmaps"
    nes_heatmap_dir = plot_dir / "NES_heatmaps"
    leading_edge_dir = (report_dir / "leading_edge")
    leading_edge_by_cell_dir = (leading_edge_dir / "by_cell_type")
    leading_edge_heatmap_dir = plot_dir / "leading_edge_heatmaps"

    for path in [
        report_dir,
        gene_dir,
        pathway_dir,
        plot_dir,
        de_plot_dir,
        gsea_plot_dir,
        gene_heatmap_dir,
        nes_heatmap_dir,
        leading_edge_dir,
        leading_edge_by_cell_dir,
        leading_edge_heatmap_dir
    ]:
        path.mkdir(parents=True, exist_ok=True)

    manifest = {
        "comparison_dir": str(comparison_dir),
        "report_name": args.report_name,
        "reference_analysis": args.reference_analysis,
        "analyses": args.analyses,
        "fdr": args.fdr,
        "nominal_p": args.nominal_p,
        "min_abs_log2fc": args.min_abs_log2fc,
        "min_abs_wald": args.min_abs_wald,
        "min_delta_log2fc": args.min_delta_log2fc,
        "min_abs_nes": args.min_abs_nes,
        "min_delta_nes": args.min_delta_nes,
        "min_supporting_analyses": args.min_supporting_analyses,
        "top_genes": args.top_genes,
        "top_pathways": args.top_pathways,
        "heatmap_genes": args.heatmap_genes,
        "heatmap_pathways": args.heatmap_pathways,
        "make_plots": args.make_plots,
        "leading_edge_min_analyses": args.leading_edge_min_analyses,
        "leading_edge_min_pathways": args.leading_edge_min_pathways,
        "top_leading_edge_genes": args.top_leading_edge_genes,
        "leading_edge_nonribosomal_genes": args.leading_edge_nonribosomal_genes,
        "leading_edge_ribosomal_genes": args.leading_edge_ribosomal_genes,
    }

    with open(report_dir / "report_manifest.json", "w") as handle:
        json.dump(manifest, handle, indent=2)

    # ------------------------------------------------------------------
    # Genes
    # ------------------------------------------------------------------
    gene_frames = []

    for path in de_files:
        df = classify_gene_table(
            pd.read_csv(path),
            args.analyses,
            args.reference_analysis,
            args,
        )
        df.insert(0, "cell_type", path.stem)
        gene_frames.append(df)

    genes = (
        pd.concat(gene_frames, ignore_index=True, sort=False)
        if gene_frames
        else pd.DataFrame()
    )

    if not genes.empty:
        genes.to_csv(
            report_dir / "gene_interpretation_master.csv",
            index=False,
        )

        gene_outputs = {
            "reference_concordant.csv": "REFERENCE_HIT_CONCORDANT",
            "meaningful_discordance.csv": "MEANINGFUL_DISCORDANCE",
            "supported_subset_hits.csv": "SUPPORTED_SUBSET_HIT",
            "potential_subset_specific.csv": "POTENTIAL_SUBSET_SPECIFIC",
        }

        for filename, category in gene_outputs.items():
            subset = genes[
                genes["report_category"] == category
            ].sort_values(
                "report_priority",
                ascending=False,
            )
            subset.to_csv(gene_dir / filename, index=False)

    # ------------------------------------------------------------------
    # Pathways
    # ------------------------------------------------------------------
    pathway_frames = []

    for cell_dir in gsea_cell_dirs:

        for path in sorted(
            cell_dir.glob("*.csv")
        ):

            df = classify_pathway_table(
                pd.read_csv(path),
                args.analyses,
                args.reference_analysis,
                args,
            )

            df.insert(
                0,
                "gene_set",
                path.stem,
            )

            df.insert(
                0,
                "cell_type",
                cell_dir.name,
            )

            pathway_frames.append(
                df
            )

    pathways = (
        pd.concat(pathway_frames, ignore_index=True, sort=False)
        if pathway_frames
        else pd.DataFrame()
    )

    if not pathways.empty:
        pathways.to_csv(
            report_dir / "pathway_interpretation_master.csv",
            index=False,
        )

        pathway_outputs = {
            "reference_concordant.csv": "REFERENCE_HIT_CONCORDANT",
            "meaningful_discordance.csv": "MEANINGFUL_DISCORDANCE",
            "supported_subset_hits.csv": "SUPPORTED_SUBSET_HIT",
            "potential_subset_specific.csv": "POTENTIAL_SUBSET_SPECIFIC",
        }

        for filename, category in pathway_outputs.items():
            subset = pathways[pathways["report_category"] == category].sort_values(
                "report_priority", ascending=False,)
            subset.to_csv(pathway_dir / filename, index=False)

    # ------------------------------------------------------------------
    # Cross-analysis leading-edge synthesis
    # ------------------------------------------------------------------

    (
        leading_edge_occurrences,
        recurrent_leading_edge,
        leading_edge_pathways,
    ) = synthesize_leading_edges(
        pathways=pathways,
        analyses=args.analyses,
        reference_analysis=(
            args.reference_analysis
        ),
        args=args,
    )

    if not leading_edge_occurrences.empty:

        leading_edge_occurrences.to_csv(
            leading_edge_dir
            / "leading_edge_occurrences.csv",
            index=False,
        )

    if not leading_edge_pathways.empty:

        leading_edge_pathways.to_csv(
            leading_edge_dir
            / "concordant_pathway_leading_edge_summary.csv",
            index=False,
        )

    if not recurrent_leading_edge.empty:

        recurrent_leading_edge.to_csv(
            leading_edge_dir
            / "recurrent_genes_all.csv",
            index=False,
        )

        recurrent_positive = (
            recurrent_leading_edge[
                (
                    recurrent_leading_edge[
                        "direction"
                    ]
                    == "POSITIVE"
                )
                &
                recurrent_leading_edge[
                    "reportable"
                ]
            ]
            .copy()
        )

        recurrent_negative = (
            recurrent_leading_edge[
                (
                    recurrent_leading_edge[
                        "direction"
                    ]
                    == "NEGATIVE"
                )
                &
                recurrent_leading_edge[
                    "reportable"
                ]
            ]
            .copy()
        )

        recurrent_positive.to_csv(
            leading_edge_dir
            / "recurrent_genes_positive.csv",
            index=False,
        )

        recurrent_negative.to_csv(
            leading_edge_dir
            / "recurrent_genes_negative.csv",
            index=False,
        )

        # --------------------------------------------------
        # Per-cell-type outputs
        # --------------------------------------------------

        for cell_type, sub in (
            recurrent_leading_edge
            .groupby(
                "cell_type",
                sort=True,
            )
        ):

            sub = sub[
                sub["reportable"]
            ].copy()

            sub.to_csv(
                leading_edge_by_cell_dir
                / (
                    f"{safe_filename(cell_type)}.csv"
                ),
                index=False,
            )

    else:

        recurrent_positive = (
            pd.DataFrame()
        )

        recurrent_negative = (
            pd.DataFrame()
        )

    # ------------------------------------------------------------------
    # Cell-type summary
    # ------------------------------------------------------------------
    cells = sorted(
        set(genes.get("cell_type", []))
        | set(pathways.get("cell_type", []))
    )

    summary_rows = []

    for cell_type in cells:
        gene_sub = (
            genes[genes["cell_type"] == cell_type]
            if not genes.empty
            else pd.DataFrame()
        )
        pathway_sub = (
            pathways[pathways["cell_type"] == cell_type]
            if not pathways.empty
            else pd.DataFrame()
        )
        leading_edge_sub = (
            recurrent_leading_edge[
                (
                    recurrent_leading_edge[
                        "cell_type"
                    ]
                    == cell_type
                )
                &
                recurrent_leading_edge[
                    "reportable"
                ]
            ]
            if not recurrent_leading_edge.empty
            else pd.DataFrame()
        )

        summary_rows.append(
            {
                "cell_type": cell_type,
                "reference_sig_genes": int(
                    gene_sub.get(
                        "reference_significant",
                        pd.Series(dtype=bool),
                    ).sum()
                ),
                "reportable_concordant_reference_genes": int(
                    (
                        gene_sub.get(
                            "report_category",
                            pd.Series(dtype=str),
                        )
                        == "REFERENCE_HIT_CONCORDANT"
                    ).sum()
                ),
                "meaningfully_discordant_genes": int(
                    (
                        gene_sub.get(
                            "report_category",
                            pd.Series(dtype=str),
                        )
                        == "MEANINGFUL_DISCORDANCE"
                    ).sum()
                ),
                "supported_subset_genes": int(
                    (
                        gene_sub.get(
                            "report_category",
                            pd.Series(dtype=str),
                        )
                        == "SUPPORTED_SUBSET_HIT"
                    ).sum()
                ),
                "potential_subset_specific_genes": int(
                    (
                        gene_sub.get(
                            "report_category",
                            pd.Series(dtype=str),
                        )
                        == "POTENTIAL_SUBSET_SPECIFIC"
                    ).sum()
                ),
                "reference_sig_pathways": int(
                    pathway_sub.get(
                        "reference_significant",
                        pd.Series(dtype=bool),
                    ).sum()
                ),
                "reportable_concordant_reference_pathways": int(
                    (
                        pathway_sub.get(
                            "report_category",
                            pd.Series(dtype=str),
                        )
                        == "REFERENCE_HIT_CONCORDANT"
                    ).sum()
                ),
                "meaningfully_discordant_pathways": int(
                    (
                        pathway_sub.get(
                            "report_category",
                            pd.Series(dtype=str),
                        )
                        == "MEANINGFUL_DISCORDANCE"
                    ).sum()
                ),
                "supported_subset_pathways": int(
                    (
                        pathway_sub.get(
                            "report_category",
                            pd.Series(dtype=str),
                        )
                        == "SUPPORTED_SUBSET_HIT"
                    ).sum()
                ),
                "potential_subset_specific_pathways": int(
                    (
                        pathway_sub.get(
                            "report_category",
                            pd.Series(dtype=str),
                        )
                        == "POTENTIAL_SUBSET_SPECIFIC"
                    ).sum()
                ),
                "recurrent_leading_edge_genes": (
                    len(
                        leading_edge_sub
                    )
                ),

                "recurrent_positive_leading_edge_genes": int(
                    (
                        leading_edge_sub[
                            "direction"
                        ]
                        == "POSITIVE"
                    ).sum()
                    if not leading_edge_sub.empty
                    else 0
                ),

                "recurrent_negative_leading_edge_genes": int(
                    (
                        leading_edge_sub[
                            "direction"
                        ]
                        == "NEGATIVE"
                    ).sum()
                    if not leading_edge_sub.empty
                    else 0
                ),
            }
        )

    summary = pd.DataFrame(summary_rows)
    summary.to_csv(
        report_dir / "cell_type_interpretation_summary.csv",
        index=False,
    )

    # ------------------------------------------------------------------
    # Plots
    # ------------------------------------------------------------------
    if args.make_plots:
        for cell_type in cells:
            gene_sub = (
                genes[genes["cell_type"] == cell_type].copy()
                if not genes.empty
                else pd.DataFrame()
            )

            if not gene_sub.empty:
                cell_plot_dir = de_plot_dir / safe_filename(cell_type)
                cell_plot_dir.mkdir(parents=True, exist_ok=True)

                plot_de_concordance(
                    gene_sub,
                    cell_type,
                    args.analyses,
                    args.reference_analysis,
                    cell_plot_dir,
                    args,
                )

                plot_gene_heatmap(
                    gene_sub,
                    cell_type,
                    args.analyses,
                    args.reference_analysis,
                    gene_heatmap_dir
                    / f"{safe_filename(cell_type)}.png",
                    args,
                )

            pathway_sub = (
                pathways[pathways["cell_type"] == cell_type].copy()
                if not pathways.empty
                else pd.DataFrame()
            )

            if not pathway_sub.empty:
                for gene_set, gs_df in pathway_sub.groupby(
                    "gene_set",
                    sort=True,
                ):
                    gs_plot_dir = (
                        gsea_plot_dir
                        / safe_filename(cell_type)
                        / safe_filename(gene_set)
                    )
                    gs_plot_dir.mkdir(parents=True, exist_ok=True)

                    plot_nes_concordance(
                        gs_df,
                        cell_type,
                        gene_set,
                        args.analyses,
                        args.reference_analysis,
                        gs_plot_dir,
                        args,
                    )

                    # Heatmaps for every collection are capped and therefore
                    # remain interpretable
                    plot_pathway_heatmap(
                        gs_df,
                        cell_type,
                        gene_set,
                        args.analyses,
                        args.reference_analysis,
                        nes_heatmap_dir
                        / (
                            f"{safe_filename(cell_type)}__"
                            f"{safe_filename(gene_set)}.png"
                        ),
                        args,
                    )

            if (
                not leading_edge_occurrences.empty
                and
                not recurrent_leading_edge.empty
            ):

                for direction in [
                    "POSITIVE",
                    "NEGATIVE",
                ]:

                    plot_leading_edge_heatmap(
                        occurrences=(
                            leading_edge_occurrences
                        ),
                        recurrent=(
                            recurrent_leading_edge
                        ),
                        cell_type=cell_type,
                        direction=direction,
                        analyses=args.analyses,
                        outpath=(
                            leading_edge_heatmap_dir
                            / (
                                f"{safe_filename(cell_type)}__"
                                f"{direction.lower()}.png"
                            )
                        ),
                        args=args,
                    )

    # ------------------------------------------------------------------
    # HTML report
    # ------------------------------------------------------------------
    css = """
    <style>
      body {font-family: Arial, sans-serif; color:#263238; margin:34px;
            max-width:1500px; line-height:1.45;}
      h1,h2 {color:#244a73;} h3 {color:#334e68; border-bottom:1px solid #d9e2ec;
            padding-bottom:4px; margin-top:26px;}
      .note,.guide {background:#eef4f8; padding:14px 18px;
            border-left:4px solid #4472a3; margin:14px 0 24px;}
      .data {border-collapse:collapse; font-size:12px; width:100%; margin:8px 0 22px;}
      .data th {background:#244a73; color:white; padding:6px; position:sticky; top:0;}
      .data td {border:1px solid #d9e2ec; padding:5px;}
      .data tr:nth-child(even) {background:#f7f9fb;}
      .empty {color:#66788a; font-style:italic;}
      .plot-grid {display:grid; grid-template-columns:repeat(auto-fit,minmax(420px,1fr)); gap:16px;}
      .plot-card {background: #ffffff; border: 1px solid #d9e2ec; padding: 10px;}
      .plot-card h4 {color: #334e68; margin: 0 0 8px 0; font-size: 14px;}
      .plot-card img {width: 100%; height: auto; border: none;}
      dt {font-weight:bold; color:#334e68; margin-top:6px;}
      dd {margin-bottom:5px;}
    </style>
    """

    html_parts = [
        "<html><head><meta charset='utf-8'>",
        f"<title>{html.escape(args.report_name)}</title>",
        css,
        "</head><body>",
        f"<h1>{html.escape(args.report_name)}</h1>",
        (
            "<div class='note'>"
            f"<b>Reference analysis:</b> {html.escape(args.reference_analysis)}<br>"
            f"<b>Compared analyses:</b> {html.escape(', '.join(args.analyses))}"
            "</div>"
        ),
        explanatory_html(args),
    ]

    if not summary.empty:
        html_parts.extend(
            [
                "<h2>Cell-type overview</h2>",
                html_table(
                    summary,
                    list(summary.columns),
                    len(summary),
                ),
            ]
        )

    gene_columns = [
        "gene",
        f"{args.reference_analysis}_log2FoldChange",
        f"{args.reference_analysis}_lfcSE",
        f"{args.reference_analysis}_FDR",
    ]
    for analysis in args.analyses:
        if analysis != args.reference_analysis:
            gene_columns.extend(
                [
                    f"{analysis}_log2FoldChange",
                    f"{analysis}_FDR",
                ]
            )
    gene_columns.extend(
        [
            "n_significant",
            "n_supporting_secondary",
            "median_abs_delta_from_reference",
            "report_category",
        ]
    )

    pathway_columns = [
        "Term",
        "gene_set",
        f"{args.reference_analysis}_NES",
        f"{args.reference_analysis}_FDR q-val",
    ]
    for analysis in args.analyses:
        if analysis != args.reference_analysis:
            pathway_columns.extend(
                [
                    f"{analysis}_NES",
                    f"{analysis}_FDR q-val",
                ]
            )
    pathway_columns.extend(
        [
            "n_significant",
            "n_supporting_secondary",
            "median_abs_delta_from_reference",
            "report_category",
        ]
    )

    for cell_type in cells:
        html_parts.append(f"<h2>{html.escape(cell_type)}</h2>")

        gene_sub = (
            genes[genes["cell_type"] == cell_type].copy()
            if not genes.empty
            else pd.DataFrame()
        )
        pathway_sub = (
            pathways[pathways["cell_type"] == cell_type].copy()
            if not pathways.empty
            else pd.DataFrame()
        )

        if not gene_sub.empty:
            for title, category in [
                ("Concordant reference DE hits", "REFERENCE_HIT_CONCORDANT"),
                ("Meaningfully discordant genes", "MEANINGFUL_DISCORDANCE"),
                ("Supported subset DE hits", "SUPPORTED_SUBSET_HIT"),
                ("Potential subset-specific DE hits", "POTENTIAL_SUBSET_SPECIFIC"),
            ]:
                display = gene_sub[
                    gene_sub["report_category"] == category
                ].sort_values("report_priority", ascending=False)
                html_parts.append(f"<h3>{title}</h3>")
                html_parts.append(
                    html_table(
                        display,
                        gene_columns,
                        args.top_genes,
                    )
                )

        if not pathway_sub.empty:
            for title, category in [
                ("Concordant reference pathways", "REFERENCE_HIT_CONCORDANT"),
                ("Meaningfully discordant pathways", "MEANINGFUL_DISCORDANCE"),
                ("Supported subset pathway hits", "SUPPORTED_SUBSET_HIT"),
                ("Potential subset-specific pathways", "POTENTIAL_SUBSET_SPECIFIC"),
            ]:
                display = pathway_sub[
                    pathway_sub["report_category"] == category
                ].sort_values("report_priority", ascending=False)
                html_parts.append(f"<h3>{title}</h3>")
                html_parts.append(
                    html_table(
                        display,
                        pathway_columns,
                        args.top_pathways,
                    )
                )

        if not recurrent_leading_edge.empty:

            leading_edge_sub = (
                recurrent_leading_edge[
                    (
                        recurrent_leading_edge[
                            "cell_type"
                        ]
                        == cell_type
                    )
                    &
                    recurrent_leading_edge[
                        "reportable"
                    ]
                ]
                .copy()
            )

            if not leading_edge_sub.empty:

                html_parts.append(
                    "<h3>"
                    "Recurrent leading-edge genes"
                    "</h3>"
                )

                html_parts.append(
                    "<p>"
                    "These genes recur in the leading edges of "
                    "concordant, FDR-significant pathways across "
                    "multiple analyses. Recurrence can identify "
                    "genes contributing repeatedly to coordinated "
                    "pathway-level responses even when individual "
                    "gene-level DE does not reach FDR significance."
                    "</p>"
                )

                leading_edge_columns = [
                    "gene",
                    "direction",
                    "n_analyses",
                    "n_pathways",
                    "n_gene_set_collections",
                    "n_occurrences",
                    "mean_abs_NES",
                    "min_FDR",
                    "analyses",
                    "gene_set_collections",
                ]

                html_parts.append(
                    html_table(
                        leading_edge_sub,
                        leading_edge_columns,
                        args.top_leading_edge_genes,
                    )
                )

        if args.make_plots:
            images = []

            # --------------------------------------------------
            # Gene log2FC heatmap
            # --------------------------------------------------

            gene_heatmap = (
                gene_heatmap_dir
                / f"{safe_filename(cell_type)}.png"
            )

            if gene_heatmap.exists():
                images.append(
                    (
                        "Selected gene effect sizes",
                        gene_heatmap,
                    )
                )

            # --------------------------------------------------
            # Recurrent leading-edge gene heatmap
            # --------------------------------------------------

            positive_leading_edge_heatmap = (
                leading_edge_heatmap_dir
                / (
                    f"{safe_filename(cell_type)}__"
                    "positive.png"
                )
            )

            if positive_leading_edge_heatmap.exists():
                images.append(
                    (
                        "Positive recurrent leading-edge genes",
                        positive_leading_edge_heatmap,
                    )
                )

            negative_leading_edge_heatmap = (
                leading_edge_heatmap_dir
                / (
                    f"{safe_filename(cell_type)}__"
                    "negative.png"
                )
            )

            if negative_leading_edge_heatmap.exists():
                images.append(
                    (
                        "Negative recurrent leading-edge genes",
                        negative_leading_edge_heatmap,
                    )
                )

            # --------------------------------------------------
            # Selected pathway NES heatmaps
            # --------------------------------------------------

            for gene_set in [
                "Hallmark",
                "GO_BP",
            ]:

                image = (
                    nes_heatmap_dir
                    / (
                        f"{safe_filename(cell_type)}__"
                        f"{safe_filename(gene_set)}.png"
                    )
                )

                if image.exists():
                    images.append(
                        (
                            f"{gene_set} pathway enrichment",
                            image,
                        )
                    )

            # --------------------------------------------------
            # Embed plots in HTML
            # --------------------------------------------------

            if images:

                html_parts.append(
                    "<h3>Selected cross-analysis heatmaps</h3>"
                )

                html_parts.append(
                    "<div class='plot-grid'>"
                )

                for title, image in images:

                    rel = image.relative_to(
                        report_dir
                    )

                    html_parts.append(
                        "<div class='plot-card'>"
                    )

                    html_parts.append(
                        f"<h4>{html.escape(title)}</h4>"
                    )

                    html_parts.append(
                        f"<img "
                        f"src='{html.escape(str(rel))}' "
                        f"alt='{html.escape(title)}'>"
                    )

                    html_parts.append(
                        "</div>"
                    )

                html_parts.append(
                    "</div>"
                )

    html_parts.append("</body></html>")

    with open(
        report_dir / "interpretation_report.html",
        "w",
        encoding="utf-8",
    ) as handle:
        handle.write("\n".join(html_parts))

    print(f"Interpretation outputs written to: {report_dir}")


if __name__ == "__main__":
    main()
