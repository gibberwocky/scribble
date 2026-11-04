#!/usr/bin/env python3
"""
Collate DE and GSEA outputs across multiple analysis directories.

Expected input layout
---------------------
<DGE_DIR>/<analysis>/DE/<annotation>/*_all_genes.csv
<DGE_DIR>/<analysis>/GSEA/<annotation>/*_ranking.csv
<DGE_DIR>/<analysis>/GSEA/<annotation>/<cell_type>/<gene_set>/all_pathways.csv
<DGE_DIR>/<analysis>/<annotation>_DE_GSEA_summary.csv

Outputs
-------
<DGE_DIR>/comparison/<annotation>/
    analysis_summary_long.csv
    analysis_summary_wide.csv
    DE/all_DE_long.csv
    DE/pairwise_log2FC_concordance.csv
    DE/by_cell_type/<cell_type>.csv
    GSEA/all_GSEA_long.csv
    GSEA/pairwise_NES_concordance.csv
    GSEA/by_cell_type/<cell_type>/<gene_set>.csv

The script does not alter source analysis outputs.
"""

import argparse
import re
from itertools import combinations
from pathlib import Path
import numpy as np
import pandas as pd
import json

def safe_filename(value):
    value = str(value).strip()
    value = re.sub(r"[^\w.-]+", "_", value)
    value = re.sub(r"_+", "_", value)
    return value.strip("_.")


def parse_args():
    p = argparse.ArgumentParser()
    p.add_argument("--analysis_dirs", nargs="+", required=True,
        help="Full paths to analysis directories, e.g. .../DGE/exvivo_major .../DGE/organoid_major")
    p.add_argument("--analysis_names", nargs="+", required=True,
        help="Names used in outputs corresponding to analysis_dirs")
    p.add_argument("--annotation", required=True,
        help="Annotation level, e.g. cell_type_major or cell_type_minor")
    p.add_argument("--reference_analysis", required=True,
        help="Primary/reference analysis used for reference-relative fields")
    p.add_argument("--comparison_name", required=True,
        help="Name for this comparison, e.g. 'biological_subsets' or 'organoid_model_sensitivity'.")
    p.add_argument("--outdir", required=True, help="Output directory for results")
    return p.parse_args()


def require_columns(df, columns, path):
    missing = [c for c in columns if c not in df.columns]
    if missing:
        raise ValueError(f"Missing columns {missing} in {path}")


def find_de_files(analysis_dir, annotation):
    root = Path(analysis_dir) / "DE" / annotation
    if not root.exists():
        return []
    return sorted(root.glob("*_all_genes.csv"))


def cell_type_from_de_filename(path):
    suffix = "_all_genes.csv"
    if not path.name.endswith(suffix):
        raise ValueError(f"Unexpected DE filename: {path}")
    return path.name[:-len(suffix)]


def load_summary(analysis_dir, analysis_name, annotation):
    path = (Path(analysis_dir) / f"{annotation}_DE_GSEA_summary.csv")
    if not path.exists():
        print(f"WARNING: summary missing: {path}")
        return None
    df = pd.read_csv(path)
    require_columns(df, ["cell_type"], path)
    if "analysis_name" not in df.columns:
        df.insert(0, "analysis_name", analysis_name)
    else:
        df["analysis_name"] = analysis_name
    df["analysis_dir"] = analysis_name
    return df


def collect_summaries(analysis_dirs, analysis_names, annotation):
    frames = []

    for analysis_dir, analysis_name in zip(
        analysis_dirs,
        analysis_names,
    ):

        df = load_summary(
            analysis_dir,
            analysis_name,
            annotation,
        )

        if df is not None:
            frames.append(df)

    return (
        pd.concat(
            frames,
            ignore_index=True,
            sort=False,
        )
        if frames
        else pd.DataFrame()
    )


def make_summary_wide(summary_long, analyses):
    if summary_long.empty:
        return pd.DataFrame()

    id_cols = {
        "analysis_name", "analysis_dir", "method", "cell_type",
        "formula", "contrast_variable", "contrast_test", "contrast_reference"
    }
    value_cols = [c for c in summary_long.columns if c not in id_cols]
    parts = []

    cell_types = sorted(summary_long["cell_type"].dropna().astype(str).unique())
    for cell_type in cell_types:
        row = {"cell_type": cell_type}
        sub = summary_long[summary_long["cell_type"].astype(str) == cell_type]
        for analysis in analyses:
            a = sub[sub["analysis_dir"] == analysis]
            if a.empty:
                continue
            rec = a.iloc[0]
            for col in value_cols:
                row[f"{analysis}_{col}"] = rec.get(col, np.nan)
        parts.append(row)
    return pd.DataFrame(parts)


def collect_de(analysis_dirs, analysis_names, annotation):
    frames = []
    for analysis_dir, analysis in zip(analysis_dirs, analysis_names):
        for path in find_de_files(analysis_dir, annotation):
            df = pd.read_csv(path)
            require_columns(
                df,
                [
                    "gene",
                    "baseMean",
                    "log2FoldChange",
                    "lfcSE",
                    "stat",
                    "pval",
                    "FDR",
                ],
                path,
            )

            df = df[
                [
                    "gene",
                    "baseMean",
                    "log2FoldChange",
                    "lfcSE",
                    "stat",
                    "pval",
                    "FDR",
                ]
            ].copy()
            df.insert(0, "cell_type", cell_type_from_de_filename(path))
            df.insert(0, "analysis", analysis)
            frames.append(df)
    return pd.concat(frames, ignore_index=True, sort=False) if frames else pd.DataFrame()


def classify_direction(values):
    vals = pd.Series(values, dtype=float).dropna()
    if vals.empty:
        return "NOT_TESTED"
    signs = np.sign(vals[vals != 0])
    if signs.empty:
        return "ZERO_ONLY"
    if (signs > 0).all():
        return "CONSISTENT_POSITIVE"
    if (signs < 0).all():
        return "CONSISTENT_NEGATIVE"
    return "MIXED_DIRECTION"


def make_de_wide_for_cell(de_cell, analyses, reference_analysis):
    metrics = ["baseMean", "log2FoldChange", "lfcSE", "stat", "pval", "FDR"]
    merged = None

    for analysis in analyses:
        a = de_cell[de_cell["analysis"] == analysis][["gene"] + metrics].copy()
        if a["gene"].duplicated().any():
            dups = a.loc[a["gene"].duplicated(keep=False), "gene"].unique().tolist()
            raise ValueError(
                f"Duplicate gene identifiers in {analysis}: {dups[:10]}"
            )
        a = a.rename(columns={m: f"{analysis}_{m}" for m in metrics})
        merged = a if merged is None else merged.merge(a, on="gene", how="outer")

    if merged is None:
        return pd.DataFrame()

    tested_cols = []
    lfc_cols = []
    for analysis in analyses:
        stat = f"{analysis}_stat"
        lfc = f"{analysis}_log2FoldChange"
        fdr = f"{analysis}_FDR"
        tested = f"{analysis}_tested"
        sig = f"{analysis}_FDR_005"
        merged[tested] = merged[stat].notna() if stat in merged.columns else False
        merged[sig] = merged[fdr].lt(0.05) if fdr in merged.columns else False
        tested_cols.append(tested)
        if lfc in merged.columns:
            lfc_cols.append(lfc)

    merged["n_analyses_tested"] = merged[tested_cols].sum(axis=1)
    merged["n_FDR_005"] = merged[[f"{a}_FDR_005" for a in analyses]].sum(axis=1)
    for analysis in analyses:
        sig = f"{analysis}_FDR_005"
        if sig in merged.columns:
            merged[f"{analysis}_significant"] = merged[sig]
    merged["direction_class"] = merged[lfc_cols].apply(classify_direction, axis=1)
    merged["all_same_direction"] = merged["direction_class"].isin(
        ["CONSISTENT_POSITIVE", "CONSISTENT_NEGATIVE"]
    )

    for a, b in combinations(analyses, 2):

        shared_sig = (
            merged[f"{a}_FDR_005"]
            &
            merged[f"{b}_FDR_005"]
        )

        merged[f"{a}_{b}_shared_sig"] = shared_sig

        merged[f"{a}_{b}_shared_sig_same_direction"] = (
            shared_sig
            &
            (
                np.sign(
                    merged[f"{a}_log2FoldChange"]
                )
                ==
                np.sign(
                    merged[f"{b}_log2FoldChange"]
                )
            )
        )

    ref_lfc = f"{reference_analysis}_log2FoldChange"
    if ref_lfc in merged.columns:
        for analysis in analyses:
            if analysis == reference_analysis:
                continue
            col = f"{analysis}_log2FoldChange"
            if col in merged.columns:
                merged[f"{analysis}_minus_{reference_analysis}_log2FC"] = merged[col] - merged[ref_lfc]
                merged[f"{analysis}_same_direction_as_{reference_analysis}"] = (
                    merged[col].notna()
                    & merged[ref_lfc].notna()
                    & (np.sign(merged[col]) == np.sign(merged[ref_lfc]))
                )

    sort_cols = []
    if f"{reference_analysis}_FDR" in merged.columns:
        sort_cols.append(f"{reference_analysis}_FDR")
    if f"{reference_analysis}_pval" in merged.columns:
        sort_cols.append(f"{reference_analysis}_pval")
    if sort_cols:
        merged = merged.sort_values(sort_cols, na_position="last")
    return merged


def stat_concordance(
    de_long,
    analyses,
    reference_analysis
):
    rows = []

    if de_long.empty:
        return pd.DataFrame()

    for cell_type, sub in (
        de_long.groupby(
            "cell_type",
            sort=True,
        )
    ):

        wide = make_de_wide_for_cell(
            sub,
            analyses,
            reference_analysis,
        )

        for a, b in combinations(
            analyses,
            2,
        ):

            x = f"{a}_stat"
            y = f"{b}_stat"

            if (
                x not in wide.columns
                or
                y not in wide.columns
            ):
                continue

            stats = correlation_pair(
                wide,
                x,
                y,
            )

            rows.append(
                {
                    "cell_type": cell_type,
                    "analysis_a": a,
                    "analysis_b": b,
                    **stats,
                }
            )

    return pd.DataFrame(rows)

def correlation_pair(df, x, y):
    tmp = df[[x, y]].replace([np.inf, -np.inf], np.nan).dropna()
    out = {"n_common": len(tmp), "pearson": np.nan, "spearman": np.nan}
    if len(tmp) >= 3 and tmp[x].nunique() > 1 and tmp[y].nunique() > 1:
        out["pearson"] = tmp[x].corr(tmp[y], method="pearson")
        out["spearman"] = tmp[x].corr(tmp[y], method="spearman")
    return out


def de_concordance(de_long, analyses, reference_analysis):
    rows = []
    if de_long.empty:
        return pd.DataFrame()
    for cell_type, sub in de_long.groupby("cell_type", sort=True):
        wide = make_de_wide_for_cell(sub, analyses, reference_analysis)
        for a, b in combinations(analyses, 2):
            x, y = f"{a}_log2FoldChange", f"{b}_log2FoldChange"
            if x not in wide.columns or y not in wide.columns:
                continue
            stats = correlation_pair(wide, x, y)
            sig_a = (
                wide[f"{a}_FDR"] < 0.05
            )

            sig_b = (
                wide[f"{b}_FDR"] < 0.05
            )

            shared_sig = (
                sig_a
                & sig_b
            )

            same_direction = (
                shared_sig
                &
                (
                    np.sign(
                        wide[x]
                    )
                    ==
                    np.sign(
                        wide[y]
                    )
                )
            )

            stats["n_sig_shared"] = (
                shared_sig.sum()
            )

            stats["n_sig_same_direction"] = (
                same_direction.sum()
            )

            stats["fraction_sig_same_direction"] = (
                same_direction.sum()
                / shared_sig.sum()
                if shared_sig.sum() > 0
                else np.nan
            )
            rows.append({"cell_type": cell_type, "analysis_a": a, "analysis_b": b, **stats})
    return pd.DataFrame(rows)


def find_gsea_files(analysis_dir, annotation):
    root = (Path(analysis_dir) / "GSEA" / annotation)
    if not root.exists():
        return []
    return sorted(root.glob("*/*/all_pathways.csv"))


def parse_gsea_path(path, analysis_dir, annotation):
    root = (Path(analysis_dir) / "GSEA" / annotation)
    rel = path.relative_to(root)
    if len(rel.parts) != 3:
        raise ValueError(f"Unexpected GSEA path: {path}")
    cell_type, gene_set, _ = rel.parts
    return cell_type, gene_set


def collect_gsea(analysis_dirs, analysis_names, annotation):
    frames = []
    keep = ["Term", "ES", "NES", "NOM p-val", "FDR q-val", "FWER p-val", "Tag %", "Gene %", "Lead_genes"]
    for analysis_dir, analysis in zip(analysis_dirs, analysis_names):
        for path in find_gsea_files(analysis_dir, annotation):
            cell_type, gene_set = parse_gsea_path(path, analysis_dir, annotation)
            df = pd.read_csv(path)
            require_columns(df, ["Term", "NES", "FDR q-val"], path)
            cols = [c for c in keep if c in df.columns]
            df = df[cols].copy()
            df.insert(0, "gene_set", gene_set)
            df.insert(0, "cell_type", cell_type)
            df.insert(0, "analysis", analysis)
            frames.append(df)
    return pd.concat(frames, ignore_index=True, sort=False) if frames else pd.DataFrame()


def make_gsea_wide(group, analyses, reference_analysis):
    metrics = ["ES", "NES", "NOM p-val", "FDR q-val", "FWER p-val", "Tag %", "Gene %", "Lead_genes"]
    merged = None
    for analysis in analyses:
        cols = ["Term"] + [m for m in metrics if m in group.columns]
        a = group[group["analysis"] == analysis][cols].copy()
        if a["Term"].duplicated().any():
            dups = a.loc[a["Term"].duplicated(keep=False), "Term"].unique().tolist()
            raise ValueError(f"Duplicate GSEA terms in {analysis}: {dups[:10]}")
        a = a.rename(columns={m: f"{analysis}_{m}" for m in cols if m != "Term"})
        merged = a if merged is None else merged.merge(a, on="Term", how="outer")
    if merged is None:
        return pd.DataFrame()

    tested_cols, nes_cols, sig_cols = [], [], []
    for analysis in analyses:
        nes = f"{analysis}_NES"
        fdr = f"{analysis}_FDR q-val"
        tested = f"{analysis}_tested"
        sig = f"{analysis}_FDR_005"
        merged[tested] = merged[nes].notna() if nes in merged.columns else False
        merged[sig] = merged[fdr].lt(0.05) if fdr in merged.columns else False
        tested_cols.append(tested); sig_cols.append(sig)
        if nes in merged.columns:
            nes_cols.append(nes)

    merged["n_analyses_tested"] = merged[tested_cols].sum(axis=1)
    merged["n_FDR_005"] = merged[sig_cols].sum(axis=1)
    merged["direction_class"] = merged[nes_cols].apply(classify_direction, axis=1)
    merged["all_same_direction"] = merged["direction_class"].isin(
        ["CONSISTENT_POSITIVE", "CONSISTENT_NEGATIVE"]
    )

    ref_nes = f"{reference_analysis}_NES"
    if ref_nes in merged.columns:
        for analysis in analyses:
            if analysis == reference_analysis:
                continue
            col = f"{analysis}_NES"
            if col in merged.columns:
                merged[f"{analysis}_minus_{reference_analysis}_NES"] = merged[col] - merged[ref_nes]
                merged[f"{analysis}_same_direction_as_{reference_analysis}"] = (
                    merged[col].notna()
                    & merged[ref_nes].notna()
                    & (np.sign(merged[col]) == np.sign(merged[ref_nes]))
                )

    sort_cols = []
    if f"{reference_analysis}_FDR q-val" in merged.columns:
        sort_cols.append(f"{reference_analysis}_FDR q-val")
    if f"{reference_analysis}_NES" in merged.columns:
        sort_cols.append(f"{reference_analysis}_NES")
    if sort_cols:
        asc = [True] + ([False] if len(sort_cols) > 1 else [])
        merged = merged.sort_values(sort_cols, ascending=asc, na_position="last")
    return merged


def gsea_concordance(gsea_long, analyses, reference_analysis):
    rows = []
    if gsea_long.empty:
        return pd.DataFrame()
    for (cell_type, gene_set), sub in gsea_long.groupby(["cell_type", "gene_set"], sort=True):
        wide = make_gsea_wide(sub, analyses, reference_analysis)
        for a, b in combinations(analyses, 2):
            x, y = f"{a}_NES", f"{b}_NES"
            if x not in wide.columns or y not in wide.columns:
                continue
            stats = correlation_pair(wide, x, y)
            sig_a = (
                wide[f"{a}_FDR q-val"] < 0.05
            )

            sig_b = (
                wide[f"{b}_FDR q-val"] < 0.05
            )

            shared_sig = (
                sig_a
                & sig_b
            )

            same_direction = (
                shared_sig
                &
                (
                    np.sign(
                        wide[x]
                    )
                    ==
                    np.sign(
                        wide[y]
                    )
                )
            )

            stats["n_pathways_shared"] = (
                shared_sig.sum()
            )

            stats["n_pathways_same_direction"] = (
                same_direction.sum()
            )

            stats["fraction_pathways_same_direction"] = (
                same_direction.sum()
                / shared_sig.sum()
                if shared_sig.sum() > 0
                else np.nan
            )
            rows.append({"cell_type": cell_type, "gene_set": gene_set, "analysis_a": a, "analysis_b": b, **stats})
    return pd.DataFrame(rows)


def main():
    args = parse_args()
    analysis_dirs = [
        Path(x).resolve()
        for x in args.analysis_dirs
    ]

    analyses = args.analysis_names

    if len(analysis_dirs) != len(analyses):

        raise ValueError(
            "--analysis_dirs and "
            "--analysis_names must "
            "have the same length"
        )

    if args.reference_analysis not in analyses:
        raise ValueError("--reference_analysis must also be listed in --analyses")

    missing_analysis_dirs = [
        str(x)
        for x in analysis_dirs
        if not x.exists()
    ]
    if missing_analysis_dirs:
        raise FileNotFoundError(f"Analysis directories not found: {missing_analysis_dirs}")

    comparison_name = safe_filename(
        args.comparison_name
    )

    outdir = Path(args.outdir)
    de_out = outdir / "DE"
    de_by_cell = de_out / "by_cell_type"
    gsea_out = outdir / "GSEA"
    gsea_by_cell = gsea_out / "by_cell_type"
    for p in [outdir, de_out, de_by_cell, gsea_out, gsea_by_cell]:
        p.mkdir(parents=True, exist_ok=True)

    # Run manifest
    comparison_manifest = {
        "comparison_name": args.comparison_name,
        "annotation": args.annotation,
        "reference_analysis": args.reference_analysis,
        "analysis_names": args.analysis_names,
        "analysis_dirs": args.analysis_dirs,
        "outdir": str(args.outdir)
    }

    with open(
        outdir / "comparison_manifest.json",
        "w",
    ) as handle:
        json.dump(
            comparison_manifest,
            handle,
            indent=2,
        )

    # Analysis summary
    summary_long = collect_summaries(
        analysis_dirs,
        analyses,
        args.annotation,
    )
    if not summary_long.empty:
        summary_long.to_csv(outdir / "analysis_summary_long.csv", index=False)
        summary_wide = make_summary_wide(summary_long, analyses)
        summary_wide.to_csv(outdir / "analysis_summary_wide.csv", index=False)

    # DE
    de_long = collect_de(
        analysis_dirs,
        analyses,
        args.annotation,
    )
    if not de_long.empty:
        de_long.to_csv(de_out / "all_DE_long.csv", index=False)
        for cell_type, sub in de_long.groupby("cell_type", sort=True):
            wide = make_de_wide_for_cell(sub, analyses, args.reference_analysis)
            wide.to_csv(de_by_cell / f"{safe_filename(cell_type)}.csv", index=False)
        de_corr = de_concordance(de_long, analyses, args.reference_analysis)
        de_corr.to_csv(de_out / "pairwise_log2FC_concordance.csv", index=False)
        stat_corr = stat_concordance(de_long, analyses, args.reference_analysis)
        stat_corr.to_csv(de_out / "pairwise_stat_concordance.csv", index=False)
    else:
        print("WARNING: no DE all_genes files found")

    # GSEA
    gsea_long = collect_gsea(
        analysis_dirs,
        analyses,
        args.annotation,
    )
    if not gsea_long.empty:
        gsea_long.to_csv(gsea_out / "all_GSEA_long.csv", index=False)
        for (cell_type, gene_set), sub in gsea_long.groupby(["cell_type", "gene_set"], sort=True):
            target = gsea_by_cell / safe_filename(cell_type)
            target.mkdir(parents=True, exist_ok=True)
            wide = make_gsea_wide(sub, analyses, args.reference_analysis)
            wide.to_csv(target / f"{safe_filename(gene_set)}.csv", index=False)
        gsea_corr = gsea_concordance(gsea_long, analyses, args.reference_analysis)
        gsea_corr.to_csv(gsea_out / "pairwise_NES_concordance.csv", index=False)
    else:
        print("WARNING: no GSEA all_pathways files found")

    print(f"Comparison outputs written to: {outdir}")


if __name__ == "__main__":
    main()
