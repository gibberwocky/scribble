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


def safe_filename(value):
    value = str(value).strip()
    value = re.sub(r"[^\w.-]+", "_", value)
    value = re.sub(r"_+", "_", value)
    return value.strip("_.")


def parse_args():
    p = argparse.ArgumentParser()
    p.add_argument("--dge_dir", required=True,
                   help="Directory containing analysis subdirectories, e.g. project/scribble/DGE")
    p.add_argument("--analyses", nargs="+", required=True,
                   help="Analysis directory names to compare")
    p.add_argument("--annotation", required=True,
                   help="Annotation level, e.g. cell_type_major or cell_type_minor")
    p.add_argument("--reference_analysis", required=True,
                   help="Primary/reference analysis used for reference-relative fields")
    p.add_argument("--output_dir", default=None,
                   help="Optional output directory; default: <dge_dir>/comparison/<annotation>")
    return p.parse_args()


def require_columns(df, columns, path):
    missing = [c for c in columns if c not in df.columns]
    if missing:
        raise ValueError(f"Missing columns {missing} in {path}")


def find_de_files(dge_dir, analysis, annotation):
    root = dge_dir / analysis / "DE" / annotation
    if not root.exists():
        return []
    return sorted(root.glob("*_all_genes.csv"))


def cell_type_from_de_filename(path):
    suffix = "_all_genes.csv"
    if not path.name.endswith(suffix):
        raise ValueError(f"Unexpected DE filename: {path}")
    return path.name[:-len(suffix)]


def load_summary(dge_dir, analysis, annotation):
    path = dge_dir / analysis / f"{annotation}_DE_GSEA_summary.csv"
    if not path.exists():
        print(f"WARNING: summary missing: {path}")
        return None
    df = pd.read_csv(path)
    require_columns(df, ["cell_type"], path)
    if "analysis_name" not in df.columns:
        df.insert(0, "analysis_name", analysis)
    else:
        df["analysis_name"] = analysis
    df["analysis_dir"] = analysis
    return df


def collect_summaries(dge_dir, analyses, annotation):
    frames = []
    for analysis in analyses:
        df = load_summary(dge_dir, analysis, annotation)
        if df is not None:
            frames.append(df)
    return pd.concat(frames, ignore_index=True, sort=False) if frames else pd.DataFrame()


def make_summary_wide(summary_long, analyses):
    if summary_long.empty:
        return pd.DataFrame()

    id_cols = {
        "analysis_name", "analysis_dir", "cell_type",
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


def collect_de(dge_dir, analyses, annotation):
    frames = []
    for analysis in analyses:
        for path in find_de_files(dge_dir, analysis, annotation):
            df = pd.read_csv(path)
            require_columns(
                df,
                ["gene", "baseMean", "log2FoldChange", "lfcSE", "Wald_stat", "pval", "FDR"],
                path,
            )
            df = df[["gene", "baseMean", "log2FoldChange", "lfcSE", "Wald_stat", "pval", "FDR"]].copy()
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
    metrics = ["baseMean", "log2FoldChange", "lfcSE", "Wald_stat", "pval", "FDR"]
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
        wald = f"{analysis}_Wald_stat"
        lfc = f"{analysis}_log2FoldChange"
        fdr = f"{analysis}_FDR"
        tested = f"{analysis}_tested"
        sig = f"{analysis}_FDR_005"
        merged[tested] = merged[wald].notna() if wald in merged.columns else False
        merged[sig] = merged[fdr].lt(0.05) if fdr in merged.columns else False
        tested_cols.append(tested)
        if lfc in merged.columns:
            lfc_cols.append(lfc)

    merged["n_analyses_tested"] = merged[tested_cols].sum(axis=1)
    merged["n_FDR_005"] = merged[[f"{a}_FDR_005" for a in analyses]].sum(axis=1)
    merged["direction_class"] = merged[lfc_cols].apply(classify_direction, axis=1)
    merged["all_same_direction"] = merged["direction_class"].isin(
        ["CONSISTENT_POSITIVE", "CONSISTENT_NEGATIVE"]
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


def correlation_pair(df, x, y):
    tmp = df[[x, y]].replace([np.inf, -np.inf], np.nan).dropna()
    out = {"n_common": len(tmp), "pearson": np.nan, "spearman": np.nan}
    if len(tmp) >= 3 and tmp[x].nunique() > 1 and tmp[y].nunique() > 1:
        out["pearson"] = tmp[x].corr(tmp[y], method="pearson")
        out["spearman"] = tmp[x].corr(tmp[y], method="spearman")
    return out


def de_concordance(de_long, analyses):
    rows = []
    if de_long.empty:
        return pd.DataFrame()
    for cell_type, sub in de_long.groupby("cell_type", sort=True):
        wide = make_de_wide_for_cell(sub, analyses, analyses[0])
        for a, b in combinations(analyses, 2):
            x, y = f"{a}_log2FoldChange", f"{b}_log2FoldChange"
            if x not in wide.columns or y not in wide.columns:
                continue
            stats = correlation_pair(wide, x, y)
            rows.append({"cell_type": cell_type, "analysis_a": a, "analysis_b": b, **stats})
    return pd.DataFrame(rows)


def find_gsea_files(dge_dir, analysis, annotation):
    root = dge_dir / analysis / "GSEA" / annotation
    if not root.exists():
        return []
    return sorted(root.glob("*/*/all_pathways.csv"))


def parse_gsea_path(path, dge_dir, analysis, annotation):
    root = dge_dir / analysis / "GSEA" / annotation
    rel = path.relative_to(root)
    if len(rel.parts) != 3:
        raise ValueError(f"Unexpected GSEA path: {path}")
    cell_type, gene_set, _ = rel.parts
    return cell_type, gene_set


def collect_gsea(dge_dir, analyses, annotation):
    frames = []
    keep = ["Term", "ES", "NES", "NOM p-val", "FDR q-val", "FWER p-val", "Tag %", "Gene %", "Lead_genes"]
    for analysis in analyses:
        for path in find_gsea_files(dge_dir, analysis, annotation):
            cell_type, gene_set = parse_gsea_path(path, dge_dir, analysis, annotation)
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


def gsea_concordance(gsea_long, analyses):
    rows = []
    if gsea_long.empty:
        return pd.DataFrame()
    for (cell_type, gene_set), sub in gsea_long.groupby(["cell_type", "gene_set"], sort=True):
        wide = make_gsea_wide(sub, analyses, analyses[0])
        for a, b in combinations(analyses, 2):
            x, y = f"{a}_NES", f"{b}_NES"
            if x not in wide.columns or y not in wide.columns:
                continue
            stats = correlation_pair(wide, x, y)
            rows.append({"cell_type": cell_type, "gene_set": gene_set, "analysis_a": a, "analysis_b": b, **stats})
    return pd.DataFrame(rows)


def main():
    args = parse_args()
    dge_dir = Path(args.dge_dir).resolve()
    analyses = args.analyses

    if args.reference_analysis not in analyses:
        raise ValueError("--reference_analysis must also be listed in --analyses")

    missing_analysis_dirs = [a for a in analyses if not (dge_dir / a).exists()]
    if missing_analysis_dirs:
        raise FileNotFoundError(f"Analysis directories not found: {missing_analysis_dirs}")

    outdir = Path(args.output_dir).resolve() if args.output_dir else dge_dir / "comparison" / args.annotation
    de_out = outdir / "DE"
    de_by_cell = de_out / "by_cell_type"
    gsea_out = outdir / "GSEA"
    gsea_by_cell = gsea_out / "by_cell_type"
    for p in [outdir, de_out, de_by_cell, gsea_out, gsea_by_cell]:
        p.mkdir(parents=True, exist_ok=True)

    # Analysis summary
    summary_long = collect_summaries(dge_dir, analyses, args.annotation)
    if not summary_long.empty:
        summary_long.to_csv(outdir / "analysis_summary_long.csv", index=False)
        summary_wide = make_summary_wide(summary_long, analyses)
        summary_wide.to_csv(outdir / "analysis_summary_wide.csv", index=False)

    # DE
    de_long = collect_de(dge_dir, analyses, args.annotation)
    if not de_long.empty:
        de_long.to_csv(de_out / "all_DE_long.csv", index=False)
        for cell_type, sub in de_long.groupby("cell_type", sort=True):
            wide = make_de_wide_for_cell(sub, analyses, args.reference_analysis)
            wide.to_csv(de_by_cell / f"{safe_filename(cell_type)}.csv", index=False)
        de_corr = de_concordance(de_long, analyses)
        de_corr.to_csv(de_out / "pairwise_log2FC_concordance.csv", index=False)
    else:
        print("WARNING: no DE all_genes files found")

    # GSEA
    gsea_long = collect_gsea(dge_dir, analyses, args.annotation)
    if not gsea_long.empty:
        gsea_long.to_csv(gsea_out / "all_GSEA_long.csv", index=False)
        for (cell_type, gene_set), sub in gsea_long.groupby(["cell_type", "gene_set"], sort=True):
            target = gsea_by_cell / safe_filename(cell_type)
            target.mkdir(parents=True, exist_ok=True)
            wide = make_gsea_wide(sub, analyses, args.reference_analysis)
            wide.to_csv(target / f"{safe_filename(gene_set)}.csv", index=False)
        gsea_corr = gsea_concordance(gsea_long, analyses)
        gsea_corr.to_csv(gsea_out / "pairwise_NES_concordance.csv", index=False)
    else:
        print("WARNING: no GSEA all_pathways files found")

    print(f"Comparison outputs written to: {outdir}")


if __name__ == "__main__":
    main()
