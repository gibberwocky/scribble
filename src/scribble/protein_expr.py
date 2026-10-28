#!/usr/bin/env python
from pathlib import Path

def run_protein_expr(args):

    from pathlib import Path
    import scanpy as sc
    import pandas as pd
    import numpy as np
    from scribble.refine import restore_counts
    from scribble.utils import write_run_manifest


    PROJECT_DIR = Path(args.project_dir)
    MANIFEST_DIR = PROJECT_DIR / "scribble/manifests"
    OUT_DIR = PROJECT_DIR / "scribble/protein"

    OUT_DIR.mkdir(parents=True, exist_ok=True)

    mean_expr_file = (OUT_DIR / f"{args.output_prefix}_mean_expression.tsv")
    frac_expr_file = (OUT_DIR / f"{args.output_prefix}_fraction_expressing.tsv")
    score_expr_file = (OUT_DIR / f"{args.output_prefix}_expression_score.tsv")
    top_sources_file = (OUT_DIR / f"{args.output_prefix}_top_sources.tsv")
    protein_counts_file = (OUT_DIR / f"{args.output_prefix}_counts_by_celltype.tsv")
    manifest_file = (MANIFEST_DIR / f"{args.output_prefix}_proteins_manifest.json")

    # Import adata
    adata = sc.read_h5ad(args.input_file)

    # Check samples requested are present
    if "sample" not in adata.obs.columns:
        raise ValueError(
            "Required column 'sample' "
            "is not present in adata.obs."
        )

    # Restore raw counts
    adata_counts = restore_counts(adata)

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

    # Subset samples
    adata_counts = adata_counts[
        adata_counts.obs["sample"].isin(
            args.samples
        )
    ].copy()

    # Normalise + log transform
    sc.pp.normalize_total(adata_counts)
    sc.pp.log1p(adata_counts)

    # Load proteins
    proteins = pd.read_csv(args.proteins)
    genes = []
    for value in proteins["genes"]:
        value = str(value).strip()
        if not value:
            continue
        # Split semicolon-delimited entries
        for gene in value.split(";"):
            gene = gene.strip()
            if gene:
                genes.append(gene)

    # retain genes present in AnnData
    genes = sorted(set(genes))
    genes = [
        g
        for g in genes
        if g in adata_counts.var_names
    ]
    print(f"{len(genes)} proteins found in AnnData")

    # Expression matrices
    celltypes = sorted(adata_counts.obs[args.annotation].dropna().unique())
    mean_expr = pd.DataFrame(index=genes, columns=celltypes, dtype=float)
    frac_expr = pd.DataFrame(index=genes, columns=celltypes, dtype=float)
    score_expr = pd.DataFrame(index=genes, columns=celltypes, dtype=float)

    for cell_type in celltypes:
        subset = adata_counts[adata_counts.obs[args.annotation] == cell_type]
        X = subset[:, genes].X
        mean_values = np.asarray(X.mean(axis=0)).ravel()
        frac_values = np.asarray((X > 0).mean(axis=0)).ravel()
        score_values = (mean_values * frac_values)
        mean_expr[cell_type] = mean_values
        frac_expr[cell_type] = frac_values
        score_expr[cell_type] = score_values

    mean_expr.to_csv(mean_expr_file, sep="\t")
    frac_expr.to_csv(frac_expr_file, sep="\t")
    score_expr.to_csv(score_expr_file, sep="\t")

    # --------------------------------------
    # Top 3 likely source cell types
    # --------------------------------------
    top_rows = []
    for gene in score_expr.index:
        valid = (
            (mean_expr.loc[gene] > args.expression_threshold)
            &
            (frac_expr.loc[gene] > args.fraction_threshold)
        )
        ranked = (score_expr.loc[gene][valid].sort_values(ascending=False))
        top = ranked.head(3)
        row = {
            "gene": gene,
            "n_sources": len(ranked),
        }
        for i, (ct, score) in enumerate(
            top.items(),
            start=1,
        ):
            row[f"source_{i}"] = ct
            row[f"score_{i}"] = score
        top_rows.append(row)

    top_df = pd.DataFrame(top_rows)
    top_df.to_csv(top_sources_file, sep="\t", index=False)

    # --------------------------------------
    # Protein count per cell type
    # --------------------------------------
    counts_rows = []
    for cell_type in celltypes:
        n_detected = (
            (mean_expr[cell_type] > args.expression_threshold)
            &
            (frac_expr[cell_type] > args.fraction_threshold)
        ).sum()
        counts_rows.append({
            "cell_type": cell_type,
            "n_proteins_detected": int(n_detected),
        })

    counts_df = pd.DataFrame(counts_rows).sort_values("n_proteins_detected", ascending=False)
    counts_df.to_csv(protein_counts_file, sep="\t", index=False)

    # Summary
    print("\nTop source examples")
    print(top_df.head(10))
    print("\nProteins detected per cell type")
    print(counts_df)

    # Write manifest
    write_run_manifest(
        manifest_file=manifest_file,
        module="proteins",
        args=args,
        inputs={
            "input_file": str(args.input_file),
        },
        outputs={
            "mean_expr tsv": str(mean_expr_file),
            "frac_expr tsv": str(frac_expr_file),
            "score_expr tsv": str(score_expr_file),
            "top_sources tsv": str(top_sources_file),
            "protein_counts tsv": str(protein_counts_file),
        }
    )
