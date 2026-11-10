#!/usr/bin/env python

from pathlib import Path


def run_cellxgene(args):
    import scanpy as sc
    import anndata as ad

    input_file = Path(args.input)
    output_file = Path(args.output)

    if Path(output_file).exists():
        continue
    else:
        print(f"{input_file} not found")

    if input_file == output_file:
        raise ValueError("Output file must be different to input file")

    if Path(output_file).suffix.lower() != ".h5ad":
        raise ValueError("Output file must have the .h5ad extentions")

    # Import data
    adata = sc.read(input_file)

    ad.settings.allow_write_nullable_strings = False
    adata_cxg = adata.copy()
    adata_cxg.write_h5ad(output_file, compression="gzip")

    print(f"cellxgene compatible output written to {output_file}")
