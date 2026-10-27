#!/usr/bin/env python

from pathlib import Path

def update_metadata(
    adata,
    metadata_file,
    fields,
):
    import pandas as pd

    meta = pd.read_excel(
        metadata_file,
        sheet_name="meta",
        keep_default_na=False,
    )

    if "sample" not in meta.columns:
        raise ValueError(
            "Metadata must contain a "
            "'sample' column"
        )

    adata.obs["sample"] = (
        adata.obs["sample"]
        .astype(str)
    )

    meta["sample"] = (
        meta["sample"]
        .astype(str)
    )

    meta = meta.set_index("sample")

    updated_fields = []

    for field in fields:

        if field not in meta.columns:

            print(
                f"WARNING: {field} not present "
                f"in metadata workbook"
            )

            continue

        mapping = meta[field].to_dict()

        adata.obs[field] = (
            adata.obs["sample"]
            .map(mapping)
        )

        updated_fields.append(field)

    # ---------------------------------------
    # Fix mixed datatypes
    # ---------------------------------------

    for field in updated_fields:

        types = (
            adata.obs[field]
            .dropna()
            .map(type)
            .unique()
        )

        if len(types) > 1:

            print(
                f"{field}: mixed types "
                f"{[x.__name__ for x in types]}"
                f" -> converting to string"
            )

            adata.obs[field] = (
                adata.obs[field]
                .astype("string")
            )

    return adata


def run_metadata(args):
    import scanpy as sc
    import pandas as pd
    import numpy as np
    from datetime import datetime
    from scribble.import_data import (
        setup_environment
    )

    from scribble.utils import (
        serializable_args,
        write_run_manifest,
    )

    PROJECT_DIR = Path(args.project_dir)
    MANIFEST_DIR = (PROJECT_DIR / "scribble/manifests")

    input_file = Path(args.input)
    metadata_file = Path(args.metadata_file)

    tag = datetime.now().strftime("%Y%m%d_%H%M%S")
    manifest_file = (MANIFEST_DIR / f"{input_file.stem}_metadata_{tag}_manifest.json")

    print(f"Loading AnnData → {input_file}")
    adata = sc.read(input_file)

    adata = update_metadata(
        adata,
        metadata_file,
        args.fields,
    )

    print(f"Saving -> {input_file}")
    adata.write(input_file)

    write_run_manifest(
        manifest_file=manifest_file,
        module="metadata",
        args=args,
        inputs={
            "input_file": str(input_file),
            "metadata_file": str(metadata_file),
        },
        outputs={
            "adata": str(input_file),
        }
    )
