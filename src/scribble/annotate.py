#!/usr/bin/env python
from pathlib import Path

def normalise_refine_label(label):

    label = str(label)

    if label.startswith("L1_"):
        return label.replace("L1_", "")

    if label.startswith("L2_"):
        return label.replace("L2_", "")

    if label.startswith("L3_"):
        return label.replace("L3_", "")

    return label


def run_annotate(args):
    import scanpy as sc
    import pandas as pd
    import numpy as np
    import random
    import re
    import matplotlib.pyplot as plt

    from pathlib import Path
    from scribble.import_data import setup_environment

    PROJECT_DIR = Path(args.project_dir)
    PLOT_DIR = PROJECT_DIR / "scribble/plots"
    TABLE_DIR = PROJECT_DIR / "scribble/tables"

    setup_environment(sc, np, random, PLOT_DIR)

    input_file = Path(args.input)
    annotation_file = Path(args.annotations)

    output_file = input_file.with_name(
        f"{input_file.stem}_annotated.h5ad"
    )

    print(f"Loading {input_file}")
    adata = sc.read(input_file)

    print(f"Loading annotations → {annotation_file}")

    annotations = pd.read_excel(
        annotation_file,
        sheet_name="annotations"
    )

    if "refine_cluster" not in annotations.columns:
        raise ValueError(
            "annotations worksheet must contain "
            "'refine_cluster'"
        )


    # --------------------------------------------------
    # Map annotations onto refine labels
    # --------------------------------------------------

    annotations["refine_cluster"] = (
        annotations["refine_cluster"]
        .astype(str)
        .apply(normalise_refine_label)
    )

    anno_index = annotations.set_index("refine_cluster")

    # normalise atlas labels
    adata.obs["refine_label"] = (
        adata.obs["refine_label"]
        .astype(str)
        .apply(normalise_refine_label)
    )

    for col in annotations.columns:

        if col == "refine_cluster":
            continue

        mapping = anno_index[col].to_dict()

        adata.obs[col] = adata.obs["refine_label"].map(mapping)

    print("\nAnnotation coverage")

    n_total = adata.n_obs

    for col in annotations.columns:

        if col == "refine_cluster":
            continue

        n_missing = adata.obs[col].isna().sum()

        print(
            f"{col}: "
            f"{n_total - n_missing:,}/{n_total:,} "
            f"annotated "
            f"({100*(n_total-n_missing)/n_total:.1f}%)"
        )

    mapped = set(annotations["refine_cluster"])

    observed = set(
        adata.obs["refine_label"]
        .unique()
    )

    unmatched = sorted(observed - mapped)

    if unmatched:

        print("\nUnmatched refine labels")

        for x in unmatched:
            print(x)

    # --------------------------------------------------
    # UMAP cell type major
    # --------------------------------------------------

    print("Generating cell_type_major UMAP")

    sc.pl.umap(
        adata,
        color="cell_type_major",
        legend_loc="on data",
        legend_fontsize=8,
        show=False
    )

    plt.savefig(
        PLOT_DIR / "UMAP_cell_type_major.png",
        dpi=300,
        bbox_inches="tight"
    )

    plt.close()

    # --------------------------------------------------
    # Marker visualisation (optional)
    # --------------------------------------------------

    if args.plot_markers is not None:

        if "cell_type_major" not in annotations.columns:
            raise ValueError(
                "cell_type_major column required for "
                "marker visualisation"
            )

        if args.plot_markers not in annotations.columns:
            raise ValueError(
                f"{args.plot_markers} not found in "
                "annotations worksheet"
            )

        def safe_name(x):
            return re.sub(r"[^A-Za-z0-9._-]+", "_", str(x))

        unique_types = (
            annotations["cell_type_major"]
            .dropna()
            .unique()
        )

        marker_dict = {}

        # ------------------------------------------
        # Build marker dictionary
        # ------------------------------------------

        for cell_type in sorted(unique_types):

            subset = annotations.loc[
                annotations["cell_type_major"] == cell_type
            ]

            marker_set = set()

            for markers in subset[args.plot_markers].dropna():

                genes = [
                    g.strip()
                    for g in str(markers).split(";")
                    if g.strip()
                ]

                marker_set.update(genes)

            markers = sorted(
                g for g in marker_set
                if g in adata.var_names
            )

            if len(markers) == 0:

                print(
                    f"No markers found for {cell_type}"
                )

                continue

            marker_dict[cell_type] = markers

        # ------------------------------------------
        # Single dotplot
        # ------------------------------------------

        if len(marker_dict) > 0:

            print("Generating marker dotplot")

            sc.pl.dotplot(
                adata,
                var_names=marker_dict,
                groupby="cell_type_major",
                show=False
            )

            plt.savefig(
                PLOT_DIR / "dotplot_cell_type_major.png",
                dpi=300,
                bbox_inches="tight"
            )

            plt.close()

        # ------------------------------------------
        # Marker UMAPs
        # ------------------------------------------

        for cell_type, markers in marker_dict.items():

            print(f"Processing {cell_type}")

            for gene in markers:

                sc.pl.umap(
                    adata,
                    color=gene,
                    cmap="Reds",
                    title=f"{cell_type} | {gene}",
                    show=False
                )

                plt.savefig(
                    PLOT_DIR /
                    (
                        f"UMAP_"
                        f"{safe_name(cell_type)}_"
                        f"{safe_name(gene)}.png"
                    ),
                    dpi=300,
                    bbox_inches="tight"
                )

                plt.close()

    # --------------------------------------------------
    # Save
    # --------------------------------------------------

    print(f"Saving → {output_file}")

    adata.write(output_file)
