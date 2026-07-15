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
        legend_loc="right margin",
        frameon=False,
        show=False
    )

    ax = plt.gca()


    # Get Scanpy-generated legend
    leg = ax.get_legend()

    if leg is not None:

        handles = leg.legend_handles
        labels = [t.get_text() for t in leg.get_texts()]

        leg.remove()

        ax.legend(
            handles,
            labels,
            loc="center left",
            bbox_to_anchor=(1.01, 0.5),
            ncol=1,
            fontsize=4,
            frameon=False,
            markerscale=0.25,
            handletextpad=0.3,
            labelspacing=0.2
        )

    for txt in ax.texts:
        txt.set_fontweight("normal")
        txt.set_alpha(1.0)

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
                    for g in re.split(r"[;,]", str(markers))
                    if g.strip()
                ]

                marker_set.update(genes)

            markers = sorted(
                g for g in marker_set
                if g in adata.var_names
            )

            missing_markers = sorted(
                g for g in marker_set
                if g not in adata.var_names
            )

            if len(markers) == 0:
                print(
                    f"No markers found for {cell_type}"
                )
                continue

            print(
                f"{cell_type}: "
                f"{len(markers)} markers found"
            )
            print(markers[:10])

            if len(missing_markers) > 0:
                print(
                    f"{cell_type}: "
                    f"{len(missing_markers)} markers not found "
                    f"in adata"
                )
                print(missing_markers[:10])

            marker_dict[cell_type] = markers

        # ------------------------------------------
        # Single dotplot
        # ------------------------------------------

        if len(marker_dict) > 0:

            print("Generating marker dotplot")

            # ------------------------------------------
            # Remove duplicated genes across cell types
            # for dotplot only
            # ------------------------------------------

            used_genes = set()
            dotplot_dict = {}

            for cell_type, markers in marker_dict.items():

                unique_markers = []

                for gene in markers:

                    if gene in used_genes:
                        continue

                    unique_markers.append(gene)
                    used_genes.add(gene)

                if len(unique_markers) > 0:
                    dotplot_dict[cell_type] = unique_markers

            sc.pl.dotplot(
                adata,
                var_names=dotplot_dict,
                groupby="cell_type_major",
                show=False
            )

            plt.savefig(
                PLOT_DIR / f"dotplot_cell_type_major_{args.plot_markers}.png",
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
