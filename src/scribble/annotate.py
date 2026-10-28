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
    from scribble.refine import restore_counts
    from scribble.utils import write_run_manifest

    PROJECT_DIR = Path(args.project_dir)
    PLOT_DIR = PROJECT_DIR / "scribble/plots"
    MANIFEST_DIR = PROJECT_DIR / "scribble/manifests"

    setup_environment(sc, np, random, PLOT_DIR)

    input_file = Path(args.input)
    annotation_file = Path(args.annotations)
    manifest_file = (MANIFEST_DIR / f"{input_file.stem}_annotate_manifest.json")
    cluster_field = args.cluster_field
    generated_plots = []
    plot_prefix = (
        f"{args.plot_prefix}_"
        if args.plot_prefix
        else ""
    )
    umap_file = PLOT_DIR / f"{plot_prefix}UMAP_{args.label}.png"

    output_file = input_file.with_name(
        f"{input_file.stem}_annotated.h5ad"
    )

    print(f"Loading {input_file}")
    adata = sc.read(input_file)
    print(adata.shape)

    if cluster_field not in adata.obs:
        raise ValueError(
            f"{cluster_field} not found in adata.obs"
        )

    adata.obs[cluster_field] = (
        adata.obs[cluster_field]
        .astype(str)
    )

    # --------------------------------------------------
    # Generate sequential cluster IDs
    # --------------------------------------------------

    unique_clusters = (
        pd.Series(
            adata.obs[cluster_field]
            .astype(str)
        )
        .drop_duplicates()
        .tolist()
    )

    cluster_map = {
        cluster: str(i + 1)
        for i, cluster in enumerate(
            unique_clusters
        )
    }

    adata.uns[
        f"{args.cluster_destination}_map"
    ] = cluster_map

    adata.obs[
        args.cluster_destination
    ] = (
        adata.obs[cluster_field]
        .map(cluster_map)
    )

    print(
        f"\nCreated "
        f"{args.cluster_destination} "
        f"field"
    )

    print(cluster_map)



    try:
        adata_plot = restore_counts(adata)

        # preserve annotations/embeddings
        adata_plot.obsm = adata.obsm.copy()
        adata_plot.obs = adata.obs.copy()
        adata_plot.uns = adata.uns.copy()

        adata = adata_plot
        adata.uns.pop("log1p", None)

        sc.pp.normalize_total(adata)
        sc.pp.log1p(adata)
        adata.raw = adata

        print(
            f"Restored full feature space: "
            f"{adata.n_vars:,} genes"
        )

    except Exception as e:

        print(
            f"Could not restore counts: {e}"
        )

    print(f"Loading annotations → {annotation_file}")

    annotations = pd.read_excel(
        annotation_file,
        sheet_name="annotations",
        keep_default_na=False
    )

    if cluster_field not in annotations.columns:
        raise ValueError(
            f"annotations worksheet must contain "
            f"'{cluster_field}'"
        )

    if (
        args.cluster_hex is not None
        and args.cluster_hex not in annotations.columns
    ):
        raise ValueError(
            f"{args.cluster_hex} not found in "
            "annotations worksheet"
        )


    # --------------------------------------------------
    # Map annotations onto clusters
    # --------------------------------------------------

    annotations[cluster_field] = (
        annotations[cluster_field]
        .astype(str)
    )

    anno_index = annotations.set_index(cluster_field)

    # normalise atlas labels
    for col in annotations.columns:
        if col == cluster_field:
            continue

        mapping = anno_index[col].to_dict()

        adata.obs[col] = (
            adata.obs[cluster_field]
            .map(mapping)
        )

        # Convert any columns with more than one type (e.g. str and int) to str
        types = adata.obs[col].dropna().map(type).unique()
        if len(types) > 1:
            print(
                f"WARNING: {col!r} contains mixed types: "
                f"{[t.__name__ for t in types]}. Converting to string."
            )
            adata.obs[col] = adata.obs[col].astype("string")

    print("\nAnnotation coverage")

    n_total = adata.n_obs

    for col in annotations.columns:

        if col == cluster_field:
            continue

        n_missing = adata.obs[col].isna().sum()

        print(
            f"{col}: "
            f"{n_total - n_missing:,}/{n_total:,} "
            f"annotated "
            f"({100*(n_total-n_missing)/n_total:.1f}%)"
        )

    mapped = set(annotations[cluster_field])

    observed = set(
        adata.obs[cluster_field]
        .unique()
    )

    unmatched = sorted(observed - mapped)

    if unmatched:

        print("\nUnmatched cluster labels")

        for x in unmatched:
            print(x)

    # --------------------------------------------------
    # UMAP cell type
    # --------------------------------------------------
    plot_adata = adata
    if args.plot_samples is not None:
        plot_adata = adata[
            adata.obs["sample"].isin(
                args.plot_samples
            )
        ].copy()
        print(f"Plotting subset contains {plot_adata.n_obs:,} cells")

    palette = None

    if args.cluster_hex is not None:

        colour_df = (
            annotations[
                [args.label, args.cluster_hex]
            ]
            .dropna()
            .drop_duplicates()
        )

        n_celltypes = (
            colour_df[args.label]
            .nunique()
        )

        n_hex = (
            colour_df[args.cluster_hex]
            .nunique()
        )

        if n_celltypes == n_hex:

            palette = dict(
                zip(
                    colour_df[args.label],
                    colour_df[args.cluster_hex]
                )
            )

            print(
                f"Using {n_hex} colours from "
                f"{args.cluster_hex}"
            )

        else:

            print(
                "WARNING: "
                f"{n_celltypes} unique cell types but "
                f"{n_hex} unique colours. "
                "Using default Scanpy palette."
            )

    print(f"Generating {args.label} UMAP")

    if args.combined_label:
        plot_adata.obs["Annotation"] = (
            "[" + plot_adata.obs[cluster_field].astype(str) + "] "
            + plot_adata.obs[args.label].astype(str)
        )
    else:
        plot_adata.obs["Annotation"] = plot_adata.obs[args.label].astype(str)

    if palette is not None:

        # Scanpy expects categorical values
        plot_adata.obs["Annotation"] = (
            plot_adata.obs["Annotation"]
            .astype("category")
        )

        categories = (
            plot_adata.obs["Annotation"]
            .cat.categories
        )

        if set(categories).issubset(
            set(palette.keys())
        ):

            plot_adata.uns[
                "Annotation_colors"
            ] = [
                palette[x]
                for x in categories
            ]

    sc.pl.umap(
        plot_adata,
        color="Annotation",
        legend_loc=args.legend_loc,
        legend_fontsize=args.legend_fontsize,
        legend_fontoutline=args.legend_fontoutline,
        add_outline=args.add_outline,
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
        umap_file,
        dpi=300,
        bbox_inches="tight"
    )

    plt.close()

    generated_plots.append(str(umap_file))

    # --------------------------------------------------
    # Marker visualisation (optional)
    # --------------------------------------------------

    if args.plot_markers is not None:

        if args.label not in annotations.columns:
            raise ValueError(
                f"{args.label} column required for "
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
            annotations[args.label]
            .dropna()
            .unique()
        )

        marker_dict = {}

        # ------------------------------------------
        # Build marker dictionary
        # ------------------------------------------

        for cell_type in sorted(unique_types):

            subset = annotations.loc[
                annotations[args.label] == cell_type
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
                if g in plot_adata.var_names
            )

            missing_markers = sorted(
                g for g in marker_set
                if g not in plot_adata.var_names
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
                plot_adata,
                var_names=dotplot_dict,
                groupby=args.label,
                show=False
            )

            dotplot_file = PLOT_DIR / f"{plot_prefix}dotplot_{args.label}_{args.plot_markers}.png"
            plt.savefig(
                dotplot_file,
                dpi=300,
                bbox_inches="tight"
            )

            plt.close()

            generated_plots.append(str(dotplot_file))

        # ------------------------------------------
        # Marker UMAPs
        # ------------------------------------------

        for cell_type, markers in marker_dict.items():

            print(f"Processing {cell_type}")

            for gene in markers:

                sc.pl.umap(
                    plot_adata,
                    color=gene,
                    cmap="Reds",
                    title=f"{cell_type} | {gene}",
                    show=False
                )

                marker_umap_file = (
                    PLOT_DIR /
                    (
                        f"{plot_prefix}UMAP_"
                        f"{safe_name(cell_type)}_"
                        f"{safe_name(gene)}.png"
                    )
                )

                plt.savefig(
                    marker_umap_file,
                    dpi=300,
                    bbox_inches="tight"
                )

                plt.close()

                generated_plots.append(str(marker_umap_file))

    # --------------------------------------------------
    # Save
    # --------------------------------------------------

    if args.skip_write:
        print("Skipping writing adata file")
        output_file = "Skipped"
    else:
        print(f"Saving → {output_file}")
        adata.write(output_file)

    outputs = {
        "output adata": str(output_file),
        "plots": generated_plots,
    }

    write_run_manifest(
        manifest_file=manifest_file,
        module="annotate",
        args=args,
        inputs={
            "input_file": str(input_file),
        },
        outputs=outputs,
    )
