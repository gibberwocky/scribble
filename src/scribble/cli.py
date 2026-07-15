#!/usr/bin/env python

__author__ = "David Wragg"
__license__ = "GNU GPL v3.0"

import argparse
#from pandas.core.arrays.base import InterpolateOptions

def main():
    parser = argparse.ArgumentParser(prog="scribble",
        description=f"""
\033[34m
╭─╮╭─╴╭─╮╷╭╮ ╭╮ ╷  ╭─╴
╰─╮│  ├┬╯│├┴╮├┴╮│  ├╴
╰─╯╰─╴╵╰╴╵╰─╯╰─╯╰─╴╰─╴
\033[32m
A single-cell RNA-seq workflow CLI built on Scanpy.

Expected project structure:
- cellranger outputs:
    <project>/cellranger/<sample>/outs/filtered_feature_bc_matrix
- velocyto outputs:
    <project>/velocyto/<sample>/*.loom
- metadata Excel file:
    must include sheet 'meta' with column 'sample'

\033[0m
GitHub: https://github.com/gibberwocky/scribble
""",
    formatter_class=argparse.RawTextHelpFormatter,
)
    subparsers = parser.add_subparsers(
        dest="command",
        required=True,
        title="Commands",
        description="Available processing steps"
    )

    # ---------- IMPORT COMMAND ----------
    import_parser = subparsers.add_parser("import")
    import_parser.add_argument("--project_dir", required=True,
        help='Root project directory')
    import_parser.add_argument("--cellranger_dir", required=True,
        help='Directory containing cellranger/ sample subdirectories')
    import_parser.add_argument("--velocyto_dir", required=True,
        help='Directory containing velocyto/ sample subdirectories')
    import_parser.add_argument("--samples", nargs="+",
        help="Samples to process")
    import_parser.add_argument("--samples-file",
        help="One sample per line")
    import_parser.add_argument("--metadata_file", required=True,
        help="Excel metadata file with sheet 'meta' containing a 'sample' column matching input folder names")
    import_parser.add_argument("--inflection_lower", type=int, default=100,
        help="Lower bound for total cell UMIs when calculation inflection point")

    # -------------- MT QC ---------------
    mt_parser = subparsers.add_parser("mt")
    mt_parser.add_argument("--project_dir", required=True,
        help="Project directory containing scribble output folder")
    mt_parser.add_argument("--input", required=True,
        help="Input AnnData (.h5ad) file (generated with scribble)")
    mt_parser.add_argument("--nmads", type=int, default=3,
        help="Number of median absolute deviations used to threshold mitochondrial content")

    # ------------ Doublet QC ------------
    dbl_parser = subparsers.add_parser("dbl")
    dbl_parser.add_argument("--project_dir", required=True,
        help="Project directory containing scribble output folder")
    dbl_parser.add_argument("--input", required=True,
        help="Input AnnData (.h5ad) file (generated with scribble)")
    dbl_parser.add_argument("--expected", type=float, default=0.07,
        help=("Expected doublet rate (used for Scrublet prior and quantile-based thresholding). "
              "Actual doublet fraction may differ due to hybrid thresholding."))
    dbl_parser.add_argument("--mode", choices=["hybrid", "scrublet", "quantile"], default="hybrid",
        help="Doublet calling strategy")
    dbl_parser.add_argument("--normalize", action="store_true",
        help="Apply normalization before Scrublet")
    dbl_parser.add_argument("--min_cells", type=int, default=200,
        help="Warn if sample has fewer cells than this for Scrublet")

    # -------------- PCA QC --------------
    pca_parser = subparsers.add_parser("pca")
    pca_parser.add_argument("--project_dir", required=True,
        help="Project directory containing scribble output folder")
    pca_parser.add_argument("--input", required=True,
        help="Input AnnData (.h5ad) file (generated with scribble)")
    pca_parser.add_argument("--mingenes", type=int, default=200,
        help="Minimum number of detected genes required for a cell")
    pca_parser.add_argument("--maxgenes", type=int, default=6000,
        help="Maximum number of detected genes allowed for a cell")
    pca_parser.add_argument("--hvgs", type=int, default=3000,
        help="Number of highly variable genes used for PCA")
    pca_parser.add_argument("--vmax", type=float, default=0.99,
        help="Upper limit for PCA variance ratio plot scaling")

    # -------------- Filter --------------
    filter_parser = subparsers.add_parser("filter")
    filter_parser.add_argument("--project_dir", required=True,
        help="Project directory containing scribble output folder")
    filter_parser.add_argument("--input", required=True,
        help="Input AnnData (.h5ad) file (generated with scribble)")
    filter_parser.add_argument("--mingenes", type=int, default=200,
        help="Minimum number of genes with non-zero counts in a given cell")
    filter_parser.add_argument("--maxgenes", type=int, default=6000,
        help="Maximum number of genes with non-zero counts in a given cell")
    filter_parser.add_argument("--filter_xlsx",
        help="Optional Excel file with per-sample filtering thresholds (sheet 'filters', must include 'sample' column)")

    # ---------- Preintegration ----------
    preintegration_parser = subparsers.add_parser("preintegration")
    preintegration_parser.add_argument("--project_dir", required=True,
        help="Project directory containing scribble output folder")
    preintegration_parser.add_argument("--input", required=True,
        help="Input AnnData (.h5ad) file (generated with scribble)")
    preintegration_parser.add_argument("--hvgs", type=int, default=3000,
        help="Number of highly variable genes")
    preintegration_parser.add_argument("--npcs", type=int, default=30,
        help="Number of principal components to compute")
    preintegration_parser.add_argument("--neighbors", type=int, default=15,
        help="Number of nearest neighbours used to build the kNN graph")
    preintegration_parser.add_argument("--batch", type=str, default="sample",
        help="Column in adata.obs used to define batches for HVG selection")
    preintegration_parser.add_argument("--vars", nargs="+", default=["sample"],
        help="Observation columns to colour PCA/UMAP plots")
    preintegration_parser.add_argument("--no-scale", action="store_true",
        help="Disable data scaling prior to PCA (enabled by default)")
    preintegration_parser.add_argument("--regress", nargs="+", default=None,
        help="Variables to regress out (e.g. total_counts pct_counts_mt)")
    preintegration_parser.add_argument("--min_cells_per_gene", "--min_cells", dest="min_cells_per_gene", type=int, default=3,
        help="Minimum number of cells a gene must be expressed in")

    # -------------- Harmony -------------
    harmony_parser = subparsers.add_parser("harmony")
    harmony_parser.add_argument("--project_dir", required=True,
        help="Project directory containing scribble output folder")
    harmony_parser.add_argument("--input", required=True,
        help="Input AnnData (.h5ad) file (generated with scribble)")
    harmony_parser.add_argument("--npcs", type=int, default=50,
        help="Number of principal components used for Harmony integration")
    harmony_parser.add_argument("--neighbors", type=int, default=15,
        help="Number of nearest neighbours used to build the integrated graph")
    harmony_parser.add_argument("--batch", type=str, default="sample",
        help="Batch key in adata.obs (e.g. 'sample')")
    harmony_parser.add_argument("--theta", type=int, default=2,
        help="Harmony diversity penalty (higher values increase batch mixing but risk overcorrection)")
    harmony_parser.add_argument("--vars", nargs="+", default=["sample"],
        help="Observation columns to colour UMAP plots")

    # -------------- Cluster -------------
    cluster_parser = subparsers.add_parser("cluster",
        help="Cluster cells and compute stability + marker genes")
    cluster_parser.add_argument("--project_dir", required=True,
        help="Project directory containing scribble output folder")
    cluster_parser.add_argument("--input", required=True,
        help="Input AnnData (.h5ad) file (generated with scribble)")
    cluster_parser.add_argument("--output_prefix", required=True,
        help="Prefix for output files (defaults to args.input.stem)")
    cluster_parser.add_argument("--embedding", default="X_pca",
        help="Embedding to use for clustering ('X_pca' or 'X_pca_harmony')")
    cluster_parser.add_argument("--resolution", type=float, default=1.0,
        help="Leiden resolution parameter controlling cluster granularity")
    cluster_parser.add_argument("--neighbors", type=int, default=15,
        help="Number of nearest neighbours used to build the clustering graph")
    cluster_parser.add_argument("--vars", nargs="+", default=["sample"],
        help="Observation columns to colour UMAP plots alongside cluster labels")
    cluster_parser.add_argument("--n_repeats", type=int, default=10,
        help="Number of repeated clustering runs used to estimate cluster stability")
    cluster_parser.add_argument("--auto_resolution", action="store_true",
        help="Automatically select optimal Leiden resolution based on silhouette score")
    cluster_parser.add_argument("--res_min", type=float, default=0.2,
        help="Lower bound for automatic resolution scan")
    cluster_parser.add_argument("--res_max", type=float, default=2.0,
        help="Upper bound for automatic resolution scan")
    cluster_parser.add_argument("--res_steps", type=int, default=10,
        help="Number of resolution values tested")
    cluster_parser.add_argument("--fine_width", type=float, default=0.3,
        help="Width of fine search window around best coarse resolution")
    cluster_parser.add_argument("--nmarkers", type=int, default=100,
        help="Number of top marker genes exported per cluster")

    # -------------- Evaluate ------------
    evaluate_parser = subparsers.add_parser("evaluate",
        help="Evaluate clusters and generate refinement decisions (subset / merge / keep)")
    evaluate_parser.add_argument("--project_dir", required=True,
        help="Project directory containing scribble output folder")
    evaluate_parser.add_argument("--input", required=True, nargs="+",
        help="One or more cluster_summary.tsv files generated by scribble cluster")
    evaluate_parser.add_argument("--min_cells", type=int, default=200,
        help="Minimum cluster size below which clusters are ignored for refinement")
    evaluate_parser.add_argument("--large_cells", type=int, default=800,
        help="Minimum size threshold for considering a cluster as large and potentially heterogeneous")
    evaluate_parser.add_argument("--low_stability", type=float, default=0.75,
        help="Threshold below which clusters are considered unstable (potential continua)")
    evaluate_parser.add_argument("--high_stability", type=float, default=0.95,
        help="Threshold above which clusters are considered highly stable")
    evaluate_parser.add_argument("--low_entropy", type=float, default=0.5,
        help="Threshold below which clusters are considered poorly mixed across samples")
    evaluate_parser.add_argument("--merge_size_ratio", type=float, default=2.5,
        help="Maximum size ratio allowed when merging clusters")
    evaluate_parser.add_argument("--merge_stability_tol", type=float, default=0.1,
        help="Maximum difference in stability allowed when merging clusters")
    evaluate_parser.add_argument("--merge_entropy_tol", type=float, default=0.2,
        help="Maximum difference in entropy allowed when merging clusters")

    # -------------- Evaluate ------------
    refine_parser = subparsers.add_parser("refine",
        description="""
        Refine selected clusters by:
        - subsetting cells using cluster decision file
        - re-running HVG selection, PCA, and Harmony integration
        - optimising clustering resolution
        - computing refined clusters and stability
        - identifying marker genes

        Outputs:
        - updated AnnData with hierarchical labels (leiden_L2)
        - marker gene Excel files
        """)
    refine_parser.add_argument("--project_dir", required=True,
        help="Project directory containing scribble/ subfolders")
    refine_parser.add_argument("--input", required=True,
        help="Input AnnData (.h5ad) file (generated with scribble)")
    refine_parser.add_argument("--decisions", required=True,
        help="Cluster decisions TSV file")
    refine_parser.add_argument("--batch", type=str, default="sample",
        help="Column in adata.obs used as batch key for Harmony integration")
    refine_parser.add_argument("--hvgs", type=int, default=3000,
        help="Number of highly variable genes selected during reprocessing")
    refine_parser.add_argument("--no_scale", action="store_true",
        help="Disable scaling prior to PCA (enabled by default)")
    refine_parser.add_argument("--npcs", type=int, default=50,
        help="Number of principal components used for refinement")
    refine_parser.add_argument("--neighbors", type=int, default=15,
        help="Number of nearest neighbours used to build the graph")
    refine_parser.add_argument("--theta", type=float, default=1,
        help="Harmony theta parameter (use lower values for subsets to avoid overcorrection)")
    refine_parser.add_argument("--resolution", type=float, default=0.5,
        help="Fixed Leiden resolution for refined clustering")
    refine_parser.add_argument("--auto_resolution", action="store_true",
        help="Automatically determine optimal Leiden resolution")
    refine_parser.add_argument("--res_min", type=float, default=0.1,
        help="Lower bound for automatic resolution scan")
    refine_parser.add_argument("--res_max", type=float, default=1.5,
        help="Upper bound for automatic resolution scan")
    refine_parser.add_argument("--res_steps", type=int, default=10,
        help="Number of resolution values tested")
    refine_parser.add_argument("--fine_width", type=float, default=0.2,
        help="Width of fine search window around best coarse resolution")
    refine_parser.add_argument("--nmarkers", type=int, default=100,
        help="Number of top marker genes exported per cluster")
    refine_parser.add_argument("--n_repeats", type=int, default=10,
        help="Number of repeated clustering runs used to estimate cluster stability")
    refine_parser.add_argument("--min_cells_per_group", type=int, default=500,
        help="Minimum number of cells required to refine a cluster group")
    refine_parser.add_argument("--min_cells_per_gene", "--min_cells", dest="min_cells_per_gene", type=int, default=3,
        help="Minimum number of cells a gene must be expressed in")
    refine_parser.add_argument("--max_refine_depth", type=int, default=2,
        help="Maximum recusrive refinement depth")
    refine_parser.add_argument("--stability_threshold", type=float, default=0.9,
        help="Threshold for unstable clusters")
    refine_parser.add_argument("--min_cells_per_cluster", type=int, default=50,
        help="Minimum number of cells permitted per cluster")
    refine_parser.add_argument("--marker_strength_threshold", type=float, default=1.0,
        help="logFC threshold")
    refine_parser.add_argument("--skip-integration", action="store_true",
        help="Skip Harmony and use PCA directly.")

    # ------------- Map to atlas ---------
    map_parser = subparsers.add_parser("map",
        description="""
        Map query (e.g. organoid) data onto a reference atlas using scANVI.

        This command performs:
          - gene space harmonisation between reference and query
          - semi-supervised label transfer using scANVI
          - confidence estimation (max probability and entropy)
          - rejection of low-confidence predictions
          - graph-based label smoothing
          - propagation of lineage annotations
          - generation of diagnostic plots (UMAP + confidence)

        Outputs:
          - mapped AnnData object with:
              * raw, smoothed and final predicted labels
              * prediction confidence and entropy scores
              * mapped lineage annotations
          - diagnostic plots in scribble/plots/:
              * predicted label UMAP
              * confidence UMAP
              * entropy (uncertainty) UMAP
              * raw prediction UMAP
          - summary tables in scribble/tables/
        """
    )
    map_parser.add_argument("--project_dir", required=True,
        help="Project directory containing scribble/ subfolders")
    map_parser.add_argument("--reference", required=True,
        help="Reference atlas AnnData (.h5ad), annotated using Scribble")
    map_parser.add_argument("--query", required=True,
        help="Query AnnData (.h5ad), e.g. organoid data processed with Scribble")
    map_parser.add_argument("--label_key", type=str, default="cell_type_major",
        help="Reference annotation key to transfer (default: 'cell_type_major')")
    map_parser.add_argument("--lineage_key", type=str, default="lineage",
        help="Reference lineage key to propagate (default: 'lineage')")
    map_parser.add_argument("--confidence_threshold", type=float, default=0.5,
        help=("Minimum prediction confidence for label assignment. "
            "Cells below this threshold are labelled as 'Unassigned' "))
    map_parser.add_argument("--nieghbors", type=int, default=15,
        help="Number of nearest neighbours used to smoothing labels")

    # -----------Add annotations ---------
    annotate_parser = subparsers.add_parser(
        "annotate",
        description="""
    Annotate a refined Scribble atlas using a worksheet called
    'annotations' from an Excel workbook.

    The worksheet must contain a 'refine_cluster' column matching
    the cluster labels stored in adata.obs['refine_label'].

    All additional columns are imported automatically and appended
    to adata.obs.

    An annotated AnnData file is written with the suffix
    '_annotated.h5ad'.

    A UMAP coloured by 'cell_type_major' is generated if that
    annotation column is present.

    Optional marker visualisation can be generated by specifying
    an annotation column containing semicolon-separated marker
    genes (e.g. 'key_markers' or 'canonical_markers').
    """
    )
    annotate_parser.add_argument("--project_dir", required=True,
        help="Project directory containing the Scribble output folders")
    annotate_parser.add_argument("--input", required=True,
        help="Input AnnData (.h5ad) file generated by Scribble")
    annotate_parser.add_argument("--annotations", required=True,
        help=("Excel workbook containing an 'annotations' worksheet "
            "with a required 'refine_cluster' column"))
    annotate_parser.add_argument("--plot_markers", default=None,
        help=("Optional annotation column containing semicolon-separated "
            "marker genes. Generates dot plots grouped by "
            "'cell_type_major' and UMAP expression plots for each marker "
            "(e.g. key_markers or canonical_markers)."))

    args = parser.parse_args()

    # ---------- ROUTING ----------
    COMMANDS = {
        "import": ("scribble.import_data", "run_import"),
        "mt": ("scribble.mt_qc", "run_mt_qc"),
        "dbl": ("scribble.dbl_qc", "run_dbl_qc"),
        "pca": ("scribble.pca_qc", "run_pca_qc"),
        "filter": ("scribble.filter", "run_filter"),
        "preintegration": ("scribble.preintegration", "run_preintegration"),
        "harmony": ("scribble.harmony", "run_harmony"),
        "cluster": ("scribble.cluster", "run_cluster"),
        "evaluate": ("scribble.evaluate", "run_evaluate"),
        "refine": ("scribble.refine", "run_refine"),
        "map": ("scribble.map", "run_map"),
        "annotate": ("scribble.annotate", "run_annotate")
    }
    module_path, func_name = COMMANDS[args.command]

    module = __import__(module_path, fromlist=[func_name])
    func = getattr(module, func_name)
    func(args)

if __name__ == "__main__":
    main()
