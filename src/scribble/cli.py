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
A single cell workflow CLI based on scapy for processing 10X data. Scribble expects the project directory to include:
- A directory containing cellranger count outputs (e.g.\033[34m <project>/cellranger/<sample>/outs/filtered_feature_bc_matrix\033[32m)
- A directory containing velocyto loom outputs (e.g.\033[34m <project>/velocyto/<sample>/<filename>.loom\033[32m)
- An Excel file containing sample metadata, which must include worksheet 'meta' containing column 'sample' whose values match the velocyto and cellranger <sample> folder names
\033[0m
GitHub: https://github.com/gibberwocky/scribble
""",
    formatter_class=argparse.RawTextHelpFormatter,
)
    subparsers = parser.add_subparsers(dest="command", required=True)

    # ---------- IMPORT COMMAND ----------
    import_parser = subparsers.add_parser("import")
    import_parser.add_argument("--project_dir", required=True,
        help='Directory containing subdirectories of cellranger and velocyto outputs')
    import_parser.add_argument("--metadata_file", required=True,
        help="Excel file containing metadata for samples, must include worksheet 'meta' with a 'sample' column")

    # -------------- MT QC ---------------
    mt_parser = subparsers.add_parser("mt")
    mt_parser.add_argument("--project_dir", required=True)
    mt_parser.add_argument("--input", required=True)
    mt_parser.add_argument("--nmads", type=int, default=3)

    # ------------ Doublet QC ------------
    dbl_parser = subparsers.add_parser("dbl")
    dbl_parser.add_argument("--project_dir", required=True)
    dbl_parser.add_argument("--input", required=True)
    dbl_parser.add_argument("--expected", type=float, default=0.07)

    # -------------- PCA QC --------------
    pca_parser = subparsers.add_parser("pca")
    pca_parser.add_argument("--project_dir", required=True)
    pca_parser.add_argument("--input", required=True)
    pca_parser.add_argument("--mingenes", type=int, default=200)
    pca_parser.add_argument("--maxgenes", type=int, default=6000)
    pca_parser.add_argument("--hvgs", type=int, default=3000)
    pca_parser.add_argument("--vmax", type=float, default=0.99)

    # -------------- Filter --------------
    filter_parser = subparsers.add_parser("filter")
    filter_parser.add_argument("--project_dir", required=True)
    filter_parser.add_argument("--input", required=True)
    filter_parser.add_argument("--mingenes", type=int, default=200,
        help="Minimum number of genes with non-zero counts in a given cell")
    filter_parser.add_argument("--maxgenes", type=int, default=6000,
        help="Maximum number of genes with non-zero counts in a given cell")
    filter_parser.add_argument("--filter_xlsx",
        help="Excel file containing filters for samples, must include worksheet 'filters' with a 'sample' column")

    # ---------- Preintegration ----------
    preintegration_parser = subparsers.add_parser("preintegration")
    preintegration_parser.add_argument("--project_dir", required=True)
    preintegration_parser.add_argument("--input", required=True)
    preintegration_parser.add_argument("--hvgs", type=int, default=3000)
    preintegration_parser.add_argument("--npcs", type=int, default=30)
    preintegration_parser.add_argument("--neighbors", type=int, default=15)
    preintegration_parser.add_argument("--batch", type=str, default="sample")
    preintegration_parser.add_argument("--vars", nargs="+", default=["sample"],
        help="UMAP colour variable(s)")
    preintegration_parser.add_argument("--no-scale", action="store_true",
        help="Apply scaling before PCA (default: on)")
    preintegration_parser.add_argument("--regress", nargs="+", default=None,
        help="Variables to regress out (e.g. total_counts pct_counts_mt)")

    # -------------- Harmony -------------
    harmony_parser = subparsers.add_parser("harmony")
    harmony_parser.add_argument("--project_dir", required=True)
    harmony_parser.add_argument("--input", required=True)
    harmony_parser.add_argument("--npcs", type=int, default=30)
    harmony_parser.add_argument("--neighbors", type=int, default=15)
    harmony_parser.add_argument("--batch", type=str, default="sample")
    harmony_parser.add_argument("--theta", type=int, default=2)
    harmony_parser.add_argument("--vars", nargs="+", default=["sample"],
        help="UMAP colour variable(s)")

    # -------------- Cluster -------------
    cluster_parser = subparsers.add_parser("cluster")
    cluster_parser.add_argument("--project_dir", required=True)
    cluster_parser.add_argument("--input", required=True)
    cluster_parser.add_argument("--embedding", default="X_pca",
        help="X_pca if no integration with Harmomy, otherwise X_pca_harmony")
    cluster_parser.add_argument("--resolution", type=float, default=1.0)
    cluster_parser.add_argument("--neighbors", type=int, default=15)
    cluster_parser.add_argument("--vars", nargs="+", default=["sample"],
        help="UMAP colour variable(s) to include in addition to Leiden cluster")
    cluster_parser.add_argument("--n_repeats", type=int, default=10)


    args = parser.parse_args()

    # ---------- ROUTING ----------
    if args.command == "import":
        from scribble.import_data import run_import
        run_import(args)

    elif args.command == "mt":
        from scribble.mt_qc import run_mt_qc
        run_mt_qc(args)

    elif args.command == "dbl":
        from scribble.dbl_qc import run_dbl_qc
        run_dbl_qc(args)

    elif args.command == "pca":
        from scribble.pca_qc import run_pca_qc
        run_pca_qc(args)

    elif args.command == "filter":
        from scribble.filter import run_filter
        run_filter(args)

    elif args.command == "preintegration":
        from scribble.preintegration import run_preintegration
        run_preintegration(args)

    elif args.command == "harmony":
        from scribble.harmony import run_harmony
        run_harmony(args)

    elif args.command == "cluster":
        from scribble.cluster import run_cluster
        run_cluster(args)

if __name__ == "__main__":
    main()
