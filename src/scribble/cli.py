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

if __name__ == "__main__":
    main()
