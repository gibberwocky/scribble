#!/usr/bin/env python

import argparse

def main():
    parser = argparse.ArgumentParser(prog="scribble")
    subparsers = parser.add_subparsers(dest="command", required=True)

    # -------- IMPORT COMMAND --------
    import_parser = subparsers.add_parser("import")
    import_parser.add_argument("--project_dir", required=True)
    import_parser.add_argument("--metadata_file", required=True)

    # ------------- MT QC ------------
    mt_parser = subparsers.add_parser("mt")
    mt_parser.add_argument("--project_dir", required=True)
    mt_parser.add_argument("--input", required=True)
    mt_parser.add_argument("--nmads", type=float, default=3)

    args = parser.parse_args()

    # ---------- ROUTING ----------
    if args.command == "import":
        from scribble.import_data import run_import
        run_import(args)

    elif args.command == "mt":
        from scribble.mt_qc import run_mt_qc
        run_mt_qc(args)

if __name__ == "__main__":
    main()
