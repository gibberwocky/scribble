#!/usr/bin/env python

import argparse


def main():
    parser = argparse.ArgumentParser(prog="scribble")

    subparsers = parser.add_subparsers(dest="command", required=True)

    # -------- IMPORT COMMAND --------
    import_parser = subparsers.add_parser("import")

    import_parser.add_argument("--project_dir", required=True)
    import_parser.add_argument("--metadata_file", required=True)

    # later:
    # qc_parser = subparsers.add_parser("qc")
    # integrate_parser = subparsers.add_parser("integrate")

    args = parser.parse_args()

    # ---------- ROUTING ----------
    if args.command == "import":
        from src.import_data import run_import
        run_import(args)


if __name__ == "__main__":
    main()
