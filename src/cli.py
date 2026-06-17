#!/usr/bin/env python

import argparse
from scribble.import_data import run_import

def main():
    parser = argparse.ArgumentParser(prog="scribble")
    subparsers = parser.add_subparsers(dest="command", required=True)

    # -------- IMPORT COMMAND --------
    import_parser = subparsers.add_parser("import")
    import_parser.add_argument("--project_dir", required=True)
    import_parser.add_argument("--metadata_file", required=True)

    args = parser.parse_args()

    # ---------- ROUTING ----------
    if args.command == "import":
        from import_data import run_import
        run_import(args)


if __name__ == "__main__":
    main()
