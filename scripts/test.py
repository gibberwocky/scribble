#!/usr/bin/env python
from pathlib import Path
import scanpy as sc
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.image as mpimg
import random
from scribble.import_data import setup_environment
import os
import argparse

parser = argparse.ArgumentParser()

parser.add_argument("--input_file", type=str, required=True)
parser.add_argument("--plot_file", type=str, required=True)
parser.add_argument("--markers", nargs="+", required=True)
parser.add_argument("--dotplot", action="store_true")
parser.add_argument("--umap", action="store_true")


args = parser.parse_args()


# Import adata
adata = sc.read(args.input_file)

if args.dotplot:
    sc.pl.dotplot(
        adata,
        args.markers,
        standard_scale="var",
        groupby="refine_label",
        dendrogram=False
    )
    plt.savefig(args.plot_file, dpi=300, bbox_inches="tight")
    plt.close()

if args.umap:
    for gene in args.markers:
        try:
            sc.pl.umap(
                adata,
                color=gene,
                cmap="Reds"
            )
        except KeyError:
            print(f"Warning: {gene} not found. Skipping.")
            continue

        out_file = Path(args.plot_file).with_name(
            f"{Path(args.plot_file).stem}_{gene}.png"
        )
        plt.savefig(out_file, dpi=300, bbox_inches="tight")
        plt.close()
