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

args = parser.parse_args()

#marker_dict = ["TTR", "CP", "PAX8", "HNF1B", "DNAH11", "MAP2"]



# Import adata
adata = sc.read(args.input_file)



sc.pl.dotplot(
    adata,
    args.markers,
    standard_scale="var",
    groupby="refine_label",
    dendrogram=False
)


plt.savefig(args.plot_file, dpi=300, bbox_inches="tight")
plt.close()
