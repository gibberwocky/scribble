#!/usr/bin/env python

from pathlib import Path

def classify_cluster(row, thresholds):
    n = row["n_cells"]
    stability = row["mean_stability"]
    entropy = row["sample_entropy"]

    if n < thresholds["min_cells"]:
        return "keep", "small_cluster", "low_priority"

    if stability < thresholds["low_stability"]:
        return "trajectory", "low_stability_continuum", "medium"

    if n > thresholds["large_cells"] and stability < thresholds["high_stability"]:
        return "subset", "heterogeneous_large_cluster", "high"

    if entropy < thresholds["low_entropy"]:
        return "flag_bias", "sample_specific_cluster", "medium"

    return "keep", "well_defined_cluster", "low"


def run_evaluate(args):
    import scanpy as sc
    import numpy as np
    import pandas as pd
    import random
    from scribble.import_data import setup_environment

    PROJECT_DIR = Path(args.project_dir)
    PLOT_DIR = PROJECT_DIR / "scribble/plots"
    TABLE_DIR = PROJECT_DIR / "scribble/tables"
    setup_environment(sc, np, random, PLOT_DIR)

    input_file = Path(args.input)
    output_file = input_file.with_name(f"{input_file.stem}_decisions.tsv")

    print(f"Loading cluster summary: {input_file}")
    df = pd.read_csv(input_file, sep="\t")

    # ----------------------------
    # Thresholds (configurable)
    # ----------------------------
    thresholds = {
        "min_cells": args.min_cells,          # 200
        "large_cells": args.large_cells,      # 800
        "low_stability": args.low_stability,  # 0.75
        "high_stability": args.high_stability,# 0.95
        "low_entropy": args.low_entropy       # 0.5
    }

    decisions = []

    for _, row in df.iterrows():
        action, reason, priority = classify_cluster(row, thresholds)

        detail = (
            f"n={row['n_cells']}; "
            f"stability={row['mean_stability']:.2f}; "
            f"entropy={row['sample_entropy']:.2f}"
        )

        decisions.append({
            "cluster": row["cluster"],
            "action": action,
            "reason": reason,
            "detail": detail,
            "priority": priority
        })

    out_df = pd.DataFrame(decisions)

    print(out_df["action"].value_counts())

    print(f"Saving decisions → {output_file}")
    out_df.to_csv(output_file, sep="\t", index=False)
