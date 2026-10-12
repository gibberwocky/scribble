#!/usr/bin/env python
#
import argparse
import traceback
from pathlib import Path
import re
import scanpy as sc
import pertpy as pt
import pandas as pd
import patsy
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import json

def safe_filename(value):
    value = str(value).strip()
    value = re.sub(r"[^\w.-]+", "_", value)
    value = re.sub(r"_+", "_", value)
    return value.strip("_.")


def get_formula_variables(formula):
    """
    Extract simple variable names from a Formulaic-style
    design formula.

    Intended for validation of formulas such as:
      ~ donor + treatment
      ~ treatment
      ~ donor + treatment + batchInfo

    Formulaic itself remains responsible for actual parsing.
    """

    tokens = re.findall(
        r"\b[A-Za-z_]\w*\b",
        formula
    )

    reserved = {
        "C",
        "I",
    }

    return {
        token
        for token in tokens
        if token not in reserved
    }

def plot_da_volcano(
    milo_results,
    outpath,
    fdr_threshold=0.05,
    nominal_fdr_threshold=0.25,
    min_abs_logFC=0.5,
):
    """
    Volcano plot showing:

      grey   = not nominal
      orange = nominal
      red    = significant
    """

    df = milo_results.copy()

    df["minus_log10_SpatialFDR"] = (
        -np.log10(
            np.clip(
                df["SpatialFDR"],
                1e-300,
                None,
            )
        )
    )

    significant = (
        (
            df["SpatialFDR"]
            < fdr_threshold
        )
        &
        (
            df["logFC"]
            .abs()
            >= min_abs_logFC
        )
    )

    nominal = (
        (
            df["SpatialFDR"]
            < nominal_fdr_threshold
        )
        &
        (
            df["logFC"]
            .abs()
            >= min_abs_logFC
        )
        &
        (~significant)
    )

    background = (
        ~(significant | nominal)
    )

    fig, ax = plt.subplots(
        figsize=(8, 6)
    )

    ax.scatter(
        df.loc[
            background,
            "logFC",
        ],
        df.loc[
            background,
            "minus_log10_SpatialFDR",
        ],
        s=8,
        alpha=0.3,
        color="lightgrey",
        label=(
            f"FDR ≥ "
            f"{nominal_fdr_threshold}"
        ),
        rasterized=True,
    )

    ax.scatter(
        df.loc[
            nominal,
            "logFC",
        ],
        df.loc[
            nominal,
            "minus_log10_SpatialFDR",
        ],
        s=12,
        alpha=0.8,
        color="orange",
        label=(
            f"Nominal "
            f"(FDR < "
            f"{nominal_fdr_threshold})"
        ),
        rasterized=True,
    )

    ax.scatter(
        df.loc[
            significant,
            "logFC",
        ],
        df.loc[
            significant,
            "minus_log10_SpatialFDR",
        ],
        s=18,
        alpha=0.9,
        color="red",
        label=(
            f"Significant "
            f"(FDR < "
            f"{fdr_threshold})"
        ),
        rasterized=True,
    )

    ax.axvline(
        min_abs_logFC,
        linestyle="--",
        color="black",
    )

    ax.axvline(
        -min_abs_logFC,
        linestyle="--",
        color="black",
    )

    # Significant threshold
    ax.axhline(
        -np.log10(
            fdr_threshold
        ),
        linestyle="--",
        color="red",
        linewidth=1,
    )

    # Nominal threshold
    ax.axhline(
        -np.log10(
            nominal_fdr_threshold
        ),
        linestyle=":",
        color="orange",
        linewidth=1,
    )

    ax.set_xlabel(
        "Neighborhood logFC"
    )

    ax.set_ylabel(
        "-log10(SpatialFDR)"
    )

    ax.set_title(
        "Milo differential abundance"
    )

    ax.legend()

    fig.tight_layout()

    fig.savefig(
        outpath,
        dpi=300,
    )

    plt.close(fig)


def plot_nhood_graph(
    milo,
    mdata,
    outpath,
    fdr_threshold,
    min_abs_logFC,
):
    """
    Save Milo neighbourhood graph.
    """

    fig = milo.plot_nhood_graph(
        mdata,
        padj_threshold=fdr_threshold,
        min_logFC=min_abs_logFC,
        return_fig=True,
    )

    fig.savefig(
        outpath,
        dpi=300,
        bbox_inches="tight",
    )

    plt.close(fig)


def plot_celltype_summary(
    cell_type_summary,
    outpath,
    positive_col,
    negative_col,
    title,
    xlabel,
):
    """
    Horizontal barplot of positive and negative
    DA neighbourhood counts per annotation.

    positive_col:
        column containing positive counts

    negative_col:
        column containing negative counts
    """

    plot_df = cell_type_summary.copy()

    total = (
        plot_df[positive_col]
        +
        plot_df[negative_col]
    )

    plot_df = (
        plot_df.assign(
            total=total
        )
        .sort_values(
            "total",
            ascending=False,
        )
        .drop(
            columns="total"
        )
    )

    if plot_df.empty:
        return

    fig_height = max(
        4,
        0.35 * len(plot_df),
    )

    fig, ax = plt.subplots(
        figsize=(8, fig_height)
    )

    y = plot_df.iloc[::-1]

    ax.barh(
        y["nhood_annotation"],
        y[positive_col],
        label=(
            f"{CONTRAST_TEST}-higher"
        ),
        color="tab:red",
    )

    ax.barh(
        y["nhood_annotation"],
        -y[negative_col],
        label=(
            f"{CONTRAST_REFERENCE}-higher"
        ),
        color="tab:blue",
    )

    ax.axvline(
        0,
        color="black",
        linewidth=1,
    )

    ax.set_xlabel(
        xlabel
    )

    ax.set_ylabel(
        ""
    )

    ax.set_title(
        title
    )

    ax.legend()

    fig.tight_layout()

    fig.savefig(
        outpath,
        dpi=300,
    )

    plt.close(fig)


def plot_annotation_logFC(
    milo_results,
    outpath,
):
    """
    Boxplots of neighbourhood logFC
    grouped by annotation.
    """

    plot_df = milo_results.copy()

    counts = (
        plot_df["nhood_annotation"]
        .value_counts()
    )

    order = counts.index.tolist()

    fig_height = max(
        5,
        0.35 * len(order),
    )

    fig, ax = plt.subplots(
        figsize=(10, fig_height)
    )

    sns.boxplot(
        data=plot_df,
        x="logFC",
        y="nhood_annotation",
        order=order,
        fliersize=0,
        ax=ax,
    )

    ax.axvline(
        0,
        linestyle="--",
    )

    ax.set_xlabel(
        "Neighborhood logFC"
    )

    ax.set_ylabel(
        ""
    )

    ax.set_title(
        "Neighborhood abundance shifts by annotation"
    )

    fig.tight_layout()

    fig.savefig(
        outpath,
        dpi=300,
    )

    plt.close(fig)


def plot_annotation_median_logFC(
    cell_type_summary,
    outpath,
):
    """
    Plot median neighbourhood logFC for each annotation.

    Positive:
        neighbourhoods enriched in TEST condition

    Negative:
        neighbourhoods depleted in TEST condition
        (or enriched in REFERENCE condition)
    """

    plot_df = (
        cell_type_summary
        .sort_values(
            "median_logFC",
            ascending=True,
        )
        .copy()
    )

    if plot_df.empty:
        return

    plot_df["label"] = [
        f"{annotation} (n={n})"
        for annotation, n in zip(
            plot_df["nhood_annotation"].astype(str),
            plot_df["n_neighborhoods"],
        )
    ]

    fig_height = max(
        4,
        0.35 * len(plot_df),
    )

    fig, ax = plt.subplots(
        figsize=(8, fig_height)
    )

    colors = np.where(
        plot_df["median_logFC"] > 0,
        "tab:red",
        "tab:blue",
    )

    ax.barh(
        plot_df["label"],
        plot_df["median_logFC"],
        color=colors,
    )

    ax.axvline(
        0,
        color="black",
        linestyle="--",
        linewidth=1,
    )

    ax.set_xlabel(
        (
            f"Median neighbourhood logFC "
            f"({CONTRAST_TEST} vs "
            f"{CONTRAST_REFERENCE})"
        )
    )

    ax.set_ylabel("")

    ax.set_title(
        "Median neighbourhood abundance shift by annotation"
    )

    fig.tight_layout()

    fig.savefig(
        outpath,
        dpi=300,
        bbox_inches="tight",
    )

    plt.close(fig)

# --------------------------------------------------
# CLI parameters
# --------------------------------------------------
parser = argparse.ArgumentParser(
    description=("Perform Milo differential abundance testing using PertPy.")
)
parser.add_argument("--project_dir", type=str, required=True)
parser.add_argument("--input_file", type=str, required=True)
parser.add_argument("--annotation", type=str, required=True)
parser.add_argument("--analysis_name", type=str, required=True,
    help="Unique name for this Milo analysis.")
parser.add_argument("--samples", nargs="+", default=None,
    help=("Optional list of sample IDs to include. "
        "If omitted, all samples containing the contrast "
        "levels are considered.")
)
parser.add_argument("--formula", type=str, required=True,
    help="Milo design formula, e.g. '~ donor + treatment' or '~ treatment'.")
parser.add_argument("--contrast", nargs=3, required=True,
    metavar=("VARIABLE", "TEST_LEVEL", "REFERENCE_LEVEL",),
    help=("Contrast description. Example: --contrast treatment IL17A veh")
)
parser.add_argument("--solver", choices=["pydeseq2", "edger"], default="pydeseq2",
    help="Differential abundance solver.")
parser.add_argument("--embedding", required=True,
    help=(
        "Embedding used for neighbour graph "
        "construction, e.g. X_pca_harmony."
    ),
)
parser.add_argument("--basis", default="X_umap",
    help="Embedding used for Milo graph visualisation.")
parser.add_argument("--n_neighbors", type=int, default=30,
    help="Number of neighbours for Scanpy graph construction.")
parser.add_argument("--nhood_prop", type=float, default=0.1,
    help="Proportion of cells used as neighbourhood seeds.")
parser.add_argument("--fdr_threshold", type=float, default=0.05,
    help="SpatialFDR threshold for significant neighbourhood exports.")
parser.add_argument("--nominal_fdr_threshold", type=float, default=0.25,
    help="Relaxed SpatialFDR threshold used for exploratory / nominal neighbourhood summaries.")
parser.add_argument("--min_abs_logFC", type=float, default=0.5,
    help="Minimum absolute neighbourhood logFC for significant exports.")
parser.add_argument("--random_seed", type=int, default=0,
    help="Random seed used by Milo.")
args = parser.parse_args()


CONTRAST_VARIABLE = args.contrast[0]
CONTRAST_TEST = args.contrast[1]
CONTRAST_REFERENCE = args.contrast[2]
FORMULA_VARIABLES = get_formula_variables(args.formula)
DONOR_IN_MODEL = ("donor" in FORMULA_VARIABLES)
REQUIRED_MODEL_COLUMNS = (FORMULA_VARIABLES | {CONTRAST_VARIABLE, "sample"})

if (
    CONTRAST_VARIABLE
    not in FORMULA_VARIABLES
):
    raise ValueError(
        f"Contrast variable "
        f"'{CONTRAST_VARIABLE}' "
        "is not present in design formula "
        f"'{args.formula}'."
    )

if (
    CONTRAST_TEST
    == CONTRAST_REFERENCE
):
    raise ValueError(
        "Contrast test and reference "
        "levels must differ."
    )



RUN_INFO = {
    "analysis_name": args.analysis_name,
    "formula": args.formula,
    "contrast_variable": (
        CONTRAST_VARIABLE
    ),
    "contrast_test": (
        CONTRAST_TEST
    ),
    "contrast_reference": (
        CONTRAST_REFERENCE
    ),
}


# Create directory structure
ANALYSIS_NAME = safe_filename(args.analysis_name)
OUTDIR = (Path(args.project_dir) / "scribble" / "DA" / ANALYSIS_NAME)
NHOOD_DIR = (OUTDIR / "neighbourhoods" / args.annotation)
CELLTYPE_DIR = (OUTDIR / "cell_types" / args.annotation)
PLOTS_DIR = (OUTDIR / "plots" / args.annotation)
INSPECTION_DIR = (OUTDIR / "inspection" / args.annotation)
NHOOD_DIR.mkdir(exist_ok=True, parents=True)
CELLTYPE_DIR.mkdir(exist_ok=True, parents=True)
PLOTS_DIR.mkdir(exist_ok=True, parents=True)
INSPECTION_DIR.mkdir(parents=True, exist_ok=True)


# Run manifest
run_manifest = {
    "analysis_name": args.analysis_name,
    "input_file": str(
        Path(args.input_file)
        .resolve()
    ),
    "annotation": args.annotation,
    "samples": args.samples,
    "formula": args.formula,
    "formula_variables": sorted(FORMULA_VARIABLES),
    "donor_in_model": DONOR_IN_MODEL,
    "contrast_variable": (CONTRAST_VARIABLE),
    "contrast_test": (CONTRAST_TEST),
    "contrast_reference": (CONTRAST_REFERENCE),
    "solver": args.solver,
    "embedding": args.embedding,
    "basis": args.basis,
    "n_neighbors": args.n_neighbors,
    "nhood_prop": args.nhood_prop,
    "fdr_threshold": args.fdr_threshold,
    "nominal_fdr_threshold": args.nominal_fdr_threshold,
    "min_abs_logFC": args.min_abs_logFC,
    "random_seed": args.random_seed
}

with open(
    OUTDIR /
    "run_manifest.json",
    "w",
) as handle:

    json.dump(
        run_manifest,
        handle,
        indent=2,
    )






# Import adata
adata = sc.read_h5ad(args.input_file)

# Check samples requested are present
if "sample" not in adata.obs.columns:
    raise ValueError(
        "Required column 'sample' "
        "is not present in adata.obs."
    )


# --------------------------------------------------
# Select requested samples
# --------------------------------------------------

if args.samples is not None:

    requested_samples = set(
        args.samples
    )

    available_samples = set(
        adata.obs[
            "sample"
        ]
        .astype(str)
        .unique()
    )

    missing_samples = (
        requested_samples
        -
        available_samples
    )

    if missing_samples:
        raise ValueError(
            "Requested samples not present "
            "in AnnData: "
            f"{sorted(missing_samples)}"
        )

    adata_da = adata[
        adata.obs[
            "sample"
        ]
        .astype(str)
        .isin(
            requested_samples
        )
    ].copy()

else:

    adata_da = (
        adata.copy()
    )




# --------------------------------------------------
# Keep only requested contrast levels
# --------------------------------------------------
if (
    CONTRAST_VARIABLE
    not in adata_da.obs.columns
):
    raise ValueError(
        f"Contrast variable "
        f"'{CONTRAST_VARIABLE}' "
        "not present in adata.obs."
    )

missing_columns = sorted(
    REQUIRED_MODEL_COLUMNS
    - set(adata_da.obs.columns)
)

if missing_columns:
    raise ValueError(
        "Missing required columns in "
        "adata.obs:\n"
        + "\n".join(missing_columns)
    )

adata_da = adata_da[
    adata_da.obs[
        CONTRAST_VARIABLE
    ]
    .astype(str)
    .isin(
        [
            CONTRAST_TEST,
            CONTRAST_REFERENCE,
        ]
    )
].copy()

# --------------------------------------------------
# Validate contrast levels
# --------------------------------------------------

observed_contrast_levels = set(
    adata_da.obs[
        CONTRAST_VARIABLE
    ]
    .astype(str)
)

expected_contrast_levels = {
    CONTRAST_TEST,
    CONTRAST_REFERENCE,
}

unexpected_levels = (
    observed_contrast_levels
    -
    expected_contrast_levels
)

if unexpected_levels:
    raise ValueError(
        f"Unexpected {CONTRAST_VARIABLE} levels "
        "after contrast subsetting: "
        f"{sorted(unexpected_levels)}"
    )

missing_levels = (
    expected_contrast_levels
    -
    observed_contrast_levels
)

if missing_levels:
    raise ValueError(
        f"Missing required {CONTRAST_VARIABLE} levels: "
        f"{sorted(missing_levels)}"
    )

if adata_da.n_obs == 0:
    raise ValueError(
        "No observations remain after sample "
        "and contrast filtering."
    )

if (
    args.annotation
    not in adata_da.obs.columns
):
    raise ValueError(
        f"Annotation column "
        f"'{args.annotation}' "
        "not found in adata.obs."
    )

# Get sample cell counts
sample_counts = (
    adata_da.obs
    .groupby([args.annotation, "sample"], observed=True)
    .size()
    .unstack(fill_value=0)
)
file = (OUTDIR / f"{Path(args.input_file).stem}_{args.annotation}_counts.csv")
sample_counts.to_csv(file)


sample_metadata = (
    adata_da.obs[
        [
            "sample",
            CONTRAST_VARIABLE,
            *sorted(FORMULA_VARIABLES),
        ]
    ]
    .drop_duplicates()
    .sort_values("sample")
)

sample_metadata.to_csv(
    INSPECTION_DIR /
    "sample_metadata.csv",
    index=False,
)



# --------------------------------------------------
# Verify metadata present for all samples
# --------------------------------------------------

sample_design_cols = [
    "sample"
]

for col in sorted(
    REQUIRED_MODEL_COLUMNS
):

    if (
        col in adata_da.obs.columns
        and
        col not in sample_design_cols
    ):
        sample_design_cols.append(
            col
        )

sample_design = (
    adata_da.obs[
        sample_design_cols
    ]
    .drop_duplicates()
    .sort_values(
        "sample"
    )
)

print(
    sample_design.to_string(
        index=False
    )
)

sample_design.to_csv(
    INSPECTION_DIR /
    "sample_design.csv",
    index=False,
)

if (
    sample_design["sample"]
    .duplicated()
    .any()
):
    raise ValueError(
        "Sample design table contains "
        "duplicate sample metadata rows."
    )

print()
print("Sample counts per contrast level:")

print(
    sample_design[
        CONTRAST_VARIABLE
    ]
    .value_counts(
        dropna=False
    )
)

if DONOR_IN_MODEL:

    donor_summary = (
        sample_design
        .groupby(
            "donor",
            observed=True,
        )[CONTRAST_VARIABLE]
        .agg(
            lambda x:
            "|".join(
                sorted(
                    x.astype(str)
                    .unique()
                )
            )
        )
        .reset_index()
    )

    donor_summary.to_csv(
        INSPECTION_DIR /
        "donor_design_summary.csv",
        index=False,
    )


# --------------------------------------------------
# Validate sample-level metadata
# --------------------------------------------------

metadata_cols = [
    "sample",
] + sorted(
    REQUIRED_MODEL_COLUMNS
    - {"sample"}
)


for col in metadata_cols:

    if col not in adata_da.obs.columns:
        continue

    n_values = (
        adata_da.obs
        .groupby(
            "sample",
            observed=True,
        )[col]
        .nunique(
            dropna=False
        )
    )

    bad = n_values[
        n_values != 1
    ]

    if len(bad) > 0:
        raise ValueError(
            f"Inconsistent {col} metadata for samples: "
            f"{bad.index.tolist()}"
        )


# --------------------------------------------------
# Check required model metadata is not missing
# --------------------------------------------------

for col in sorted(
    REQUIRED_MODEL_COLUMNS
):

    if col not in adata_da.obs.columns:
        raise ValueError(
            f"Required model metadata column "
            f"'{col}' is absent from adata.obs."
        )

    missing_samples = (
        adata_da.obs
        .groupby(
            "sample",
            observed=True,
        )[col]
        .apply(
            lambda x:
                x.isna().all()
                or
                (
                    x.astype(str)
                    .str.strip()
                    .eq("")
                    .all()
                )
        )
    )

    missing_samples = (
        missing_samples[
            missing_samples
        ]
        .index
        .tolist()
    )

    if missing_samples:
        raise ValueError(
            f"Missing required model metadata "
            f"'{col}' for samples: "
            f"{missing_samples}"
        )



# --------------------------------------------------
# Update manifest
# --------------------------------------------------

selected_samples = sorted(
    adata_da.obs["sample"]
    .astype(str)
    .unique()
    .tolist()
)

run_manifest["selected_samples"] = (selected_samples)
run_manifest["n_selected_samples"] = (len(selected_samples))
run_manifest["n_cells"] = int(adata_da.n_obs)
run_manifest["n_genes"] = int(adata_da.n_vars)

with open(
    OUTDIR / "run_manifest.json",
    "w",
) as handle:

    json.dump(
        run_manifest,
        handle,
        indent=2,
    )



# --------------------------------------------------
# Full-rank design-matrix validation
# --------------------------------------------------

design_matrix = patsy.dmatrix(
    args.formula,
    sample_design,
    return_type="dataframe",
)

rank = np.linalg.matrix_rank(
    design_matrix.values
)

if rank < design_matrix.shape[1]:
    raise ValueError(
        "Design matrix is not full rank.\n"
        f"Rank: {rank}\n"
        f"Columns: {design_matrix.shape[1]}"
    )



# --------------------------------------------------
# Build neighbour graph
# --------------------------------------------------
sc.pp.neighbors(adata_da, use_rep=args.embedding)

# --------------------------------------------------
# Milo
# --------------------------------------------------
milo = pt.tl.Milo()
mdata = milo.load(adata_da)

print("Creating neighbourhoods...")
milo.make_nhoods(mdata, prop=args.nhood_prop)

print("Counting neighbourhoods...")
milo.count_nhoods(mdata, sample_col="sample")

print("Running DA...")
print()
print("Milo design:")
print(
    f"  Formula: {args.formula}"
)

print(
    f"  Contrast: "
    f"{CONTRAST_VARIABLE} "
    f"{CONTRAST_TEST} "
    f"vs "
    f"{CONTRAST_REFERENCE}"
)

print(
    f"  Solver: {args.solver}"
)
print()
milo.da_nhoods(mdata, design=args.formula, solver=args.solver)
milo.annotate_nhoods(mdata, anno_col=args.annotation)
milo_results = (mdata.mod["milo"].var.copy())
milo_results.index.name = ("neighborhood_id")
milo_results.reset_index(inplace=True)

print()
print("DA summary")
print(
    milo_results["SpatialFDR"]
    .describe()
)

print()
print(
    milo_results[
        [
            "logFC",
            "PValue",
            "SpatialFDR",
            "nhood_annotation",
        ]
    ]
    .sort_values(
        "SpatialFDR"
    )
    .head(20)
)

milo_results.to_csv(
    NHOOD_DIR /
    "all_neighborhoods.csv",
    index=False,
)

significant_nhoods = (
    milo_results[
        (
            milo_results["SpatialFDR"]
            < args.fdr_threshold
        )
        &
        (
            milo_results["logFC"]
            .abs()
            >= args.min_abs_logFC
        )
    ]
    .copy()
)

nominal_nhoods = (
    milo_results[
        (
            milo_results["SpatialFDR"]
            < args.nominal_fdr_threshold
        )
        &
        (
            milo_results["logFC"]
            .abs()
            >= args.min_abs_logFC
        )
    ]
    .copy()
)

significant_nhoods.to_csv(
    NHOOD_DIR /
    "significant_neighborhoods.csv",
    index=False,
)

nominal_nhoods.to_csv(
    NHOOD_DIR /
    "nominal_neighborhoods.csv",
    index=False,
)

significant_positive = (
    significant_nhoods[
        significant_nhoods["logFC"] > 0
    ]
)

significant_negative = (
    significant_nhoods[
        significant_nhoods["logFC"] < 0
    ]
)

nominal_positive = (
    nominal_nhoods[
        nominal_nhoods["logFC"] > 0
    ]
)

nominal_negative = (
    nominal_nhoods[
        nominal_nhoods["logFC"] < 0
    ]
)

significant_positive.to_csv(
    NHOOD_DIR /
    f"significant_{safe_filename(CONTRAST_TEST)}_up.csv",
    index=False,
)

significant_negative.to_csv(
    NHOOD_DIR /
    f"significant_{safe_filename(CONTRAST_TEST)}_down.csv",
    index=False,
)

nominal_positive.to_csv(
    NHOOD_DIR /
    f"nominal_{safe_filename(CONTRAST_TEST)}_up.csv",
    index=False,
)

nominal_negative.to_csv(
    NHOOD_DIR /
    f"nominal_{safe_filename(CONTRAST_TEST)}_down.csv",
    index=False,
)


annotated_nhoods = (
    milo_results
    .sort_values(
        [
            "nhood_annotation",
            "SpatialFDR",
        ]
    )
)

annotated_nhoods.to_csv(
    NHOOD_DIR /
    "annotated_neighborhoods.csv",
    index=False,
)

cell_type_summary = (
    milo_results
    .groupby(
        "nhood_annotation",
        observed=True,
    )
    .agg(
        n_neighborhoods=(
            "nhood_annotation",
            "size",
        ),

        n_significant=(
            "SpatialFDR",
            lambda x:
            (
                x < args.fdr_threshold
            ).sum(),
        ),

        n_significant_positive=(
            "logFC",
            lambda x:
            (
                (x > 0)
                &
                (
                    milo_results.loc[
                        x.index,
                        "SpatialFDR"
                    ]
                    < args.fdr_threshold
                )
            ).sum(),
        ),

        n_significant_negative=(
            "logFC",
            lambda x:
            (
                (x < 0)
                &
                (
                    milo_results.loc[
                        x.index,
                        "SpatialFDR"
                    ]
                    < args.fdr_threshold
                )
            ).sum(),
        ),

        n_nominal=(
            "SpatialFDR",
            lambda x:
            (
                x < args.nominal_fdr_threshold
            ).sum(),
        ),

        n_nominal_positive=(
            "logFC",
            lambda x:
            (
                (x > 0)
                &
                (
                    milo_results.loc[
                        x.index,
                        "SpatialFDR"
                    ]
                    < args.nominal_fdr_threshold
                )
            ).sum(),
        ),

        n_nominal_negative=(
            "logFC",
            lambda x:
            (
                (x < 0)
                &
                (
                    milo_results.loc[
                        x.index,
                        "SpatialFDR"
                    ]
                    < args.nominal_fdr_threshold
                )
            ).sum(),
        ),

        median_logFC=(
            "logFC",
            "median",
        ),

        median_SpatialFDR=(
            "SpatialFDR",
            "median",
        ),

        mean_annotation_fraction=(
            "nhood_annotation_frac",
            "mean",
        ),
    )
    .reset_index()
)

cell_type_summary.to_csv(
    CELLTYPE_DIR /
    "cell_type_summary.csv",
    index=False,
)


for annotation, sub in (
    milo_results.groupby(
        "nhood_annotation",
        observed=True,
        sort=True,
    )
):

    outfile = (
        CELLTYPE_DIR
        / f"{safe_filename(annotation)}.csv"
    )

    sub.sort_values(
        "SpatialFDR"
    ).to_csv(
        outfile,
        index=False,
    )

#mdata.write(OUTDIR / "milo_results.h5mu")


# --------------------------------------------------
# Plots
# --------------------------------------------------

milo.build_nhood_graph(
    mdata,
    basis=args.basis,
)

plot_da_volcano(
    milo_results=milo_results,
    outpath=(
        PLOTS_DIR /
        "DA_volcano.png"
    ),
    fdr_threshold=args.fdr_threshold,
    nominal_fdr_threshold=args.nominal_fdr_threshold,
    min_abs_logFC=args.min_abs_logFC,
)

plot_nhood_graph(
    milo=milo,
    mdata=mdata,
    outpath=(
        PLOTS_DIR /
        "nhood_graph.png"
    ),
    fdr_threshold=args.fdr_threshold,
    min_abs_logFC=args.min_abs_logFC,
)

plot_celltype_summary(
    cell_type_summary=cell_type_summary,
    outpath=(
        PLOTS_DIR /
        "celltype_summary_significant.png"
    ),
    positive_col="n_significant_positive",
    negative_col="n_significant_negative",
    title=(
        "Significant differential abundance "
        "by annotation"
    ),
    xlabel=(
        "Significant neighbourhoods"
    ),
)

plot_celltype_summary(
    cell_type_summary=cell_type_summary,
    outpath=(
        PLOTS_DIR /
        "celltype_summary_nominal.png"
    ),
    positive_col="n_nominal_positive",
    negative_col="n_nominal_negative",
    title=(
        "Nominal differential abundance "
        "by annotation"
    ),
    xlabel=(
        "Nominal neighbourhoods"
    ),
)

plot_annotation_logFC(
    milo_results=milo_results,
    outpath=(
        PLOTS_DIR /
        "annotation_logFC.png"
    ),
)


plot_annotation_median_logFC(
    cell_type_summary=cell_type_summary,
    outpath=(
        PLOTS_DIR /
        "annotation_median_logFC.png"
    ),
)

# --------------------------------------------------
# Update manifest
# --------------------------------------------------

run_manifest[
    "n_neighborhoods"
] = int(
    milo_results.shape[0]
)

run_manifest[
    "n_significant_neighborhoods"
] = int(
    significant_nhoods.shape[0]
)

run_manifest[
    "n_annotations"
] = int(
    cell_type_summary.shape[0]
)

with open(
    OUTDIR /
    "run_manifest.json",
    "w",
) as handle:

    json.dump(
        run_manifest,
        handle,
        indent=2,
    )
