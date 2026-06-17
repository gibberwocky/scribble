#!/usr/bin/env python

from pathlib import Path


def compute_mt_outliers(adata, nmads):
    import numpy as np

    adata.obs["mt_outlier"] = False
    mt_stats = {}

    for sample in adata.obs["sample"].unique():
        mask = adata.obs["sample"] == sample
        mt = adata.obs.loc[mask, "pct_counts_mt"]

        median = np.median(mt)
        mad = np.median(np.abs(mt - median))

        threshold = median + nmads * mad
        outliers = (mt - median) > nmads * mad

        adata.obs.loc[mask, "mt_outlier"] = outliers

        n_cells = int(mask.sum())
        n_outliers = int(outliers.sum())
        frac = float(n_outliers / n_cells)

        mt_stats[sample] = {
            "n_cells_pre": n_cells,
            "n_outliers": n_outliers,
            "fraction_outliers": frac,
            "threshold": float(threshold),
            "median": float(median),
            "mad": float(mad),
            "nmads": nmads,
        }

        print(f"{sample}: {frac:.3f}")

    adata.uns["mt_qc"] = mt_stats

    return adata


def mt_qc_panel(adata, outfile, nmads):
    import numpy as np
    import matplotlib.pyplot as plt
    import seaborn as sns

    mt = adata.obs["pct_counts_mt"].values
    counts = adata.obs["total_counts"].values
    outliers = adata.obs["mt_outlier"].values

    median = np.median(mt)
    mad = np.median(np.abs(mt - median))
    threshold = median + nmads * mad
    fraction = np.mean(outliers) * 100

    fig = plt.figure(figsize=(16, 10))
    gs = fig.add_gridspec(2, 3)

    ax1 = fig.add_subplot(gs[0, 0])
    ax2 = fig.add_subplot(gs[0, 1])
    ax3 = fig.add_subplot(gs[0, 2])
    ax4 = fig.add_subplot(gs[1, :])

    # HEXBIN
    hb = ax1.hexbin(
        counts,
        mt,
        gridsize=150,
        bins="log",
        xscale="log",
        cmap="viridis",
        mincnt=1,
    )
    ax1.set_xlabel("Total counts (log)")
    ax1.set_ylabel("% MT")
    ax1.set_title("MT density")
    fig.colorbar(hb, ax=ax1, label="log10(cell density)")

    # OUTLIERS
    ax2.scatter(counts[~outliers], mt[~outliers], c="lightgray", s=5)
    ax2.scatter(counts[outliers], mt[outliers], c="red", s=5)
    ax2.axhline(threshold, color="black", linestyle="--")

    ax2.set_xscale("log")
    ax2.set_title(f"MT outliers ({fraction:.1f}%)")

    # HIST
    ax3.hist(mt, bins=100, log=True)
    ax3.axvline(threshold, color="red", linestyle="--")
    ax3.set_title("MT distribution")

    # VIOLIN
    sns.violinplot(
        data=adata.obs,
        x="sample",
        y="pct_counts_mt",
        ax=ax4,
        inner=None,
    )

    for label in ax4.get_xticklabels():
        label.set_rotation(90)

    ax4.set_title("MT per sample")

    plt.tight_layout()
    plt.savefig(outfile, dpi=300)
    plt.close()


def run_mt_qc(args):
    import scanpy as sc

    input_file = Path(args.input)
    output_file = input.replace(".h5ad", f"_mtqc_nMADs-{args.nmads}.h5ad")
    plot_file = input.replace(".h5ad", f"_mtqc_nMADs-{args.nmads}.png")

    print(f"Loading {input_file}")
    adata = sc.read(input_file)

    adata = compute_mt_outliers(adata, args.nmads)

    print("Generating MT QC plot...")
    mt_qc_panel(adata, plot_file, args.nmads)

    print(f"Saving updated AnnData → {output_file}")
    adata.write(output_file)
