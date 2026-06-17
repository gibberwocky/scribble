#!/usr/bin/env python

# ---------- Hexbin density plot ---------
def qc_hexbin_ax(
    ax,
    x,
    y,
    xlabel,
    ylabel,
    log_x=True,
    log_y=False,
    gridsize=175,
    cmap="viridis"
):
    hb = ax.hexbin(
        x,
        y,
        gridsize=gridsize,
        bins="log",
        xscale="log" if log_x else "linear",
        yscale="log" if log_y else "linear",
        cmap=cmap,
        mincnt=1
    )
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    return hb

# ------------ Knee plot panel -----------
def knee_plot_ax(ax, df, inflection):
    # Rank cutoff
    rank_cutoff = df.loc[df["total"] > inflection, "rank"].max()
    # Curve
    ax.plot(df["total"], df["rank"], color="black", linewidth=1)
    # Threshold lines
    ax.axvline(inflection, linestyle="--", color="gray")
    ax.axhline(rank_cutoff, linestyle="--", color="gray")
    # Annotation
    ax.text(
        inflection,
        rank_cutoff,
        f"{int(rank_cutoff):,} cells",
        verticalalignment="bottom"
    )
    # Log scaling
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlabel("Total UMIs (log scale)")
    ax.set_ylabel("Rank")
    ax.set_title("Knee plot")

# -------------- MT QC Panel -------------
def mt_qc_panel(np, plt, adata, outfile, nmads):
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
