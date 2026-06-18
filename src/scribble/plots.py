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


# ------------ Doublet Panel -------------
def doublet_qc_panel(np, plt, threshold=0.93, adata, outfile):
    scores = adata.obs["doublet_score"].values
    counts = adata.obs["total_counts"].values
    outliers = adata.obs["doublet_outlier"].values
    frac_total = float(np.mean(outliers) * 100)
    # -------- compute per-sample thresholds --------
    sample_thresholds = {}
    for sample in adata.obs["sample"].unique():
        mask = adata.obs["sample"] == sample
        s = adata.obs.loc[mask, "doublet_score"]
        sample_thresholds[sample] = np.quantile(s, threshold)
    # -------- layout --------
    fig = plt.figure(figsize=(12,8))
    gs = fig.add_gridspec(2, 2)
    ax1 = fig.add_subplot(gs[0, 0])   # distribution
    ax2 = fig.add_subplot(gs[0, 1])   # scatter
    ax3 = fig.add_subplot(gs[1, :])   # wide violin
    # -------- 1. DISTRIBUTION --------
    ax1.hist(scores, bins=100, color="gray", alpha=0.7, log=True)
    ax1.set_xlabel("Doublet score")
    ax1.set_ylabel("Number of cells (log scale)")
    ax1.set_title("Doublet score distribution")
    # -------- 2. OUTLIER SCATTER --------
    ax2.scatter(
        counts[~outliers],
        scores[~outliers],
        c="lightgray",
        s=5,
        label="Kept"
    )
    ax2.scatter(
        counts[outliers],
        scores[outliers],
        c="red",
        s=5,
        label="Removed"
    )
    ax2.set_xscale("log")
    ax2.set_xlabel("Total counts (log scale)")
    ax2.set_ylabel("Doublet score")
    ax2.set_title(f"Doublet outliers ({frac_total:.1f}% of cells)")
    ax2.legend(markerscale=3)
    # -------- 3. VIOLIN PLOT --------
    sns.violinplot(
        data=adata.obs,
        x="sample",
        y="doublet_score",
        ax=ax3,
        inner=None
    )
    # overlay thresholds per sample
    for i, sample in enumerate(adata.obs["sample"].cat.categories):
        if sample in sample_thresholds:
            ax3.hlines(
                sample_thresholds[sample],
                i - 0.4,
                i + 0.4,
                colors="red",
                linestyles="--"
            )
    ax3.set_title("Doublet score per sample")
    ax3.set_xlabel("Sample")
    ax3.set_ylabel("Doublet score")
    # rotate labels (important!)
    for label in ax3.get_xticklabels():
        label.set_rotation(90)
    plt.tight_layout()
    plt.savefig(outfile, dpi=300, bbox_inches="tight")
    plt.close()


# ---------- PCA before/after ------------
def pca_qc_panel(
    np, sc, plt,
    adata,
    outfile,
    min_genes=200,
    max_genes=6000,
    hvg=3000,
    vmax_pc=0.99
):
    # ---------- PREP GLOBAL METRICS ----------
    # log counts stored in obs (fix for Scanpy behaviour)
    adata.obs["log10_total_counts"] = np.log10(adata.obs["total_counts"] + 1)
    # ---------- BEFORE ----------
    adata_before = adata.copy()
    sc.pp.normalize_total(adata_before)
    sc.pp.log1p(adata_before)
    sc.pp.highly_variable_genes(adata_before, n_top_genes=hvg)
    adata_before = adata_before[:, adata_before.var.highly_variable].copy()
    sc.pp.scale(adata_before, max_value=10)
    sc.tl.pca(adata_before)
    # ---------- AFTER ----------
    adata_after = adata[
        (~adata.obs["mt_outlier"]) &
        (~adata.obs["doublet_outlier"]) &
        (adata.obs["n_genes_by_counts"] > min_genes) &
        (adata.obs["n_genes_by_counts"] < max_genes)
    ].copy()
    # log counts for filtered object
    adata_after.obs["log10_total_counts"] = np.log10(
        adata_after.obs["total_counts"] + 1
    )
    sc.pp.normalize_total(adata_after)
    sc.pp.log1p(adata_after)
    sc.pp.highly_variable_genes(adata_after, n_top_genes=hvg)
    adata_after = adata_after[:, adata_after.var.highly_variable].copy()
    sc.pp.scale(adata_after, max_value=10)
    sc.tl.pca(adata_after)
    # ---------- STANDARDISE COLOR SCALES ----------
    # counts (log-transformed)
    vmin_counts = min(
        adata.obs["log10_total_counts"].min(),
        adata_after.obs["log10_total_counts"].min()
    )
    vmax_counts = max(
        np.quantile(adata.obs["log10_total_counts"], vmax_pc),
        np.quantile(adata_after.obs["log10_total_counts"], vmax_pc)
    )
    # doublet score
    vmin_dbl = min(
        adata.obs["doublet_score"].min(),
        adata_after.obs["doublet_score"].min()
    )
    vmax_dbl = max(
        np.quantile(adata.obs["doublet_score"], vmax_pc),
        np.quantile(adata_after.obs["doublet_score"], vmax_pc)
    )
    # MT %
    vmin_mt = min(
        adata.obs["pct_counts_mt"].min(),
        adata_after.obs["pct_counts_mt"].min()
    )
    vmax_mt = max(
        np.quantile(adata.obs["pct_counts_mt"], vmax_pc),
        np.quantile(adata_after.obs["pct_counts_mt"], vmax_pc)
    )
    # ---------- PLOT ----------
    fig, axes = plt.subplots(2, 3, figsize=(15, 10))
    # ----- COUNTS (log-scaled) -----
    sc.pl.pca(
        adata_before,
        color="log10_total_counts",
        vmin=vmin_counts,
        vmax=vmax_counts,
        ax=axes[0, 0],
        show=False
    )
    axes[0, 0].set_title("Before (log10 counts)")
    sc.pl.pca(
        adata_after,
        color="log10_total_counts",
        vmin=vmin_counts,
        vmax=vmax_counts,
        ax=axes[1, 0],
        show=False
    )
    axes[1, 0].set_title("After (log10 counts)")
    # ----- MT -----
    sc.pl.pca(
        adata_before,
        color="pct_counts_mt",
        vmin=vmin_mt,
        vmax=vmax_mt,
        ax=axes[0, 2],
        show=False
    )
    axes[0, 2].set_title("Before (MT %)")
    sc.pl.pca(
        adata_after,
        color="pct_counts_mt",
        vmin=vmin_mt,
        vmax=vmax_mt,
        ax=axes[1, 2],
        show=False
    )
    axes[1, 2].set_title("After (MT %)")
    # ----- DOUBLETS -----
    sc.pl.pca(
        adata_before,
        color="doublet_score",
        vmin=vmin_dbl,
        vmax=vmax_dbl,
        ax=axes[0, 1],
        show=False
    )
    axes[0, 1].set_title("Before (doublet score)")
    sc.pl.pca(
        adata_after,
        color="doublet_score",
        vmin=vmin_dbl,
        vmax=vmax_dbl,
        ax=axes[1, 1],
        show=False
    )
    axes[1, 1].set_title("After (doublet score)")
    # ---------- GLOBAL ANNOTATION ----------
    frac_kept = adata_after.n_obs / adata.n_obs * 100
    fig.suptitle(
        f"QC filtering retained {frac_kept:.1f}% of cells",
        fontsize=14
    )
    plt.tight_layout()
    plt.savefig(outfile, dpi=300, bbox_inches="tight")
    plt.close()
