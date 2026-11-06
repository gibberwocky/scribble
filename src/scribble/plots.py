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
def doublet_qc_panel(np, plt, sns, adata, outfile, threshold=0.93):
    has_scrublet = ("scrublet_score" in adata.obs.columns)
    has_scdblfinder = ("scdblfinder_score" in adata.obs.columns)
    has_score = (has_scrublet or has_scdblfinder)

    scrublet_scores = None
    scdblfinder_scores = None

    if has_scrublet:
        scrublet_scores = (adata.obs["scrublet_score"].values)

    if has_scdblfinder:
        scdblfinder_scores = (adata.obs["scdblfinder_score"].values)

    counts = adata.obs["total_counts"].values
    outliers = adata.obs["doublet_outlier"].values
    frac_total = float(np.mean(outliers) * 100)
    # -------- compute per-sample thresholds --------
    sample_thresholds = {}

    if has_scrublet:
        for sample in adata.obs["sample"].unique():
            mask = (adata.obs["sample"] == sample)
            s = (adata.obs.loc[mask, "scrublet_score"].dropna())
            if len(s) > 0:
                sample_thresholds[sample] = (np.quantile(s, threshold))

    # -------- layout --------
    fig = plt.figure(figsize=(12,8))
    gs = fig.add_gridspec(2, 2)
    ax1 = fig.add_subplot(gs[0, 0])   # distribution
    ax2 = fig.add_subplot(gs[0, 1])   # scatter
    ax3 = fig.add_subplot(gs[1, :])   # wide violin
    # -------- 1. DISTRIBUTION --------
    if has_score:

        if has_scrublet:

            ax1.hist(
                scrublet_scores,
                bins=100,
                alpha=0.5,
                log=True,
                label="Scrublet",
            )

        if (has_scdblfinder and np.unique(scdblfinder_scores).size > 1):

            ax1.hist(
                scdblfinder_scores,
                bins=100,
                alpha=0.5,
                log=True,
                label="scDblFinder",
            )

        ax1.set_xlabel("Doublet score")
        ax1.set_ylabel("Number of cells (log scale)")
        ax1.set_title("Doublet score distribution")
        ax1.legend()

    else:

        call_counts = {
            "doublet": int(
                adata.obs[
                    "doublet_outlier"
                ].sum()
            ),
            "singlet": int(
                (~adata.obs[
                    "doublet_outlier"
                ]).sum()
            ),
        }
        ax1.bar(call_counts.keys(), call_counts.values(),)
        ax1.set_ylabel("Number of cells")
        ax1.set_title("Doublet calls")


    # -------- 2. OUTLIER SCATTER --------
    if has_score:

        if has_scrublet:

            ax2.scatter(
                counts,
                scrublet_scores,
                c=np.where(
                    outliers,
                    "red",
                    "lightgray"
                ),
                s=5,
                alpha=0.5,
                label="Scrublet",
            )

        if has_scdblfinder:

            ax2.scatter(
                counts,
                scdblfinder_scores,
                c=np.where(
                    outliers,
                    "darkblue",
                    "steelblue"
                ),
                s=5,
                alpha=0.5,
                label="scDblFinder",
            )

        ax2.set_ylabel("Doublet score")

    else:

        ax2.scatter(
            counts,
            np.zeros_like(counts),
            c=np.where(
                outliers,
                "red",
                "lightgray"
            ),
            s=5,
        )

        ax2.set_ylabel("Doublet call")

    ax2.set_xscale("log")
    ax2.set_xlabel("Total counts (log scale)")
    ax2.set_title(f"Cells flagged as doublets ({frac_total:.1f}% of cells)")

    handles, labels = (ax2.get_legend_handles_labels())
    if len(handles) > 0:
        ax2.legend(
            markerscale=3
        )

    # -------- 3. VIOLIN PLOT --------
    if has_score:

        if has_scrublet:

            sns.violinplot(
                data=adata.obs,
                x="sample",
                y="scrublet_score",
                ax=ax3,
                inner=None,
                color="lightgrey",
            )

        if has_scdblfinder:

            sns.stripplot(
                data=adata.obs,
                x="sample",
                y="scdblfinder_score",
                ax=ax3,
                color="darkblue",
                size=1,
                alpha=0.3,
            )

        ax3.set_ylabel("Doublet score")

    else:

        call_fraction = (
            adata.obs
            .groupby("sample")[
                "doublet_outlier"
            ]
            .mean()
            .reset_index()
        )

        sns.barplot(
            data=call_fraction,
            x="sample",
            y="doublet_outlier",
            ax=ax3,
        )

        ax3.set_ylabel("Doublet fraction")

    # overlay thresholds per sample
    if has_scrublet:

        sample_order = (
            adata.obs["sample"]
            .astype(str)
            .unique()
        )

        for i, sample in enumerate(sample_order):
            if sample in sample_thresholds:
                ax3.hlines(
                    sample_thresholds[sample],
                    i - 0.4,
                    i + 0.4,
                    colors="red",
                    linestyles="--",
                )

        ax3.set_title("Scrublet score per sample")

    else:

        ax3.set_title("Doublet fraction per sample")

    ax3.set_xlabel("Sample")

    # rotate labels
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
    doublets,
    min_genes=200,
    max_genes=6000,
    hvg=3000,
    vmax_pc=0.99,
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
        (~adata.obs[doublets]) &
        (adata.obs["n_genes_by_counts"] > min_genes) &
        (adata.obs["n_genes_by_counts"] < max_genes)
    ].copy()
    # log counts for filtered object
    adata_after.obs["log10_total_counts"] = np.log10(
        adata_after.obs["total_counts"] + 1
    )
    # Identify doublet score column
    score_column = None

    if doublets == "scrublet_doublet":
        score_column = "scrublet_score"

    elif doublets == "scdblfinder_doublet":
        score_column = "scdblfinder_score"

    elif "scrublet_score" in adata.obs.columns:
        score_column = "scrublet_score"

    elif "scdblfinder_score" in adata.obs.columns:
        score_column = "scdblfinder_score"
    # Process data
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
    vmin_dbl = None
    vmax_dbl = None

    if score_column is not None:

        vmin_dbl = min(
            adata.obs[score_column].min(),
            adata_after.obs[score_column].min()
        )

        vmax_dbl = max(
            np.quantile(
                adata.obs[score_column].dropna(),
                vmax_pc
            ),
            np.quantile(
                adata_after.obs[score_column].dropna(),
                vmax_pc
            )
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
    if score_column is not None:

        sc.pl.pca(
            adata_before,
            color=score_column,
            vmin=vmin_dbl,
            vmax=vmax_dbl,
            ax=axes[0, 1],
            show=False
        )

        axes[0, 1].set_title(
            f"Before ({score_column})"
        )

        sc.pl.pca(
            adata_after,
            color=score_column,
            vmin=vmin_dbl,
            vmax=vmax_dbl,
            ax=axes[1, 1],
            show=False
        )

        axes[1, 1].set_title(
            f"After ({score_column})"
        )

    else:

        sc.pl.pca(
            adata_before,
            color=doublets,
            ax=axes[0, 1],
            show=False
        )

        sc.pl.pca(
            adata_after,
            color=doublets,
            ax=axes[1, 1],
            show=False
        )
    # ---------- GLOBAL ANNOTATION ----------
    frac_kept = adata_after.n_obs / adata.n_obs * 100
    fig.suptitle(
        f"QC filtering retained {frac_kept:.1f}% of cells",
        fontsize=14
    )
    plt.tight_layout()
    plt.savefig(outfile, dpi=300, bbox_inches="tight")
    plt.close()
