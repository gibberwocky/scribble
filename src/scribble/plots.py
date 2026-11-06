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
def doublet_qc_panel(
    np,
    pd,
    plt,
    sns,
    adata,
    outfile,
    threshold=0.93,
):
    has_scrublet = ("scrublet_score" in adata.obs.columns)
    has_scdblfinder = ("scdblfinder_score" in adata.obs.columns)
    counts = adata.obs["total_counts"].values
    outliers = adata.obs["doublet_outlier"].values
    frac_total = float(np.mean(outliers) * 100)

    # --------------------------------------------------
    # Scrublet thresholds
    # --------------------------------------------------

    sample_thresholds = {}

    if has_scrublet:

        for sample in (
            adata.obs["sample"]
            .astype(str)
            .unique()
        ):

            mask = (
                adata.obs["sample"]
                .astype(str)
                == sample
            )

            s = (
                adata.obs.loc[
                    mask,
                    "scrublet_score"
                ]
                .dropna()
            )

            if len(s) > 0:

                sample_thresholds[
                    sample
                ] = np.quantile(
                    s,
                    threshold
                )

    # --------------------------------------------------
    # Figure layout
    # --------------------------------------------------

    fig = plt.figure(
        figsize=(16, 12)
    )

    gs = fig.add_gridspec(
        3,
        2,
        height_ratios=[1, 1, 1.2]
    )

    ax1 = fig.add_subplot(gs[0, 0])  # scrublet hist
    ax2 = fig.add_subplot(gs[0, 1])  # scrublet scatter

    ax3 = fig.add_subplot(gs[1, 0])  # scdblfinder hist
    ax4 = fig.add_subplot(gs[1, 1])  # scdblfinder scatter

    ax5 = fig.add_subplot(gs[2, :])  # sample summary

    # --------------------------------------------------
    # SCRUBLET HISTOGRAM
    # --------------------------------------------------

    if has_scrublet:

        scrublet_scores = (
            adata.obs[
                "scrublet_score"
            ]
            .dropna()
            .values
        )

        ax1.hist(
            scrublet_scores,
            bins=100,
            color="steelblue",
            alpha=0.8,
            log=True,
        )

        global_threshold = np.quantile(
            scrublet_scores,
            threshold
        )

        ax1.axvline(
            global_threshold,
            color="red",
            linestyle="--",
            linewidth=2,
        )

        ax1.set_title(
            "Scrublet score distribution"
        )

        ax1.set_xlabel(
            "Scrublet score"
        )

        ax1.set_ylabel(
            "Cells (log scale)"
        )

    else:

        ax1.axis("off")

    # --------------------------------------------------
    # SCRUBLET SCATTER
    # --------------------------------------------------

    if has_scrublet:

        scrublet_scores = (
            adata.obs[
                "scrublet_score"
            ]
            .values
        )

        ax2.scatter(
            counts[~outliers],
            scrublet_scores[~outliers],
            c="lightgrey",
            s=4,
            alpha=0.4,
        )

        ax2.scatter(
            counts[outliers],
            scrublet_scores[outliers],
            c="red",
            s=4,
            alpha=0.4,
        )

        ax2.set_xscale("log")

        ax2.set_title(
            "Scrublet score vs counts"
        )

        ax2.set_xlabel(
            "Total counts"
        )

        ax2.set_ylabel(
            "Scrublet score"
        )

    else:

        ax2.axis("off")

    # --------------------------------------------------
    # SCDBLFINDER HISTOGRAM
    # --------------------------------------------------

    if has_scdblfinder:

        scdbl_scores = (
            adata.obs[
                "scdblfinder_score"
            ]
            .dropna()
            .values
        )

        ax3.hist(
            scdbl_scores,
            bins=100,
            color="darkorange",
            alpha=0.8,
            log=True,
        )

        ax3.set_title(
            "scDblFinder score distribution"
        )

        ax3.set_xlabel(
            "scDblFinder score"
        )

        ax3.set_ylabel(
            "Cells (log scale)"
        )

    else:

        ax3.axis("off")

    # --------------------------------------------------
    # SCDBLFINDER SCATTER
    # --------------------------------------------------

    if has_scdblfinder:

        scdbl_scores = (
            adata.obs[
                "scdblfinder_score"
            ]
            .values
        )

        ax4.scatter(
            counts[~outliers],
            scdbl_scores[~outliers],
            c="lightgrey",
            s=4,
            alpha=0.4,
        )

        ax4.scatter(
            counts[outliers],
            scdbl_scores[outliers],
            c="red",
            s=4,
            alpha=0.4,
        )

        ax4.set_xscale("log")

        ax4.set_title(
            "scDblFinder score vs counts"
        )

        ax4.set_xlabel(
            "Total counts"
        )

        ax4.set_ylabel(
            "scDblFinder score"
        )

    else:

        ax4.axis("off")

    # --------------------------------------------------
    # SAMPLE-LEVEL SUMMARY
    # --------------------------------------------------

    rows = []

    for sample in (
        adata.obs["sample"]
        .astype(str)
        .unique()
    ):

        mask = (
            adata.obs["sample"]
            .astype(str)
            == sample
        )

        row = {
            "sample": sample,
            "Union": (
                adata.obs.loc[
                    mask,
                    "doublet_outlier"
                ]
                .mean()
            )
        }

        if has_scrublet:

            row["Scrublet"] = (
                adata.obs.loc[
                    mask,
                    "scrublet_doublet"
                ]
                .mean()
            )

        if has_scdblfinder:

            row["scDblFinder"] = (
                adata.obs.loc[
                    mask,
                    "scdblfinder_doublet"
                ]
                .mean()
            )

        rows.append(row)

    frac_df = pd.DataFrame(rows)

    plot_df = frac_df.melt(
        id_vars="sample",
        var_name="method",
        value_name="fraction"
    )

    sns.barplot(
        data=plot_df,
        x="sample",
        y="fraction",
        hue="method",
        ax=ax5,
    )

    ax5.set_ylabel(
        "Fraction flagged"
    )

    ax5.set_xlabel(
        "Sample"
    )

    ax5.set_title(
        (
            f"Doublet detection summary "
            f"(union={frac_total:.1f}% cells)"
        )
    )

    for label in ax5.get_xticklabels():

        label.set_rotation(90)

    plt.tight_layout()

    plt.savefig(
        outfile,
        dpi=300,
        bbox_inches="tight",
    )

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
