#!/usr/bin/env python

from pathlib import Path

def qc_scvi_input(adata, name):

    X = adata.X.data if sparse.issparse(adata.X) else np.asarray(adata.X)
    cell_sums = np.asarray(
        adata.X.sum(axis=1)
    ).ravel()

    print(
        name,
        np.nanmin(X),
        np.nanmax(X),
        np.isnan(X).sum(),
        (X < 0).sum(),
        "zero_cells=",
        (cell_sums == 0).sum()
    )

def run_map(args):
    import anndata as ad
    import scanpy as sc
    import numpy as np
    import pandas as pd
    import matplotlib.pyplot as plt
    from scipy.stats import entropy
    from scipy import sparse
    import random
    import scvi

    from scribble.import_data import setup_environment

    # ----------------------------
    # Setup
    # ----------------------------
    PROJECT_DIR = Path(args.project_dir)
    PLOT_DIR = PROJECT_DIR / "scribble/plots"
    TABLE_DIR = PROJECT_DIR / "scribble/tables"
    DATA_DIR = PROJECT_DIR / "scribble/data"

    setup_environment(sc, np, random, PLOT_DIR)

    # ----------------------------
    # Load data
    # ----------------------------
    print("Loading datasets...")

    adata_ref = sc.read(DATA_DIR / args.reference)
    adata_query = sc.read(DATA_DIR / args.query)

    # Set counts
    if "counts" in adata_ref.layers:
        adata_ref.X = adata_ref.layers["counts"].copy()

    if "counts" in adata_query.layers:
        adata_query.X = adata_query.layers["counts"].copy()

    # Set query annotations as "Unknown"
    adata_query.obs[args.label_key] = "Unknown"

    # Check counts
    qc_scvi_input(adata_ref, "ref")
    qc_scvi_input(adata_query, "query")

    # ----------------------------
    # Harmonise gene space
    # ----------------------------
    print("Harmonising gene space...")
    common = adata_ref.var_names.intersection(adata_query.var_names)

    print(f"Reference genes: {adata_ref.n_vars}")
    print(f"Query genes: {adata_query.n_vars}")

    adata_ref = adata_ref[:, common].copy()
    adata_query = adata_query[:, common].copy()
    print(f"Shared genes: {len(common)}")

    # ----------------------------
    # Combine datasets
    # ----------------------------
    print("Combining datasets...")

    adata_combined = ad.concat(
        [adata_ref, adata_query],
        label="dataset",
        keys=["ref", "query"],
        join="inner"
    )
    print(adata_combined.obs["dataset"].value_counts())

    # ----------------------------
    # SCVI model
    # ----------------------------
    print("Training SCVI model...")

    scvi.model.SCVI.setup_anndata(
        adata_combined,
        batch_key="dataset"
    )

    model = scvi.model.SCVI(adata_combined)
    model.train(accelerator="auto")

    # ----------------------------
    # SCANVI model
    # ----------------------------
    print("Training SCANVI model...")

    print(
        adata_combined.obs[args.label_key]
        .value_counts(dropna=False)
    )

    scanvi = scvi.model.SCANVI.from_scvi_model(
        model,
        labels_key=args.label_key,
        unlabeled_category="Unknown"
    )

    scanvi.train(
        accelerator="auto"
    )

    # ----------------------------
    # Predictions
    # ----------------------------
    print("Predicting labels...")

    preds = scanvi.predict(adata_combined)
    probs = scanvi.predict_proba(adata_combined)

    adata_combined.obs["predicted_labels"] = preds
    adata_combined.obsm["prediction_probs"] = probs

    # ----------------------------
    # Split query
    # ----------------------------
    adata_query = adata_combined[
        adata_combined.obs["dataset"] == "query"
    ].copy()

    # ----------------------------
    # Confidence metrics
    # ----------------------------
    print("Calculating confidence scores...")

    max_prob = np.max(adata_query.obsm["prediction_probs"], axis=1)
    ent = entropy(adata_query.obsm["prediction_probs"].T)

    # Normalised entropy (0 = confident, 1 = uncertain)
    ent_norm = ent / np.log(adata_query.obsm["prediction_probs"].shape[1])

    adata_query.obs["prediction_confidence"] = max_prob
    adata_query.obs["prediction_entropy"] = ent_norm

    # ----------------------------
    # Apply confidence filtering
    # ----------------------------
    print("Applying confidence filtering...")

    adata_query.obs["cell_type_major_pred"] = adata_query.obs["predicted_labels"]

    low_conf = adata_query.obs["prediction_confidence"] < args.confidence_threshold

    adata_query.obs.loc[low_conf, "cell_type_major_pred"] = "Unassigned"

    # ----------------------------
    # Label smoothing
    # ----------------------------
    print("Smoothing labels...")

    sc.pp.neighbors(adata_query, n_neighbors=args.neighbors)

    conn = adata_query.obsp["connectivities"]

    labels = adata_query.obs["cell_type_major_pred"].astype(str)

    smoothed = []

    for i in range(conn.shape[0]):
        idx = conn[i].nonzero()[1]
        vote = labels.iloc[idx]
        smoothed.append(vote.value_counts().idxmax())

    adata_query.obs["cell_type_major_smooth"] = smoothed

    adata_query.obs["cell_type_major_final"] = adata_query.obs["cell_type_major_smooth"]

    # ----------------------------
    # Map lineage
    # ----------------------------
    print("Mapping lineage...")

    lineage_map = (
        adata_ref.obs[[args.label_key, args.lineage_key]]
        .drop_duplicates()
        .set_index(args.label_key)[args.lineage_key]
        .to_dict()
    )

    adata_query.obs["lineage"] = adata_query.obs["cell_type_major_final"].map(lineage_map)

    # ----------------------------
    # UMAP for query (for diagnostics)
    # ----------------------------
    print("Computing UMAP for query...")

    sc.tl.umap(adata_query)

    # ----------------------------
    # Diagnostic plots
    # ----------------------------
    print("Generating plots...")

    # Predicted labels
    sc.pl.umap(
        adata_query,
        color="cell_type_major_final",
        title="Predicted cell types",
        show=False
    )
    plt.savefig(PLOT_DIR / "map_predicted_labels_umap.png", dpi=300)
    plt.close()

    # Confidence
    sc.pl.umap(
        adata_query,
        color="prediction_confidence",
        cmap="viridis",
        title="Prediction confidence",
        show=False
    )
    plt.savefig(PLOT_DIR / "map_confidence_umap.png", dpi=300)
    plt.close()

    # Entropy
    sc.pl.umap(
        adata_query,
        color="prediction_entropy",
        cmap="magma",
        title="Prediction uncertainty",
        show=False
    )
    plt.savefig(PLOT_DIR / "map_entropy_umap.png", dpi=300)
    plt.close()

    # Raw predictions
    sc.pl.umap(
        adata_query,
        color="cell_type_major_pred",
        title="Raw predictions",
        show=False
    )
    plt.savefig(PLOT_DIR / "map_raw_predictions_umap.png", dpi=300)
    plt.close()

    # ----------------------------
    # Save results
    # ----------------------------
    print("Saving results...")

    adata_query.write(DATA_DIR / "adata_query_mapped.h5ad")

    # Save summary table
    summary = (
        adata_query.obs["cell_type_major_final"]
        .value_counts(normalize=True)
        .rename("fraction")
    )
    summary.to_csv(TABLE_DIR / "map_celltype_composition.csv")

    print("Mapping complete ✅")
