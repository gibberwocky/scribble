#!/usr/bin/env python

from pathlib import Path

def run_map(args):
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

    # ----------------------------
    # Harmonise gene space
    # ----------------------------
    print("Harmonising gene space...")

    common = adata_ref.var_names.intersection(adata_query.var_names)
    adata_ref = adata_ref[:, common].copy()
    adata_query = adata_query[:, common].copy()

    # ----------------------------
    # Combine datasets
    # ----------------------------
    print("Combining datasets...")

    adata_combined = adata_ref.concatenate(
        adata_query,
        batch_key="dataset",
        batch_categories=["ref", "query"]
    )

    # ----------------------------
    # SCVI model
    # ----------------------------
    print("Training SCVI model...")

    scvi.model.SCVI.setup_anndata(
        adata_combined,
        batch_key="dataset"
    )

    model = scvi.model.SCVI(adata_combined)
    model.train()

    # ----------------------------
    # SCANVI model
    # ----------------------------
    print("Training SCANVI model...")

    scvi.model.SCANVI.setup_anndata(
        adata_combined,
        labels_key=args.label_key
    )

    scanvi = scvi.model.SCANVI.from_scvi_model(model)
    scanvi.train()

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
