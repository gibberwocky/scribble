#!/usr/bin/env python

from pathlib import Path

def qc_scvi_input(adata, name):
    import numpy as np
    from scipy import sparse

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

def run_scanvi_mapping(
    adata_ref,
    adata_query,
    label_key,
    confidence_threshold,
    neighbors,
):
    """
    Map a query dataset to a reference using a single label_key.
    Returns dataframe of predictions indexed by query cell names.
    """

    import anndata as ad
    import scvi
    import scanpy as sc
    import numpy as np
    import pandas as pd
    from scipy.stats import entropy

    # preserve existing annotation
    adata_query_map = adata_query.copy()
    adata_query_map.obs[label_key] = "Unknown"

    adata_combined = ad.concat(
        [adata_ref, adata_query_map],
        label="dataset",
        keys=["ref", "query"],
        join="inner"
    )

    print(f"\nMapping: {label_key}")

    scvi.model.SCVI.setup_anndata(
        adata_combined,
        batch_key="dataset"
    )

    model = scvi.model.SCVI(adata_combined)
    model.train(accelerator="auto")

    scanvi = scvi.model.SCANVI.from_scvi_model(
        model,
        labels_key=label_key,
        unlabeled_category="Unknown"
    )

    scanvi.train(accelerator="auto")

    preds = scanvi.predict(adata_combined)
    probs = scanvi.predict(
        adata_combined,
        soft=True
    )

    if hasattr(probs, "values"):
        probs = probs.values

    probs = np.asarray(probs)


    adata_combined.obs["pred"] = preds

    query_mask = (
        adata_combined.obs["dataset"] == "query"
    ).values

    mapped = adata_combined[
        query_mask
    ].copy()

    query_probs = probs[query_mask]
    mapped.obsm["probs"] = query_probs

    max_prob = np.max(query_probs, axis=1)
    ent = entropy(query_probs.T)
    ent = ent / np.log(query_probs.shape[1])

    mapped.obs["confidence"] = max_prob
    mapped.obs["entropy"] = ent
    mapped.obs["pred_filtered"] = mapped.obs["pred"]

    low_conf = (mapped.obs["confidence"] < confidence_threshold)
    mapped.obs.loc[low_conf,"pred_filtered"] = "Unassigned"

    sc.pp.neighbors(mapped,n_neighbors=neighbors)

    conn = mapped.obsp["connectivities"]

    labels = mapped.obs["pred_filtered"].astype(str)

    smoothed = []

    for i in range(conn.shape[0]):
        idx = conn[i].nonzero()[1]
        vote = labels.iloc[idx]
        smoothed.append(
            vote.value_counts().idxmax()
        )

    mapped.obs["pred_final"] = smoothed

    return pd.DataFrame(
        {
            f"{label_key}_reference_pred":
                mapped.obs["pred"],
            f"{label_key}_reference_final":
                mapped.obs["pred_final"],
            f"{label_key}_reference_confidence":
                mapped.obs["confidence"],
            f"{label_key}_reference_entropy":
                mapped.obs["entropy"],
        },
        index=mapped.obs_names
    )

def run_map(args):
    import random
    import scvi
    import anndata as ad
    import scanpy as sc
    import numpy as np
    import pandas as pd
    import matplotlib.pyplot as plt
    from scipy.stats import entropy
    from scipy import sparse

    from scribble.refine import restore_counts
    from scribble.import_data import setup_environment

    # ----------------------------
    # Setup
    # ----------------------------
    PROJECT_DIR = Path(args.project_dir)
    PLOT_DIR = PROJECT_DIR / "scribble/plots"
    TABLE_DIR = PROJECT_DIR / "scribble/tables"
    DATA_DIR = PROJECT_DIR / "scribble/adata"

    PLOT_DIR.mkdir(exist_ok=True, parents=True)
    DATA_DIR.mkdir(exist_ok=True, parents=True)
    TABLE_DIR.mkdir(exist_ok=True, parents=True)

    print("PLOT_DIR :", PLOT_DIR)
    print("TABLE_DIR:", TABLE_DIR)
    print("DATA_DIR:", DATA_DIR)

    print("PLOT exists :", PLOT_DIR.exists())
    print("TABLE exists:", TABLE_DIR.exists())
    print("DATA exists:", DATA_DIR.exists())

    setup_environment(sc, np, random, PLOT_DIR)

    # ----------------------------
    # Load data
    # ----------------------------
    print("Loading datasets...")

    adata_ref = sc.read(args.reference)
    adata_query_original = sc.read(args.query)
    adata_query = adata_query_original.copy()

    # Restore raw counts if available
    try:
        adata_tmp = restore_counts(adata_ref)
        # preserve annotations/embeddings
        adata_tmp.obsm = adata_ref.obsm.copy()
        adata_tmp.obs = adata_ref.obs.copy()
        adata_tmp.uns = adata_ref.uns.copy()
        adata_ref = adata_tmp
        adata_ref.uns.pop("log1p", None)
        print("Reference counts restored")
    except Exception as e:
        print(f"Could not restore reference counts: {e}")

    try:
        adata_tmp = restore_counts(adata_query)
        # preserve annotations/embeddings
        adata_tmp.obsm = adata_query.obsm.copy()
        adata_tmp.obs = adata_query.obs.copy()
        adata_tmp.uns = adata_query.uns.copy()
        adata_query = adata_tmp
        adata_query.uns.pop("log1p", None)
        print("Query counts restored")
    except Exception as e:
        print(f"Could not restore query counts: {e}")

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
    # Major cell type mapping
    # ----------------------------

    major_results = run_scanvi_mapping(
        adata_ref=adata_ref,
        adata_query=adata_query,
        label_key="cell_type_major",
        confidence_threshold=args.confidence_threshold,
        neighbors=args.neighbors,
    )

    for col in major_results.columns:
        adata_query_original.obs.loc[
            major_results.index,
            col
        ] = major_results[col]


    if "cell_type_minor" in adata_ref.obs.columns:

        minor_results = run_scanvi_mapping(
            adata_ref=adata_ref,
            adata_query=adata_query,
            label_key="cell_type_minor",
            confidence_threshold=args.confidence_threshold,
            neighbors=args.neighbors,
        )

        for col in minor_results.columns:
            adata_query_original.obs.loc[
                minor_results.index,
                col
            ] = minor_results[col]

    # ----------------------------
    # UMAP for query (for diagnostics)
    # ----------------------------
    print("Computing UMAP for query...")

    if "X_umap" not in adata_query_original.obsm:
        sc.pp.neighbors(
            adata_query_original,
            n_neighbors=args.neighbors
        )
        sc.tl.umap(adata_query_original)

    # ----------------------------
    # Diagnostic plots
    # ----------------------------
    print("Generating plots...")

    sc.pl.umap(
        adata_query_original,
        color="cell_type_major",
        title="Original annotation",
        show=False
    )
    plt.savefig(
        PLOT_DIR / "map_original_annotation_umap.png",
        dpi=300
    )
    plt.savefig(PLOT_DIR / "umap_original_query_annotation.png", dpi=300)
    plt.close()

    sc.pl.umap(
        adata_query_original,
        color="cell_type_major_reference_final",
        title="Atlas prediction",
        show=False
    )
    plt.savefig(
        PLOT_DIR / "umap_reference_prediction_annotation.png",
        dpi=300
    )
    plt.close()

    # ----------------------------
    # Save results
    # ----------------------------
    print("Saving results...")
    adata_query_original.write(DATA_DIR / "adata_query_mapped.h5ad")



    print("Mapping complete ✅")
