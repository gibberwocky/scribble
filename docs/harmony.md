# harmony

## Purpose

The `harmony` module performs batch correction using Harmony and generates an integrated low-dimensional representation of the dataset.

Harmony is applied to the PCA embedding generated during the preintegration step. Batch-corrected PCs are then used to construct a neighbour graph and generate an integrated UMAP.

The resulting AnnData object is used for downstream clustering, cluster evaluation, and annotation.

## Workflow Position

```text
import
  ↓
mt
  ↓
dbl
  ↓
filter
  ↓
preintegration
  ↓
harmony
  ↓
cluster
```

## Inputs

- Preintegration AnnData object
- PCA embedding (`X_pca`)
- Batch annotation column
- Sample and metadata annotations

## Outputs

- Harmony-corrected AnnData object
- Integrated UMAP coloured by metadata variables
- Module manifest

## Main Functions

- Run Harmony batch correction on PCA embeddings
- Create a Harmony-corrected PCA space
- Construct a neighbour graph using Harmony embeddings
- Generate an integrated UMAP
- Save batch-corrected embeddings for downstream clustering
- Generate UMAP visualisations for batch assessment

## Parameters

### Required Parameters

| Parameter | Description |
|------------|-------------|
| `--project_dir` | Root Scribble project directory. |
| `--input` | Preintegration AnnData object containing PCA embeddings (`X_pca`). |
| `--batch` | Metadata column used for Harmony batch correction. |

### Harmony Configuration

| Parameter | Default | Description |
|------------|----------|-------------|
| `--theta` | 2 | Harmony diversity penalty controlling the strength of batch correction. Larger values increase batch mixing. |

### Neighbour Graph Construction

| Parameter | Default | Description |
|------------|----------|-------------|
| `--npcs` | 50 | Number of Harmony PCs used when constructing the neighbour graph. |
| `--neighbors` | 15 | Number of neighbours used for KNN graph construction. |

### Visualisation

| Parameter | Description |
|------------|-------------|
| `--vars` | Metadata columns to display on Harmony-corrected UMAPs. |

### Reproducibility

| Parameter | Default | Description |
|------------|----------|-------------|
| `--manifest` | None | Load parameters from a previous Scribble manifest file. |

### Notes

- The input AnnData object must contain `X_pca`. Run `preintegration` first.
- Harmony is applied to PCA coordinates, not directly to the expression matrix.
- Batch-corrected embeddings are stored in `adata.obsm["X_pca_harmony"]`.
- Neighbour graph construction and UMAP generation use the Harmony embedding.
- The resulting Harmony UMAP is intended for assessing batch mixing and cluster structure prior to clustering.
