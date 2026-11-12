# preintegration

## Purpose

The `preintegration` module performs preprocessing of the filtered dataset prior to batch correction and clustering.

It selects highly variable genes, normalises expression values, optionally regresses out technical covariates, performs PCA, constructs a neighbour graph, generates an initial UMAP, and prepares the dataset for Harmony integration.

The output is a reduced-dimension AnnData object containing log-normalised highly variable genes together with the full log-normalised gene space stored in `adata.raw`.

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

- Filtered AnnData object from the `filter` module
- Sample metadata
- QC metrics
- Raw count matrix

## Outputs

- `<input>_preintegration.h5ad`
- PCA coloured by metadata variables
- PCA coloured by QC metrics
- Initial UMAP coloured by metadata variables
- Module manifest

## Main Functions

- Filter genes by minimum detection frequency
- Select highly variable genes (HVGs)
- Normalise counts and log transform expression
- Store the full gene-space matrix in `adata.raw`
- Restrict downstream analysis to HVGs
- Optionally regress out technical variables
- Optionally scale expression values
- Perform PCA
- Construct a neighbour graph
- Generate an initial UMAP
- Create a preintegration AnnData object for batch correction

## Parameters

### Required Parameters

| Parameter | Description |
|------------|-------------|
| `--project_dir` | Root Scribble project directory. |
| `--input` | Filtered AnnData object, typically the output from the `filter` module. |

### Highly Variable Gene Selection

| Parameter | Default | Description |
|------------|----------|-------------|
| `--hvgs` | 3000 | Number of highly variable genes to retain using Seurat v3 HVG selection. |
| `--batch` | sample | Metadata column used during batch-aware HVG selection. |
| `--min_cells_per_gene` | 3 | Minimum number of cells in which a gene must be detected to be retained. |

### Dimensionality Reduction

| Parameter | Default | Description |
|------------|----------|-------------|
| `--npcs` | 30 | Number of principal components to compute and use for neighbour graph construction. |
| `--neighbors` | 15 | Number of neighbours used when constructing the KNN graph. |

### Visualisation

| Parameter | Description |
|------------|-------------|
| `--vars` | Metadata columns to display on PCA and UMAP plots. |

### Regression

| Parameter | Default | Description |
|------------|----------|-------------|
| `--regress` | None | Metadata variables to regress out prior to PCA. |

### Scaling

| Parameter | Default | Description |
|------------|----------|-------------|
| `--no_scale` | False | Skip gene scaling before PCA. |

### Reproducibility

| Parameter | Default | Description |
|------------|----------|-------------|
| `--manifest` | None | Load parameters from a previous Scribble manifest file. |

### Notes

- HVGs are selected using the Seurat v3 method.
- Counts are normalised and log-transformed before PCA.
- The full log-normalised gene space is preserved in `adata.raw`.
- PCA, neighbour graph construction, and UMAP are performed using only the selected HVGs.
- If regression variables are specified, they must exist in `adata.obs`.
- Generated PCA and UMAP plots are intended for pre-integration assessment and batch-effect inspection.
