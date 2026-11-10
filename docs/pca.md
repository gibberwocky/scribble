# pca

## Purpose

The `pca` module generates a diagnostic PCA comparison before and after QC filtering.

It is intended as an assessment tool to visualise the effects of mitochondrial filtering, doublet removal, and gene-count thresholds before proceeding to integration and clustering.

No cells are removed and no AnnData object is modified. The module only produces diagnostic figures.

## Workflow Position

```text
import
  ↓
mt
  ↓
dbl
  ↓
pca
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

- AnnData object processed by the `mt` and `dbl` modules
- Mitochondrial outlier calls (`mt_outlier`)
- Doublet classifications (user-specified column)
- Raw count matrix

## Outputs

- PCA QC comparison figure
- Module manifest

## Main Functions

- Perform PCA before filtering
- Perform PCA after applying QC thresholds
- Compare cell distributions before and after filtering
- Visualise counts, mitochondrial content, and doublet scores on PCA embeddings
- Report overall cell retention after QC filtering

## Parameters

### Required Parameters

| Parameter | Description |
|------------|-------------|
| `--project_dir` | Root Scribble project directory. |
| `--input` | Input AnnData object, typically the output from the `dbl` module. |
| `--doublets` | Name of the doublet classification column in `adata.obs` (e.g. `doublet_outlier`). |

### QC Thresholds

| Parameter | Default | Description |
|------------|----------|-------------|
| `--mingenes` | 200 | Minimum number of detected genes used when generating the post-filtering PCA. |
| `--maxgenes` | 6000 | Maximum number of detected genes used when generating the post-filtering PCA. |
| `--hvgs` | 3000 | Number of highly variable genes used for PCA calculation. |
| `--vmax` | 0.99 | Upper percentile used for colour scaling in PCA visualisations. |

### Reproducibility

| Parameter | Default | Description |
|------------|----------|-------------|
| `--manifest` | None | Load parameters from a previous Scribble manifest file. |

### Notes

- The specified doublet column must exist in `adata.obs`.
- PCA is recalculated both before and after applying QC filters for visual comparison.
- The module does not modify or save a new AnnData object.
- Generated PCA plots are intended as QC diagnostics rather than downstream analysis results.
