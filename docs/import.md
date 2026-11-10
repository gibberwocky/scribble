# import

## Purpose

The `import` module creates the initial AnnData object for a Scribble analysis.

It imports Cell Ranger or CellBender count matrices together with Velocyto loom files, standardises cell barcodes, calculates basic QC metrics, performs inflection-point based cell calling, attaches sample metadata, and generates sample-level QC plots.

The resulting dataset becomes the starting point for all downstream analyses.

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

- Cell Ranger matrices or CellBender filtered matrices
- Velocyto loom files
- Metadata Excel workbook

## Outputs

- `scribble/adata/combined.h5ad`
- Per-sample QC plots
- Import manifest

## Main Functions

- Load count matrices
- Load velocity information
- Harmonise cell barcodes
- Calculate QC metrics
- Estimate cell-calling threshold
- Filter low-count droplets
- Attach metadata
- Generate QC visualisations
- Create a combined AnnData object

## Parameters

### Required Parameters

| Parameter | Description |
|------------|-------------|
| `--project_dir` | Root Scribble project directory. |
| `--metadata_file` | Metadata workbook containing the `meta` sheet and `sample` column. |

### Count Matrix Sources

| Parameter | Description |
|------------|-------------|
| `--cellranger_dir` | Directory containing Cell Ranger outputs. |
| `--cellbender_dir` | Directory containing CellBender filtered matrices. Preferred over Cell Ranger if available. |
| `--velocyto_dir` | Directory containing Velocyto loom files. |

### Sample Selection

| Parameter | Default | Description |
|------------|----------|-------------|
| `--samples` | All samples | Import only the specified samples. |
| `--samples_file` | None | Text file containing one sample name per line. |

### Cell Calling

| Parameter | Default | Description |
|------------|----------|-------------|
| `--inflection_lower` | 100 | Minimum UMI count considered when estimating the inflection point for cell calling. |

### Reproducibility

| Parameter | Default | Description |
|------------|----------|-------------|
| `--manifest` | None | Load parameters from a previous Scribble manifest file. |

### Notes

- CellBender matrices are used in preference to Cell Ranger matrices when both are present.
- Cells below the inferred inflection-point threshold are removed during import.
- All supplied samples must have matching count matrices, loom files, and metadata entries.
