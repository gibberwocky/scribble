# annotate

## Purpose

The `annotate` module applies manually curated biological annotations to clusters and generates publication-ready visualisations.

Annotations are supplied through an Excel workbook and mapped onto cluster labels. The module can generate annotated UMAPs, marker gene dotplots, and marker expression UMAPs while preserving cluster-to-annotation mappings inside the AnnData object.

This step produces the final biologically annotated dataset used for interpretation and downstream analysis.

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
  ↓
evaluate
  ↓
refine
  ↓
curate
  ↓
annotate
```

## Inputs

- Curated or refined AnnData object
- Annotation workbook containing cluster annotations
- Cluster labels to annotate
- Optional marker gene definitions
- Optional colour palette definitions

## Outputs

- Annotated AnnData object
- Annotated UMAP
- Optional marker dotplots
- Optional marker expression UMAPs
- Annotation manifest

## Main Functions

- Map biological annotations onto clusters
- Create sequential cluster IDs
- Add annotation metadata to `adata.obs`
- Generate annotated UMAP visualisations
- Apply custom colour palettes
- Create marker gene dotplots
- Generate marker expression UMAPs
- Store cluster-to-annotation mappings
- Export an annotated AnnData object

## Parameters

### Required Parameters

| Parameter | Description |
|------------|-------------|
| `--project_dir` | Root Scribble project directory. |
| `--input` | Curated or refined AnnData object to annotate. |
| `--annotations` | Excel workbook containing cluster annotations. |
| `--cluster_field` | Cluster column in `adata.obs` used for mapping annotations (e.g. `refine_cluster`). |
| `--label` | Annotation column to use when generating annotated UMAPs. |

### Cluster Mapping

| Parameter | Default | Description |
|------------|----------|-------------|
| `--cluster_destination` | curated_cluster | Name of the sequential cluster ID field created in `adata.obs`. |

### UMAP Visualisation

| Parameter | Default | Description |
|------------|----------|-------------|
| `--legend_loc` | right margin | Location of the UMAP legend. |
| `--legend_fontsize` | 6 | Legend text size. |
| `--legend_fontoutline` | 0 | Outline width applied to legend text. |
| `--add_outline` | False | Add outlines around plotted cell groups. |
| `--plot_prefix` | None | Prefix applied to generated plot filenames. |
| `--plot_samples` | All samples | Restrict visualisation to a subset of samples. |
| `--combined_label` | False | Prefix annotation labels with cluster IDs in UMAP legends. |

### Annotation Colours

| Parameter | Default | Description |
|------------|----------|-------------|
| `--cluster_hex` | None | Column in the annotation workbook containing hexadecimal colour codes. Used to create custom UMAP palettes. |

### Marker Visualisation

| Parameter | Default | Description |
|------------|----------|-------------|
| `--plot_markers` | None | Column in the annotation workbook containing marker genes for each annotation. Enables marker dotplots and marker UMAP generation. |

### Output Behaviour

| Parameter | Default | Description |
|------------|----------|-------------|
| `--skip_write` | False | Generate plots and annotation metadata without writing a new AnnData file. |

### Reproducibility

| Parameter | Default | Description |
|------------|----------|-------------|
| `--manifest` | None | Load parameters from a previous Scribble manifest file. |

### Annotation Workbook Requirements

The workbook must contain a worksheet named:

```text
annotations
```

and include:

```text
<cluster_field>
```

corresponding to the cluster labels stored in the AnnData object.

### Notes

- Annotation columns are automatically transferred to `adata.obs`.
- Cluster IDs are converted into sequential labels and stored in `cluster_destination`.
- Full count matrices are restored before generating marker visualisations.
- Marker genes listed in `plot_markers` can be separated by commas or semicolons.
- Marker UMAPs are generated only for genes present in the AnnData object.
- A cluster-to-annotation mapping is stored in `adata.uns`.
- The final output is written as:

```text
<input>_annotated.h5ad
```
