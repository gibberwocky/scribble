# curate

## Purpose

The `curate` module applies manual curation decisions to an existing clustering result.

It supports:

- Merging biologically similar clusters
- Refining selected clusters through an additional round of subclustering
- Maintaining cluster lineage relationships
- Recalculating markers after curation
- Generating curated cluster summaries and reports

The output is a curated clustering hierarchy intended for final biological annotation.

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

- Refined or clustered AnnData object
- Cluster refinement decisions
- Existing cluster stability metrics
- Existing marker tables
- Existing lineage information

## Outputs

- Curated AnnData object
- Curated marker workbook
- Curation report
- Cluster lineage report
- Curated UMAP
- Optional reassignment and marker tables
- Module manifest

## Main Functions

- Merge selected clusters
- Refine selected clusters through recursive reclustering
- Track cluster ancestry and lineage paths
- Recalculate markers after curation
- Generate curated cluster summaries
- Record curation history
- Produce final curated cluster labels for annotation

## Parameters

### Required Parameters

| Parameter | Description |
|------------|-------------|
| `--project_dir` | Root Scribble project directory. |
| `--input` | Refined or clustered AnnData object to curate. |

### Cluster Refinement

| Parameter | Default | Description |
|------------|----------|-------------|
| `--refine` | None | Comma-separated list of clusters to subcluster (e.g. `R8,R10`). |
| `--max_refine_depth` | 5 | Maximum recursive refinement depth when refinement is performed. |
| `--min_cells_per_group` | 200 | Minimum cells required before a cluster can be refined. |
| `--min_cells_per_cluster` | 50 | Minimum size for refined subclusters. Smaller clusters are collapsed back into the parent. |

### Cluster Merging

| Parameter | Default | Description |
|------------|----------|-------------|
| `--merge` | None | Merge groups specified as `clusterA+clusterB,clusterC+clusterD`. The first cluster becomes the merge target. |

### Re-Clustering Parameters

Used only when clusters are refined.

| Parameter | Default | Description |
|------------|----------|-------------|
| `--hvgs` | 3000 | Number of highly variable genes used during refinement. |
| `--batch` | sample | Batch variable used for Harmony integration. |
| `--npcs` | 50 | Number of principal components. |
| `--neighbors` | 15 | Number of nearest neighbours. |
| `--theta` | 1 | Harmony batch-correction parameter. |
| `--resolution` | 0.5 | Leiden clustering resolution. |

### Automatic Resolution Optimisation

| Parameter | Default | Description |
|------------|----------|-------------|
| `--auto_resolution` | False | Automatically optimise Leiden resolution during refinement. |
| `--res_min` | 0.2 | Minimum resolution tested. |
| `--res_max` | 2.0 | Maximum resolution tested. |
| `--res_steps` | 10 | Number of candidate resolutions tested. |
| `--fine_width` | 0.2 | Width of secondary fine search around the optimal resolution. |

### Marker Detection

| Parameter | Default | Description |
|------------|----------|-------------|
| `--nmarkers` | 100 | Number of marker genes exported per cluster. |
| `--marker_strength_threshold` | 1.0 | Threshold used to determine whether markers are sufficiently distinct. |

### Stability Analysis

| Parameter | Default | Description |
|------------|----------|-------------|
| `--n_repeats` | 10 | Number of repeated clustering runs used for stability estimation. |
| `--stability_threshold` | 0.95 | Stability threshold above which clusters are considered sufficiently resolved. |
| `--dominant_reassignment_threshold` | 0.05 | Minimum reassignment fraction reported in summaries. |
| `--reassignment_display_threshold` | 0.01 | Minimum reassignment fraction displayed in reassignment tables. |
| `--high_reassignment` | 0.30 | Threshold used to flag potential merge candidates. |

### Integration Behaviour

| Parameter | Default | Description |
|------------|----------|-------------|
| `--skip_integration` | False | Skip Harmony integration and cluster directly in PCA space. |
| `--no_scale` | False | Skip gene scaling before PCA. |

### Saved Parameter Behaviour

| Parameter | Default | Description |
|------------|----------|-------------|
| `--ignore_saved_params` | False | Ignore clustering parameters stored in the AnnData object and use CLI values instead. |

### Reproducibility

| Parameter | Default | Description |
|------------|----------|-------------|
| `--manifest` | None | Load parameters from a previous Scribble manifest file. |

### Notes

- A cluster cannot be merged and refined in the same curate run.
- Refined clusters receive globally unique labels (e.g. `R8_S1`, `R8_S2`) during processing and are subsequently incorporated into the curated hierarchy.
- Lineage relationships are preserved and written to the curation lineage report.
- Marker genes are recalculated after curation using the final curated cluster labels.
- Curate operations are appended to the curation history stored within the AnnData object.
