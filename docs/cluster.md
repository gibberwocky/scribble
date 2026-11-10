# cluster

## Purpose

The `cluster` module performs Leiden clustering, identifies cluster marker genes, evaluates clustering stability, and generates diagnostic outputs to support downstream cluster curation.

In addition to assigning cluster labels, the module can automatically optimise Leiden resolution, merge clusters below a minimum size threshold, quantify clustering stability across repeated runs, and generate detailed cluster summaries.

The resulting clustered AnnData object forms the basis for cluster evaluation, refinement, annotation, and downstream biological interpretation.

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
annotate
```

## Inputs

- Harmony-integrated AnnData object
- Low-dimensional embedding (typically `X_pca_harmony`)
- Sample annotations
- Full expression matrix stored in `adata.raw`

## Outputs

- Clustered AnnData object
- Marker gene workbook
- Cluster stability metrics
- Cluster reassignment statistics
- Cluster summary table
- Cluster UMAP visualisations
- Resolution optimisation diagnostics
- Module manifest

## Main Functions

- Construct a neighbour graph from a chosen embedding
- Run Leiden clustering
- Optionally optimise clustering resolution automatically
- Merge clusters below a minimum size threshold
- Calculate cluster stability across repeated clustering runs
- Quantify cluster reassignment behaviour
- Identify marker genes using Wilcoxon rank-sum testing
- Generate cluster-level QC and diagnostic reports
- Create a clustered AnnData object for downstream curation and annotation

## Parameters

### Required Parameters

| Parameter | Description |
|------------|-------------|
| `--project_dir` | Root Scribble project directory. |
| `--input` | Harmony-integrated AnnData object. |
| `--embedding` | Embedding stored in `adata.obsm` used for neighbour graph construction and clustering (typically `X_pca_harmony`). |

### Leiden Clustering

| Parameter | Default | Description |
|------------|----------|-------------|
| `--resolution` | 1.0 | Leiden clustering resolution. Ignored if `--auto_resolution` is enabled. |
| `--neighbors` | 15 | Number of neighbours used when constructing the graph. |
| `--min_cluster_size` | 50 | Clusters smaller than this threshold are merged into the most strongly connected neighbouring cluster. |

### Automatic Resolution Optimisation

| Parameter | Default | Description |
|------------|----------|-------------|
| `--auto_resolution` | False | Automatically identify an optimal Leiden resolution using silhouette scores. |
| `--res_min` | 0.2 | Minimum resolution tested during optimisation. |
| `--res_max` | 2.0 | Maximum resolution tested during optimisation. |
| `--res_steps` | 10 | Number of resolution values evaluated in each optimisation pass. |
| `--fine_width` | 0.2 | Width of the fine-resolution search around the best coarse resolution. |
| `--max_cells` | 20000 | Maximum cells used when calculating silhouette scores. |
| `--max_dims` | 50 | Maximum dimensions retained from the embedding when performing optimisation. |

### Stability Analysis

| Parameter | Default | Description |
|------------|----------|-------------|
| `--n_repeats` | 10 | Number of repeated Leiden clustering runs used for stability estimation. |
| `--dominant_reassignment_threshold` | 0.05 | Minimum reassignment fraction reported in cluster summaries. |
| `--reassignment_display_threshold` | 0.01 | Minimum reassignment fraction included in reassignment summaries. |
| `--high_reassignment` | 0.30 | Threshold used to flag potential merge candidates. |

### Marker Detection

| Parameter | Default | Description |
|------------|----------|-------------|
| `--nmarkers` | 100 | Number of marker genes exported per cluster. |

### Visualisation

| Parameter | Description |
|------------|-------------|
| `--vars` | Metadata variables to display on UMAP plots. |

### Output Control

| Parameter | Default | Description |
|------------|----------|-------------|
| `--output_prefix` | Input filename stem | Prefix used when naming output files. |

### Reproducibility

| Parameter | Default | Description |
|------------|----------|-------------|
| `--manifest` | None | Load parameters from a previous Scribble manifest file. |

### Notes

- Leiden clustering is performed using the `igraph` backend.
- Small clusters are automatically merged based on graph connectivity.
- Automatic resolution optimisation uses silhouette scores calculated from the selected embedding.
- Cluster stability is estimated by repeating Leiden clustering with multiple random seeds.
- Marker genes are identified using the Wilcoxon rank-sum test.
- Marker statistics include log fold change, adjusted p-value, detection frequency, mean expression, and cluster size.
- Cluster summaries include stability, sample entropy, reassignment statistics, and merge-candidate suggestions.
