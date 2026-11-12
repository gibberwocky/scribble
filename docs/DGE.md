# DGE

## Purpose

The `DGE` script performs pseudobulk differential expression (DE) analysis and gene set enrichment analysis (GSEA) for annotated cell populations.

Cells are aggregated to the sample level, statistical testing is performed using either PyDESeq2 or edgeR, and pathway enrichment is assessed using pre-ranked GSEA with gseapy. The module also generates extensive model diagnostics, QC reports, sample-level PCA analyses, and fit-quality assessments.

This module is intended for downstream biological interpretation of annotated cell populations rather than routine preprocessing.

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
  ↓
DGE
```

## Inputs

- Annotated AnnData object
- Cell-type annotation column
- Experimental design formula
- Contrast definition
- Sample metadata
- Raw counts layer

## Outputs

- Differential expression tables
- Significant and nominal gene lists
- GSEA pathway enrichment results
- Cell-type ranking files
- Sample-level pseudobulk count matrices
- QC and model diagnostics
- PCA plots
- Sample correlation plots
- Inspection reports
- Analysis summary table
- Run manifest

## Main Functions

- Build pseudobulk count matrices for each annotated cell population
- Validate experimental design and replication structure
- Run differential expression with PyDESeq2 or edgeR
- Generate ranked gene lists
- Perform pathway enrichment analysis using GSEA
- Calculate sample-level QC metrics
- Assess model convergence and fit quality
- Generate pseudobulk PCA diagnostics
- Generate sample correlation diagnostics
- Produce cell-type specific summary statistics
- Compile a global DE/GSEA summary report

## Parameters

### Required Parameters

| Parameter | Description |
|------------|-------------|
| `--project_dir` | Root Scribble project directory. |
| `--input_file` | Annotated AnnData object containing raw counts and sample metadata. |
| `--annotation` | Annotation column used to define cell populations for pseudobulk analysis. |
| `--analysis_name` | Name used to label analysis outputs and directories. |
| `--formula` | Design formula used for differential expression modelling (e.g. `~ donor + treatment`). |
| `--contrast` | Differential expression contrast specified as `VARIABLE TEST_LEVEL REFERENCE_LEVEL`. |

### Differential Expression Method

| Parameter | Default | Description |
|------------|----------|-------------|
| `--method` | pydeseq2 | Differential expression method. Available options: `pydeseq2`, `edger`. |

### Sample Selection

| Parameter | Default | Description |
|------------|----------|-------------|
| `--samples` | None | Restrict analysis to specific samples (by default all samples are analysed). |

### Gene Filtering

| Parameter | Default | Description |
|------------|----------|-------------|
| `--min_gene_count` | 10 | Minimum pseudobulk count required for a gene to be considered expressed. |
| `--min_gene_samples` | 3 | Minimum number of samples in which a gene must satisfy `min_gene_count`. |

### Pseudobulk Inclusion

| Parameter | Default | Description |
|------------|----------|-------------|
| `--min_cells_per_sample` | 20 | Minimum cells required in a sample for a given cell type before pseudobulk generation. |
| `--min_samples_per_group` | 2 | Minimum samples required in each contrast group for DE testing. |

### Donor-Aware Design Assessment

| Parameter | Default | Description |
|------------|----------|-------------|
| `--min_informative_donors` | 3 | Minimum donors represented in both contrast groups. |
| `--require_informative_donors` | True | Require sufficient donor overlap before DE analysis proceeds. Use `--no-require_informative_donors` to disable. |

### Gene Set Enrichment Analysis

| Parameter | Default | Description |
|------------|----------|-------------|
| `--gsea_permutations` | 10000 | Number of permutations used by GSEA. |

### Parallelisation

| Parameter | Default | Description |
|------------|----------|-------------|
| `--n_cpus` | 8 | Number of CPUs used for DE and GSEA calculations. |

### Supported Gene Sets

The following pathway collections are analysed automatically:

| Collection |
|------------|
| Hallmark |
| KEGG |
| GO Biological Process |
| Reactome |
| WikiPathways |

### Notes

- Differential expression is performed independently for each annotated cell population.
- Pseudobulk counts are generated at the sample level.
- Raw counts are required and are restored automatically before analysis.
- PyDESeq2 uses Wald statistics for gene ranking during GSEA.
- Significant genes are reported using:
  - FDR < 0.05
  - optional effect-size thresholds (`|log2FC| > 1`)
- GSEA is performed using pre-ranked statistics rather than significance-filtered genes.
- Extensive QC diagnostics are generated automatically for every analysed cell type.
- Cell types that fail replication or donor-design requirements are skipped and reported in the summary output.
