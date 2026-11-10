# scribble
A single-cell RNA-seq workflow based on scanpy.


## Getting Started

- [Import](docs/import.md)
- [Mitochondrial QC](docs/mt.md)
- [Doublet Detection](docs/dbl.md)
- [PCA QC](docs/pca.md)
- [Filtering](docs/filter.md)

---

# Core Analysis Workflow

## Preprocessing and Integration

- [Preintegration](docs/preintegration.md)
- [Harmony Integration](docs/harmony.md)

## Clustering and Curation

- [Clustering](docs/cluster.md)
- [Cluster Evaluation](docs/evaluate.md)
- [Cluster Refinement](docs/refine.md)
- [Cluster Curation](docs/curate.md)

## Biological Annotation

- [Annotation](docs/annotate.md)

---

# Differential Expression and Pathway Analysis

## Differential Expression

- [DGE](docs/dge.md)

## Cross-Analysis Comparison

- [DGE Compare](docs/dge_compare.md)

## Interpretation and Reporting

- [DGE Report](docs/dge_report.md)

---


# Workflow Summary

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
  ↓
DGE compare
  ↓
DGE report
```

---

# Module Categories

| Category | Modules |
|-----------|----------|
| Data import | import |
| Quality control | mt, dbl, pca, filter |
| Integration | preintegration, harmony |
| Clustering | cluster, evaluate, refine, curate |
| Annotation | annotate |
| Differential expression | DGE |
| Comparison and reporting | DGE compare, DGE report |

---

# Typical Analysis Path

## Standard single-dataset workflow

```bash
scribble import
scribble mt
scribble dbl
scribble pca
scribble filter
scribble preintegration
scribble harmony
scribble cluster
scribble evaluate
scribble refine
scribble curate
scribble annotate
```

## Differential expression workflow

```bash
python DGE_v3.py
python DGE_compare.py
python DGE_report.py
```


<img style="float:right;width:200px;" src="./img/scribble.png" alt="scribble"/>
