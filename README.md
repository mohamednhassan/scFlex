# scFlex

**scFlex** is an R package for preservation-aware conversion among **Seurat**, **SingleCellExperiment**, **AnnData**, and **Loom**.

The goal is not merely to produce a file with a new extension. scFlex checks cell/feature alignment, distinguishes raw counts from normalized expression, preserves compatible metadata and embeddings, and reports when a target format cannot represent part of the source object.

## Overview

<p align="center">
  <img src="figures/scFlex_overview.png" width="900">
</p>

## Supported conversions

| From | To | Status |
|---|---|---|
| Seurat | AnnData | Supported |
| AnnData | Seurat | Supported |
| Seurat | SingleCellExperiment | Supported |
| SingleCellExperiment | Seurat | Supported |
| SingleCellExperiment | AnnData | Supported |
| AnnData | SingleCellExperiment | Supported |
| Seurat | Loom | Supported with Loom limitations |
| Loom | Seurat | Supported with Loom limitations |
| SingleCellExperiment | Loom | Supported with Loom limitations |
| Loom | SingleCellExperiment | Supported with Loom limitations |
| AnnData | Loom | Supported with Loom limitations |
| Loom | AnnData | Supported with Loom limitations |

> Loom has a smaller and increasingly legacy data model. scFlex supports it as an interchange format but does not claim lossless preservation of components Loom cannot represent.

## Installation

```r
# install.packages("remotes")
remotes::install_github("mohamednhassan/scFlex")
```

## Python setup

scFlex does **not** require a hard-coded Conda environment. It declares its Python requirements through `reticulate::py_require()` and lets reticulate resolve them in the user's Python configuration.

For AnnData conversion, scFlex declares `anndata>=0.10`. Loom conversion additionally declares `loompy>=3.0` only when Loom support is used.

## Basic usage

```r
library(scFlex)

convert_sc(
  input = "object.rds",
  output = "object.h5ad",
  from = "seurat",
  to = "anndata"
)
```

Format-specific wrappers are also available:

```r
convert_seurat_to_anndata("object.rds", "object.h5ad")
convert_anndata_to_seurat("object.h5ad", "object.rds")
convert_seurat_to_sce("object.rds", "object_sce.rds")
convert_sce_to_anndata("object_sce.rds", "object.h5ad")
convert_anndata_to_loom("object.h5ad", "object.loom")
```

## Inspect before converting

```r
inspect_sc("object.rds", format = "seurat")
inspect_sc("object.h5ad")
```

`inspect_sc()` reports the object structure without assigning a subjective score.

## Preservation model

Typical mappings include:

| Concept | Seurat | AnnData | SingleCellExperiment |
|---|---|---|---|
| Raw counts | `counts` | `layers["counts"]` | `assay("counts")` |
| Normalized expression | `data` | `X` | `assay("logcounts")` |
| Cell metadata | `meta.data` | `obs` | `colData` |
| Feature metadata | assay metadata | `var` | `rowData` |
| Embeddings | reductions | `obsm` | `reducedDims` |


## Seurat assay-structure helpers

```r
convert_seu_v5_to_classic(
  "object.rds",
  "object_classic.rds"
)

convert_seu_classic_to_v5(
  "object.rds",
  "object_v5.rds"
)
```

## Development

```r
devtools::document()
devtools::test()
devtools::check()
```

## Current scope

Version 0.1.0 focuses on expression matrices, cell/feature metadata, and dimensional reductions. Graphs, neighbors, Seurat command history, variable-feature state, feature loadings, multimodal `altExp`/MuData mapping, and spatial structures are not yet guaranteed to round-trip.

## License

MIT.

## Seurat v5 split layers

scFlex recognizes both canonical layers such as `counts`/`data` and split
Seurat v5 layers such as `counts.sample1`, `counts.sample2`, and `data.sample1`.
Matching split layers are joined internally on a temporary assay for conversion;
the input object is not modified.

Counts-only Seurat objects are also valid conversion inputs. When normalized
expression is absent, AnnData `X` is populated with the raw counts without
normalization, and the conversion result reports this explicitly.
