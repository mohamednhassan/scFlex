# scFlex

**scFlex** is an R package for preservation-aware conversion among
**Seurat**, **SingleCellExperiment**, **AnnData**, and **Loom**.

The goal is not merely to produce a file with a new extension. scFlex
checks cell/feature alignment, distinguishes raw counts from normalized
expression, preserves compatible metadata and embeddings, and reports
when a target format cannot represent part of the source object.

## Overview

![](figures/scFlex_overview.png)

## Supported conversions

| From                 | To                   | Status                          |
|----------------------|----------------------|---------------------------------|
| Seurat               | AnnData              | Supported                       |
| AnnData              | Seurat               | Supported                       |
| Seurat               | SingleCellExperiment | Supported                       |
| SingleCellExperiment | Seurat               | Supported                       |
| SingleCellExperiment | AnnData              | Supported                       |
| AnnData              | SingleCellExperiment | Supported                       |
| Seurat               | Loom                 | Supported with Loom limitations |
| Loom                 | Seurat               | Supported with Loom limitations |
| SingleCellExperiment | Loom                 | Supported with Loom limitations |
| Loom                 | SingleCellExperiment | Supported with Loom limitations |
| AnnData              | Loom                 | Supported with Loom limitations |
| Loom                 | AnnData              | Supported with Loom limitations |

> Loom has a smaller and increasingly legacy data model. scFlex supports
> it as an interchange format but does not claim lossless preservation
> of components Loom cannot represent.

## Installation

``` r

# install.packages("remotes")
remotes::install_github("mohamednhassan/scFlex")
```

## Python setup

scFlex does **not** require a hard-coded Conda environment. It declares
its Python requirements through
[`reticulate::py_require()`](https://rstudio.github.io/reticulate/reference/py_require.html)
and lets reticulate resolve them in the user’s Python configuration.

For AnnData conversion, scFlex declares `anndata>=0.10`. Loom conversion
additionally declares `loompy>=3.0` only when Loom support is used.

## Basic usage

### 1. Inspect the object first

Before conversion, inspect the input object to understand its structure
and available components.

``` r

inspect_sc("object.rds")
inspect_sc("object.h5ad")
inspect_sc("object.loom")
```

[`inspect_sc()`](https://mohamednhassan.github.io/scFlex/reference/inspect_sc.md)
reports the detected object structure, including information such as
assays/layers, dimensions, metadata, and dimensional reductions, without
assigning a subjective conversion score.

### 2. Convert with `convert_sc()`

After inspection, use the general
[`convert_sc()`](https://mohamednhassan.github.io/scFlex/reference/convert_sc.md)
interface:

``` r

convert_sc(
  input = "object.rds",
  output = "object.h5ad",
  source = "seurat",
  destination = "anndata"
)
```

Another example:

``` r

convert_sc(
  input = "object.h5ad",
  output = "object_sce.rds",
  source = "anndata",
  destination = "sce"
)
```

Format-specific conversion functions are also available:

``` r

convert_seurat_to_anndata("object.rds", "object.h5ad")
convert_anndata_to_seurat("object.h5ad", "object.rds")

convert_seurat_to_sce("object.rds", "object_sce.rds")
convert_sce_to_seurat("object_sce.rds", "object.rds")

convert_sce_to_anndata("object_sce.rds", "object.h5ad")
convert_anndata_to_sce("object.h5ad", "object_sce.rds")

convert_anndata_to_loom("object.h5ad", "object.loom")
convert_loom_to_anndata("object.loom", "object.h5ad")
```

## Seurat assay conversion

scFlex also provides helper functions for converting between classic
Seurat `Assay` objects and Seurat v5 `Assay5` objects.

### Seurat v5 Assay5 to classic Assay

``` r

convert_seu_v5_to_classic(
  input = "object.rds",
  output = "object_classic.rds",
  assay = "RNA"
)
```

### Classic Seurat Assay to Seurat v5 Assay5

``` r

convert_seu_classic_to_v5(
  input = "object.rds",
  output = "object_v5.rds",
  assay = "RNA"
)
```

These functions are useful when working with tools or workflows that
expect a particular Seurat assay structure.

## Preservation model

Typical mappings include:

| Concept | Seurat | AnnData | SingleCellExperiment |
|----|----|----|----|
| Raw counts | `counts` | `layers["counts"]` | `assay("counts")` |
| Normalized expression | `data` | `X` | `assay("logcounts")` |
| Cell metadata | `meta.data` | `obs` | `colData` |
| Feature metadata | assay metadata | `var` | `rowData` |
| Embeddings | reductions | `obsm` | `reducedDims` |

## Seurat v5 split layers

scFlex recognizes both canonical layers such as `counts`/`data` and
split Seurat v5 layers such as:

``` text
counts.sample1
counts.sample2
data.sample1
data.sample2
```

Matching split layers are joined internally on a temporary assay for
conversion; the input object is not modified.

Counts-only Seurat objects are also valid conversion inputs. When
normalized expression is absent, AnnData `X` is populated with the raw
counts without normalization, and the conversion result reports this
explicitly.

## Development

``` r

devtools::document()
devtools::test()
devtools::check()
```

## Current scope

Version 0.1.0 focuses on expression matrices, cell/feature metadata, and
dimensional reductions. Graphs, neighbors, Seurat command history,
variable-feature state, feature loadings, multimodal `altExp`/MuData
mapping, and spatial structures are not yet guaranteed to round-trip.

## License

MIT.
