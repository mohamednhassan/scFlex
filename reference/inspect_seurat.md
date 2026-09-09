# Inspect a Seurat object

Internal helper that reports the structure of a Seurat object, including
its dimensions, assays, default assay structure, reductions, graphs,
neighbors, metadata columns, and a sample of cell names.

## Usage

``` r
inspect_seurat(obj)
```

## Arguments

- obj:

  A Seurat object.

## Value

Invisibly returns a list containing information about the Seurat object.
