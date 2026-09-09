# Inspect an AnnData object

Internal helper that reads an H5AD file and reports its dimensions,
primary matrix, layers, cell-level matrices, graphs, feature-level
matrices, unstructured metadata, raw data, metadata columns, and a
sample of cell names.

## Usage

``` r
inspect_anndata(path_to_file)
```

## Arguments

- path_to_file:

  Path to an H5AD file.

## Value

Invisibly returns a list containing information about the AnnData
object.
