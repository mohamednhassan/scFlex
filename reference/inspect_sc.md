# Inspect a single-cell object

Inspects the structure and contents of a supported single-cell data file
without performing a conversion. Supported inputs include Seurat and
SingleCellExperiment objects stored as RDS files, AnnData H5AD files,
and Loom files.

## Usage

``` r
inspect_sc(path_to_file)
```

## Arguments

- path_to_file:

  Path to an RDS, H5AD, or Loom file.

## Value

Invisibly returns a list containing information about the inspected
object. The contents of the list depend on the input format.

## Details

The function reports basic information such as the number of cells and
features, available assays or layers, dimensional reductions, metadata,
and a sample of cell names. Additional format-specific information is
reported when available.

## Examples

``` r
if (FALSE) { # \dontrun{
inspect_sc("example.rds")
inspect_sc("example.h5ad")
inspect_sc("example.loom")
} # }
```
