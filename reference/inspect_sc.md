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
if (requireNamespace("Seurat", quietly = TRUE)) {
  counts <- matrix(
    c(1, 0, 3, 0, 2, 1),
    nrow = 2,
    dimnames = list(
      c("Gene1", "Gene2"),
      c("Cell1", "Cell2", "Cell3")
    )
  )
  obj <- Seurat::CreateSeuratObject(counts = counts)
  path <- tempfile(fileext = ".rds")
  saveRDS(obj, path)
  inspect_sc(path)
  unlink(path)
}
#> Warning: Data is of class matrix. Coercing to dgCMatrix.
#> RDS file detected.
#> Seurat object detected.
#> Version of Seurat: 5.4.0
#> Number of cells: 3
#> Number of features: 2
#> Assays [1]: RNA
#> Default assay: RNA
#> Default assay structure: Assay5
#> Default assay layers [1]: counts
#> Reductions [0]: None
#> Graphs [0]: None
#> Neighbors [0]: None
#> Metadata columns [3]: orig.ident, nCount_RNA, nFeature_RNA
#> Sample of cell names: Cell1, Cell3, Cell2
```
