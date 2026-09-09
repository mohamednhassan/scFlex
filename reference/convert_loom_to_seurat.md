# Convert Loom to Seurat

Reads a Loom file through AnnData and creates a Seurat object while
preserving supported count data, metadata, feature metadata, and
reductions.

## Usage

``` r
convert_loom_to_seurat(input, output)
```

## Arguments

- input:

  Path to an input `.loom` file.

- output:

  Path to the output Seurat `.rds` file.

## Value

The output path, returned invisibly on success. If suitable raw count
data cannot be identified, the function may return `NULL` invisibly
without writing an output object.

## Details

Raw counts are selected conservatively from an explicit `counts` layer
when available or from `X` when it is non-negative and integer-like. The
function does not guess raw counts from arbitrary named layers. Loom is
a more limited interchange format than Seurat, SingleCellExperiment, or
AnnData, so only components representable by the current mapping are
preserved.
