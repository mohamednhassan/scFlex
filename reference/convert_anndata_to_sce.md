# Convert AnnData to SingleCellExperiment

Converts an AnnData H5AD file to a SingleCellExperiment object while
preserving counts, normalized expression when available, cell metadata,
feature metadata, and dimensional reductions.

## Usage

``` r
convert_anndata_to_sce(input, output)
```

## Arguments

- input:

  Path to an AnnData H5AD file.

- output:

  Path where the converted SingleCellExperiment RDS file should be
  written.

## Value

Invisibly returns the output file path.

## Details

An explicit `layers["counts"]` matrix is preferred as the counts assay.
Fractional values are accepted in an explicit counts layer as long as
all values are finite and non-negative. If no explicit counts layer is
present, `X` is used as counts only when its values are finite,
non-negative, and integer-like.

Normalized expression is optional. If `layers["logcounts"]` is present
it is stored as the SCE `logcounts` assay. Otherwise, `X` is used as
`logcounts` when it was not already used as the counts matrix.
