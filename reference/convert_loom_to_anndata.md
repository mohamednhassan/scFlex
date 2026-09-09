# Convert Loom to AnnData

Reads a Loom file through AnnData and writes the resulting object as
H5AD while preserving components supported by the Loom-to-AnnData
reader.

## Usage

``` r
convert_loom_to_anndata(input, output)
```

## Arguments

- input:

  Path to an input `.loom` file.

- output:

  Path to the output `.h5ad` file.

## Value

The output path, returned invisibly.

## Details

Loom is a more limited interchange format than Seurat,
SingleCellExperiment, or AnnData, so only components representable by
the current mapping are preserved.
