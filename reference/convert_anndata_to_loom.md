# Convert AnnData to Loom

Reads an AnnData H5AD file and writes a Loom file using loompy,
preserving supported matrices, cell and feature attributes, and
reductions where representable.

## Usage

``` r
convert_anndata_to_loom(input, output)
```

## Arguments

- input:

  Path to an input `.h5ad` file.

- output:

  Path to the output `.loom` file.

## Value

The output path, returned invisibly.

## Details

Loom is a more limited interchange format than Seurat,
SingleCellExperiment, or AnnData, so only components representable by
the current mapping are preserved.
