# Convert SingleCellExperiment to Loom

Reads a SingleCellExperiment RDS object and writes a Loom file using
loompy, preserving supported assays, attributes, and reductions where
representable.

## Usage

``` r
convert_sce_to_loom(input, output)
```

## Arguments

- input:

  Path to an input SingleCellExperiment `.rds` file.

- output:

  Path to the output `.loom` file.

## Value

The output path, returned invisibly.

## Details

Loom is a more limited interchange format than Seurat,
SingleCellExperiment, or AnnData, so only components representable by
the current mapping are preserved.
