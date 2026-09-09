# Convert Seurat to SingleCellExperiment

Reads a Seurat RDS object and writes a SingleCellExperiment RDS object
while preserving supported counts, normalized expression, metadata,
feature metadata, and dimensional reductions.

## Usage

``` r
convert_seurat_to_sce(input, output)
```

## Arguments

- input:

  Path to an input Seurat `.rds` file.

- output:

  Path to the output SingleCellExperiment `.rds` file.

## Value

The output path, returned invisibly.
