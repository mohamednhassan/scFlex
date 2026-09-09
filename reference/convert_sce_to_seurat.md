# Convert SingleCellExperiment to Seurat

Reads a SingleCellExperiment RDS object and writes a Seurat RDS object
while preserving supported counts, normalized expression, metadata,
feature metadata, and dimensional reductions.

## Usage

``` r
convert_sce_to_seurat(input, output)
```

## Arguments

- input:

  Path to an input SingleCellExperiment `.rds` file.

- output:

  Path to the output Seurat `.rds` file.

## Value

The output path, returned invisibly.
