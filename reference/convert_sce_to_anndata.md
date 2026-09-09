# Convert SingleCellExperiment to AnnData

Reads a SingleCellExperiment RDS object and writes an H5AD AnnData file,
mapping supported assays, metadata, feature metadata, and reduced
dimensions.

## Usage

``` r
convert_sce_to_anndata(input, output)
```

## Arguments

- input:

  Path to an input SingleCellExperiment `.rds` file.

- output:

  Path to the output `.h5ad` file.

## Value

The output path, returned invisibly.
