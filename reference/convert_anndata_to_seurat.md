# Convert AnnData to Seurat

Reads an AnnData H5AD file and writes a Seurat RDS object while
preserving supported counts, normalized expression, metadata, feature
metadata, and dimensional reductions.

## Usage

``` r
convert_anndata_to_seurat(input, output)
```

## Arguments

- input:

  Path to an input `.h5ad` file.

- output:

  Path to the output Seurat `.rds` file.

## Value

The output path, returned invisibly.
