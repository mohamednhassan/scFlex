# Convert Seurat to AnnData

Reads a Seurat RDS object and writes an H5AD AnnData file while
preserving supported counts, normalized expression, metadata, feature
metadata, and dimensional reductions.

## Usage

``` r
convert_seurat_to_anndata(input, output, assay = "RNA")
```

## Arguments

- input:

  Path to an input Seurat `.rds` file.

- output:

  Path to the output `.h5ad` file.

- assay:

  Name of the Seurat assay to convert. Defaults to `"RNA"`.

## Value

The output path, returned invisibly.
