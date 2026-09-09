# Convert a Seurat Assay5 assay to a classic Assay

Converts a selected Seurat v5 `Assay5` assay to a classic Seurat `Assay`
while preserving the supported assay data and object structure.

## Usage

``` r
convert_seu_v5_to_classic(input, output, assay = "RNA")
```

## Arguments

- input:

  Path to an input Seurat `.rds` file.

- output:

  Path to the output Seurat `.rds` file.

- assay:

  Name of the assay to convert. Defaults to `"RNA"`.

## Value

The output path, returned invisibly.
