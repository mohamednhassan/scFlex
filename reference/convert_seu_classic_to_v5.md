# Convert a classic Seurat Assay to Assay5

Converts a selected classic Seurat `Assay` to a Seurat v5 `Assay5` while
preserving the supported assay data and object structure.

## Usage

``` r
convert_seu_classic_to_v5(input, output, assay = "RNA")
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
