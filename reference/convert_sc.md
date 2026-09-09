# Convert between supported single-cell formats

Dispatches a file conversion between Seurat, SingleCellExperiment,
AnnData, and Loom using the corresponding scFlex conversion function.

## Usage

``` r
convert_sc(input, source, destination, output)
```

## Arguments

- input:

  Character string giving the path to the input file.

- source:

  Character string naming the source format. Supported values are
  `"seurat"`, `"sce"`, `"anndata"`, and `"loom"`. Matching is
  case-insensitive.

- destination:

  Character string naming the destination format. Supported values are
  `"seurat"`, `"sce"`, `"anndata"`, and `"loom"`. Matching is
  case-insensitive.

- output:

  Character string giving the output file path.

## Value

The result returned by the selected conversion function, invisibly.
Current converters normally return the output path invisibly; some
conversions may return `NULL` invisibly when required count data cannot
be identified.

## Details

Source and destination must differ. Unsupported format names or
unimplemented conversion paths produce an error.
