# Inspect a SingleCellExperiment object

Internal helper that reports the structure of a SingleCellExperiment
object, including its dimensions, assays, reduced dimensions,
alternative experiments, cell metadata, feature metadata, and a sample
of cell names.

## Usage

``` r
inspect_sce(obj)
```

## Arguments

- obj:

  A SingleCellExperiment object.

## Value

Invisibly returns a list containing information about the
SingleCellExperiment object.
