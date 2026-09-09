# Inspect an R single-cell object

Internal helper that identifies whether an object read from an RDS file
is a Seurat or SingleCellExperiment object and dispatches it to the
corresponding inspection function.

## Usage

``` r
inspect_rds_object(obj)
```

## Arguments

- obj:

  An R object read from an RDS file.

## Value

Invisibly returns a list containing information about the object.
