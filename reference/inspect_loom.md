# Inspect a Loom file

Internal helper that reads the HDF5 structure of a Loom file and reports
its dimensions, layers, cell and feature attributes, graphs, and a
sample of cell names when available.

## Usage

``` r
inspect_loom(obj)
```

## Arguments

- obj:

  Path to a Loom file.

## Value

Invisibly returns a list containing information about the Loom file.
