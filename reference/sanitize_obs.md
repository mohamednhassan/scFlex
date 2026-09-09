# Sanitize cell metadata for AnnData

Checks cell metadata columns and simplifies supported one-column matrix
or data-frame columns before conversion to AnnData obs.

## Usage

``` r
sanitize_obs(obs)
```

## Arguments

- obs:

  A data frame containing cell-level metadata with cells as row names.

## Value

A sanitized data frame with the original row names preserved.
