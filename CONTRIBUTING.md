# Contributing to scFlex

Issues and pull requests are welcome. For conversion bugs, please
include:

1.  Source and target formats.
2.  [`sessionInfo()`](https://rdrr.io/r/utils/sessionInfo.html).
3.  [`reticulate::py_config()`](https://rstudio.github.io/reticulate/reference/py_config.html)
    when AnnData or Loom is involved.
4.  The output of
    [`inspect_sc()`](https://mohamednhassan.github.io/scFlex/reference/inspect_sc.md)
    when possible.
5.  A minimal reproducible object whenever data sharing permits.

New conversion features should include tests for cell order, feature
order, matrix equality, and expected metadata/reduction preservation.
