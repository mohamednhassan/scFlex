.onLoad <- function(libname, pkgname) {
  if (requireNamespace("reticulate", quietly = TRUE)) {
    reticulate::py_require("anndata>=0.10")
  }
}
