#' Sanitize cell metadata for AnnData
#'
#' Checks cell metadata columns and simplifies supported one-column matrix or data-frame columns before conversion to AnnData obs.
#'
#' @param obs A data frame containing cell-level metadata with cells as row names.
#'
#' @return A sanitized data frame with the original row names preserved.
#' @keywords internal
sanitize_obs <- function(obs) {
  cell_names <- rownames(obs)
  obs[] <- lapply(names(obs), function(col) {
    x <- obs[[col]]
    if (is.factor(x)) {
      return(x)
    }
    if (is.atomic(x) && is.null(dim(x))) {
      return(x)
    }
    if (is.matrix(x) && ncol(x) == 1) {
      return(as.vector(x[, 1]))
    }
    if (is.data.frame(x) && ncol(x) == 1) {
      return(x[[1]])
    }
    stop("Metadata column '", col, "' cannot be directly converted to AnnData obs. ",
         "Class: ", paste(class(x), collapse = ", "),
         if (!is.null(dim(x))) {
           paste0("; dimensions: ", paste(dim(x), collapse = " x "))
         } else {
           ""
         }
    )
  })
  rownames(obs) <- cell_names
  obs
}

#' Collapse values for display
#'
#' Collapses a vector into a comma-separated string and returns `"None"` for null or empty input.
#'
#' @param x A vector-like object to collapse.
#'
#' @return A single character string.
#' @keywords internal
collapse_or_none <- function(x) {
  if (is.null(x) || length(x) == 0) {
    return("None")
  }
  paste(x, collapse = ", ")
}
