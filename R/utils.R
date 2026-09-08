#' Sanitize cell metadata for AnnData
#'
#' Checks cell metadata columns and simplifies supported one-column matrix or data-frame columns before conversion to AnnData obs.
#'
#' @param obs A data frame containing cell-level metadata with cells as row names.
#'
#' @return A sanitized data frame with the original row names preserved.
#' @keywords internal
sanitize_obs <- function(obs) {
  ## keep original cell names
  cell_names <- rownames(obs)
  obs[] <- lapply(names(obs), function(col) {
    x <- obs[[col]]
    ## Preserve factors
    if (is.factor(x)) {
      return(x)
    }
    ## Preserve simple one-dimensional atomic vectors
    if (is.atomic(x) && is.null(dim(x))) {
      return(x)
    }
    ## Flatten a one-column matrix to a vector by accessing the first column
    if (is.matrix(x) && ncol(x) == 1) {
      return(as.vector(x[, 1]))
    }
    ## Extract the vector from a one-column data frame
    if (is.data.frame(x) && ncol(x) == 1) {
      return(x[[1]])
    }
    ## Otherwise, raise an error, identify the column, its class, and dimensions (if applicable)
    stop("Metadata column '", col, "' cannot be directly converted to AnnData obs. ",
         "Class: ", paste(class(x), collapse = ", "),
         if (!is.null(dim(x))) {
           paste0("; dimensions: ", paste(dim(x), collapse = " x "))
         } else {
           ""
         }
    )
  })
  ## Put the cell names back as rownames
  rownames(obs) <- cell_names
  ## Return the updated dataframe
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
  ## Check if the vector is null or empty
  if (is.null(x) || length(x) == 0) {
    return("None")
  }
  ## Otherwise return, in a comma separated string
  paste(x, collapse = ", ")
}
