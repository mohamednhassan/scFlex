.read_anndata_file <- function(input) {
  ad <- .import_anndata()
  tryCatch(
    ad$read_h5ad(input),
    error = function(e) stop("Failed to read H5AD file:\n", conditionMessage(e), call. = FALSE)
  )
}

.read_loom_file <- function(input) {
  ad <- .import_anndata(require_loom = TRUE)
  reader <- NULL
  if (!is.null(ad$read_loom)) reader <- ad$read_loom
  if (is.null(reader) && !is.null(ad$io$read_loom)) reader <- ad$io$read_loom
  if (is.null(reader)) stop("The installed anndata version does not expose a Loom reader.", call. = FALSE)

  tryCatch(
    reader(input, sparse = TRUE),
    error = function(e) stop("Failed to read Loom file:\n", conditionMessage(e), call. = FALSE)
  )
}

.write_anndata_file <- function(adata, output) {
  tryCatch(
    adata$write_h5ad(output),
    error = function(e) stop("Failed to write H5AD file:\n", conditionMessage(e), call. = FALSE)
  )
  invisible(output)
}

.write_loom_file <- function(adata, output, include_reductions = TRUE) {
  .import_anndata(require_loom = TRUE)
  tryCatch(
    adata$write_loom(output, write_obsm_varm = isTRUE(include_reductions)),
    error = function(e) stop("Failed to write Loom file:\n", conditionMessage(e), call. = FALSE)
  )
  invisible(output)
}
