#' Inspect a single-cell object or file
#'
#' Reports the major components available in Seurat, SingleCellExperiment,
#' AnnData, or Loom input without modifying the source.
#'
#' @param input Path to an `.rds`, `.h5ad`, or `.loom` file.
#' @param format Optional source format: `"seurat"`, `"sce"`, `"anndata"`, or `"loom"`.
#' @param assay Seurat assay to inspect. Default is `"RNA"`.
#' @param counts_layer Name used for raw counts in AnnData/Loom or Seurat Assay5.
#' @param data_layer Name used for normalized data in Seurat.
#' @return Invisibly returns a named list describing the source.
#' @export
inspect_sc <- function(input, format = NULL, assay = "RNA", counts_layer = "counts", data_layer = "data") {
  .validate_input_file(input)
  format <- .detect_format(input, format)

  info <- list(format = format)
  message("Source format: ", format)

  if (format == "seurat") {
    obj <- readRDS(input)
    assays <- SeuratObject::Assays(obj)
    info$cells <- ncol(obj)
    info$features <- nrow(obj)
    info$assays <- assays
    info$reductions <- SeuratObject::Reductions(obj)
    info$graphs <- SeuratObject::Graphs(obj)
    info$neighbors <- SeuratObject::Neighbors(obj)
    info$metadata_columns <- colnames(obj[[]])
    info$assay_present <- assay %in% assays

    if (info$assay_present) {
      a <- obj[[assay]]
      info$assay_class <- class(a)[1L]
      info$layers <- if (inherits(a, "Assay5")) SeuratObject::Layers(a) else methods::slotNames(a)
      if (inherits(a, "Assay5")) {
        info$count_layers <- .assay5_layer_matches(info$layers, counts_layer)
        info$data_layers <- .assay5_layer_matches(info$layers, data_layer)
      }
    }
  } else if (format == "sce") {
    .require_namespace("SingleCellExperiment")
    .require_namespace("SummarizedExperiment")
    obj <- readRDS(input)
    info$cells <- ncol(obj)
    info$features <- nrow(obj)
    info$assays <- SummarizedExperiment::assayNames(obj)
    info$reductions <- SingleCellExperiment::reducedDimNames(obj)
    info$metadata_columns <- colnames(SummarizedExperiment::colData(obj))
  } else {
    adata <- if (format == "anndata") .read_anndata_file(input) else .read_loom_file(input)
    info$cells <- as.integer(adata$n_obs)
    info$features <- as.integer(adata$n_vars)
    info$layers <- .python_keys(adata$layers)
    info$reductions <- .python_keys(adata$obsm)
    info$metadata_columns <- colnames(as.data.frame(adata$obs))
  }

  message("Cells: ", base::format(info$cells, big.mark = ","))
  message("Features: ", base::format(info$features, big.mark = ","))
  if (!is.null(info$assays)) message("Assays: ", .collapse_or_none(info$assays))
  if (!is.null(info$layers)) message("Layers/slots: ", .collapse_or_none(info$layers))
  if (!is.null(info$count_layers)) message("Detected count layers: ", .collapse_or_none(info$count_layers))
  if (!is.null(info$data_layers)) message("Detected normalized-data layers: ", .collapse_or_none(info$data_layers))
  if (!is.null(info$reductions)) message("Reductions: ", .collapse_or_none(info$reductions))
  invisible(info)
}

#' Inspect a Seurat RDS file
#'
#' Backward-compatible wrapper around [inspect_sc()] for Seurat inputs.
#'
#' @param path_to_file Path to a Seurat `.rds` file.
#' @return Invisibly returns the inspection result.
#' @export
check_input_info <- function(path_to_file) {
  warning(
    "`check_input_info()` is retained for compatibility; use `inspect_sc()` for new code.",
    call. = FALSE
  )
  inspect_sc(path_to_file, format = "seurat")
}
