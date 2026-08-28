.extract_seurat <- function(obj, assay = "RNA", counts_layer = "counts", data_layer = "data", include_reductions = TRUE) {
  .require_namespace("SeuratObject")
  if (!inherits(obj, "Seurat")) stop("Expected a Seurat object.", call. = FALSE)
  assays <- SeuratObject::Assays(obj)
  if (!assay %in% assays) {
    stop("Assay '", assay, "' was not found. Available assays: ", .collapse_or_none(assays), ".", call. = FALSE)
  }

  a <- obj[[assay]]
  counts <- data <- NULL
  notes <- character()
  assay_cells <- colnames(a)
  assay_features <- rownames(a)

  if (inherits(a, "Assay5")) {
    counts_info <- .read_assay5_layer(
      a,
      layer = counts_layer,
      cell_order = assay_cells,
      feature_order = assay_features
    )
    data_info <- .read_assay5_layer(
      a,
      layer = data_layer,
      cell_order = assay_cells,
      feature_order = assay_features
    )

    counts <- counts_info$matrix
    data <- data_info$matrix

    if (isTRUE(counts_info$joined)) {
      notes <- c(
        notes,
        paste0(
          "Joined split Seurat count layers: ",
          paste(counts_info$layers, collapse = ", "),
          "."
        )
      )
    }
    if (isTRUE(data_info$joined)) {
      notes <- c(
        notes,
        paste0(
          "Joined split Seurat normalized-data layers: ",
          paste(data_info$layers, collapse = ", "),
          "."
        )
      )
    }
  } else if (inherits(a, "Assay")) {
    counts <- tryCatch(SeuratObject::GetAssayData(a, slot = counts_layer), error = function(e) NULL)
    data <- tryCatch(SeuratObject::GetAssayData(a, slot = data_layer), error = function(e) NULL)
    if (!is.null(counts) && length(counts) == 0L) counts <- NULL
    if (!is.null(data) && length(data) == 0L) data <- NULL
  } else {
    stop("Unsupported Seurat assay class: ", paste(class(a), collapse = ", "), ".", call. = FALSE)
  }

  reference <- if (!is.null(data)) data else counts
  if (is.null(reference)) {
    stop(
      "No expression matrix was found in Seurat assay '", assay, "'. ",
      "Available layers/slots: ",
      if (inherits(a, "Assay5")) .collapse_or_none(SeuratObject::Layers(a)) else .collapse_or_none(methods::slotNames(a)),
      ".",
      call. = FALSE
    )
  }

  cell_names <- colnames(reference)
  feature_names <- rownames(reference)

  cell_md <- .sanitize_dataframe(obj[[]], "Cell metadata")
  if (!all(cell_names %in% rownames(cell_md))) {
    stop("Selected assay cells could not be matched to Seurat cell metadata.", call. = FALSE)
  }
  cell_md <- cell_md[cell_names, , drop = FALSE]

  feature_md <- tryCatch(
    .sanitize_dataframe(a[[]], "Feature metadata"),
    error = function(e) data.frame(row.names = feature_names)
  )
  if (all(feature_names %in% rownames(feature_md))) {
    feature_md <- feature_md[feature_names, , drop = FALSE]
  } else {
    feature_md <- data.frame(row.names = feature_names)
    notes <- c(notes, "Feature metadata could not be aligned to the selected expression matrix and were omitted.")
  }

  reductions <- list()
  if (isTRUE(include_reductions)) {
    for (red in SeuratObject::Reductions(obj)) {
      emb <- SeuratObject::Embeddings(obj[[red]])
      if (!all(cell_names %in% rownames(emb))) {
        notes <- c(notes, paste0("Reduction '", red, "' did not contain all selected assay cells and was skipped."))
        next
      }
      reductions[[red]] <- as.matrix(emb[cell_names, , drop = FALSE])
    }
  }

  bundle <- .new_bundle(
    counts = counts,
    data = data,
    cell_metadata = cell_md,
    feature_metadata = feature_md,
    reductions = reductions,
    source_format = "seurat",
    source_assay = assay,
    notes = notes
  )
  .validate_bundle(bundle)
  bundle
}

.extract_sce <- function(obj, counts_assay = "counts", data_assay = "logcounts", include_reductions = TRUE) {
  .require_namespace("SingleCellExperiment", "for SingleCellExperiment conversion")
  .require_namespace("SummarizedExperiment", "for SingleCellExperiment conversion")
  if (!inherits(obj, "SingleCellExperiment")) stop("Expected a SingleCellExperiment object.", call. = FALSE)

  available <- SummarizedExperiment::assayNames(obj)
  counts <- if (counts_assay %in% available) SummarizedExperiment::assay(obj, counts_assay) else NULL
  data <- if (data_assay %in% available) SummarizedExperiment::assay(obj, data_assay) else NULL

  if (is.null(counts) && is.null(data) && length(available) > 0L) {
    data <- SummarizedExperiment::assay(obj, available[[1L]])
  }

  cell_md <- .sanitize_dataframe(as.data.frame(SummarizedExperiment::colData(obj)), "Cell metadata")
  rownames(cell_md) <- colnames(obj)
  feature_md <- .sanitize_dataframe(as.data.frame(SummarizedExperiment::rowData(obj)), "Feature metadata")
  rownames(feature_md) <- rownames(obj)

  reductions <- list()
  if (isTRUE(include_reductions)) {
    for (red in SingleCellExperiment::reducedDimNames(obj)) {
      emb <- as.matrix(SingleCellExperiment::reducedDim(obj, red))
      rownames(emb) <- colnames(obj)
      reductions[[red]] <- emb
    }
  }

  bundle <- .new_bundle(
    counts = counts,
    data = data,
    cell_metadata = cell_md,
    feature_metadata = feature_md,
    reductions = reductions,
    source_format = "sce"
  )
  .validate_bundle(bundle)
  bundle
}

.extract_anndata <- function(adata, counts_layer = "counts", include_reductions = TRUE, source_format = "anndata") {
  layers <- .python_keys(adata$layers)
  cell_names <- as.character(adata$obs_names$to_list())
  feature_names <- as.character(adata$var_names$to_list())

  data <- NULL
  if (!is.null(adata$X)) {
    data <- Matrix::t(adata$X)
    rownames(data) <- feature_names
    colnames(data) <- cell_names
  }

  counts <- NULL
  if (counts_layer %in% layers) {
    counts <- Matrix::t(adata$layers[[counts_layer]])
    rownames(counts) <- feature_names
    colnames(counts) <- cell_names
  }

  cell_md <- .sanitize_dataframe(as.data.frame(adata$obs), "Cell metadata")
  rownames(cell_md) <- cell_names
  feature_md <- .sanitize_dataframe(as.data.frame(adata$var), "Feature metadata")
  rownames(feature_md) <- feature_names

  reductions <- list()
  if (isTRUE(include_reductions)) {
    for (red in .python_keys(adata$obsm)) {
      emb <- tryCatch(as.matrix(adata$obsm[[red]]), error = function(e) NULL)
      if (is.null(emb) || nrow(emb) != length(cell_names)) next
      rownames(emb) <- cell_names
      reductions[[sub("^X_", "", red)]] <- emb
    }
  }

  bundle <- .new_bundle(
    counts = counts,
    data = data,
    cell_metadata = cell_md,
    feature_metadata = feature_md,
    reductions = reductions,
    source_format = source_format
  )
  .validate_bundle(bundle)
  bundle
}
