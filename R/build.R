.build_seurat <- function(bundle, assay = "RNA", strict = TRUE) {
  .require_namespace("SeuratObject")
  .validate_bundle(bundle)

  counts <- bundle$counts
  data <- bundle$data

  if (is.null(counts)) {
    if (isTRUE(strict)) {
      stop("Raw counts were not found. Seurat creation requires counts when `strict = TRUE`.", call. = FALSE)
    }
    warning("Raw counts were not found; using the available expression matrix to initialize Seurat.", call. = FALSE)
    counts <- data
  }

  obj <- SeuratObject::CreateSeuratObject(
    counts = counts,
    meta.data = bundle$cell_metadata,
    assay = assay
  )

  if (!is.null(data)) {
    if (inherits(obj[[assay]], "Assay5")) {
      SeuratObject::LayerData(obj[[assay]], layer = "data") <- data
    } else {
      obj <- SeuratObject::SetAssayData(obj, assay = assay, slot = "data", new.data = data)
    }
  }

  if (!is.null(bundle$feature_metadata) && ncol(bundle$feature_metadata) > 0L) {
    try({
      obj[[assay]] <- SeuratObject::AddMetaData(obj[[assay]], metadata = bundle$feature_metadata)
    }, silent = TRUE)
  }

  for (red in names(bundle$reductions)) {
    emb <- bundle$reductions[[red]]
    rownames(emb) <- colnames(obj)
    key <- .make_reduction_key(red)
    colnames(emb) <- paste0(key, seq_len(ncol(emb)))
    obj[[red]] <- SeuratObject::CreateDimReducObject(
      embeddings = emb,
      assay = assay,
      key = key
    )
  }

  obj
}

.build_sce <- function(bundle, counts_assay = "counts", data_assay = "logcounts") {
  .require_namespace("SingleCellExperiment", "for SingleCellExperiment conversion")
  .require_namespace("S4Vectors", "for SingleCellExperiment conversion")
  .validate_bundle(bundle)

  assays <- list()
  if (!is.null(bundle$counts)) assays[[counts_assay]] <- bundle$counts
  if (!is.null(bundle$data)) assays[[data_assay]] <- bundle$data
  if (length(assays) == 0L) stop("No expression matrix is available to construct a SingleCellExperiment.", call. = FALSE)

  cell_md <- if (is.null(bundle$cell_metadata)) data.frame(row.names = colnames(assays[[1L]])) else bundle$cell_metadata
  feature_md <- if (is.null(bundle$feature_metadata)) data.frame(row.names = rownames(assays[[1L]])) else bundle$feature_metadata

  sce <- SingleCellExperiment::SingleCellExperiment(
    assays = assays,
    colData = S4Vectors::DataFrame(cell_md),
    rowData = S4Vectors::DataFrame(feature_md)
  )

  for (red in names(bundle$reductions)) {
    SingleCellExperiment::reducedDim(sce, red) <- bundle$reductions[[red]]
  }

  sce
}

.build_anndata <- function(bundle, counts_layer = "counts", strict = TRUE, include_reductions = TRUE) {
  .validate_bundle(bundle)
  ad <- .import_anndata()

  X_source <- bundle$data
  notes <- character()
  if (is.null(X_source)) {
    X_source <- bundle$counts
    notes <- c(
      notes,
      "Normalized data were unavailable; raw counts were written to AnnData X without normalization."
    )
  }

  X <- Matrix::t(X_source)
  obs <- if (is.null(bundle$cell_metadata)) data.frame(row.names = colnames(X_source)) else bundle$cell_metadata
  var <- if (is.null(bundle$feature_metadata)) data.frame(row.names = rownames(X_source)) else bundle$feature_metadata

  adata <- ad$AnnData(X = X, obs = obs, var = var)

  if (!is.null(bundle$counts)) {
    adata$layers$`__setitem__`(counts_layer, reticulate::r_to_py(Matrix::t(bundle$counts)))
  }

  if (isTRUE(include_reductions)) {
    for (red in names(bundle$reductions)) {
      adata$obsm$`__setitem__`(
        paste0("X_", red),
        reticulate::r_to_py(unname(bundle$reductions[[red]]))
      )
    }
  }

  list(adata = adata, notes = notes)
}
