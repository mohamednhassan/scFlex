#' Convert SingleCellExperiment to AnnData
#'
#' Reads a SingleCellExperiment RDS object and writes an H5AD AnnData file, mapping supported assays, metadata, feature metadata, and reduced dimensions.
#'
#' @param input Path to an input SingleCellExperiment `.rds` file.
#' @param output Path to the output `.h5ad` file.
#'
#' @return The output path, returned invisibly.
#' @export
convert_sce_to_anndata <- function(input, output) {
  if (!file.exists(input)) {
    stop("File does not exist. Please recheck the path.")
  }
  input <- normalizePath(input, mustWork = TRUE)
  message("Reading SingleCellExperiment object...")
  obj <- readRDS(input)
  if (!inherits(obj, "SingleCellExperiment")) {
    stop("Input RDS is not a SingleCellExperiment object.")
  }
  message("Inspecting object structure...")
  assays <- SummarizedExperiment::assayNames(obj)
  reductions <- SingleCellExperiment::reducedDimNames(obj)
  if (!"counts" %in% assays) {
    stop("counts assay not found in SingleCellExperiment.")
  }
  counts_mat <- SummarizedExperiment::assay(obj, "counts")
  data_mat <- if ("logcounts" %in% assays) {
    SummarizedExperiment::assay(obj, "logcounts")
  } else {
    NULL
  }
  message("Counts: ", nrow(counts_mat), " features x ", ncol(counts_mat), " cells")
  if (!is.null(data_mat)) {
    message("Logcounts: ", nrow(data_mat), " features x ", ncol(data_mat), " cells")
  } else {
    message("Logcounts: None")
  }
  cell_names <- colnames(obj)
  feature_names <- rownames(obj)
  obs <- as.data.frame(SummarizedExperiment::colData(obj))
  var <- as.data.frame(SummarizedExperiment::rowData(obj))
  rownames(obs) <- cell_names
  rownames(var) <- feature_names
  obs <- sanitize_obs(obs)
  counts <- Matrix::t(counts_mat)
  if (!is.null(data_mat)) {
    logcounts <- Matrix::t(data_mat)
    X <- logcounts
  } else {
    logcounts <- NULL
    X <- counts
    message("Normalized data unavailable. Using raw counts in AnnData X.")
  }
  stopifnot(
    identical(rownames(obs), rownames(X)),
    identical(rownames(var), colnames(X)),
    identical(rownames(obs), rownames(counts)),
    identical(rownames(var), colnames(counts))
  )
  stopifnot(
    !anyDuplicated(cell_names),
    !anyDuplicated(feature_names)
  )
  message("Matrix dimensions alignment checks passed!")
  reticulate::py_require("anndata>=0.10")
  ad <- reticulate::import("anndata")
  adata <- ad$AnnData(
    X = X,
    obs = obs,
    var = var
  )
  adata$layers["counts"] <- counts
  if (!is.null(logcounts)) {
    adata$layers["logcounts"] <- logcounts
  }
  if (length(reductions) > 0) {
    for (red in reductions) {
      emb <- SingleCellExperiment::reducedDim(obj, red)
      key <- paste0("X_", red)
      adata$obsm$`__setitem__`(
        key,
        reticulate::r_to_py(unname(as.matrix(emb)))
      )
      message(
        red, ": ",
        nrow(emb), " cells x ",
        ncol(emb), " dimensions"
      )
    }
  }
  message("AnnData object created")
  message("Writing h5ad output file...")
  output <- path.expand(output)
  if (!dir.exists(dirname(output))) {
    stop("Output directory does not exist: ", dirname(output))
  }
  adata$write_h5ad(output, convert_strings_to_categoricals = FALSE)
  message("H5AD created successfully.")
  invisible(output)
}

#' Convert AnnData to SingleCellExperiment
#'
#' Reads an H5AD AnnData file and creates a SingleCellExperiment while preserving supported count data, normalized data, metadata, feature metadata, and reductions.
#'
#' @details
#' Raw counts are selected conservatively from an explicit `counts` layer when available or from `X` when it is non-negative and integer-like. The function does not guess raw counts from arbitrary named layers.
#'
#' @param input Path to an input `.h5ad` file.
#' @param output Path to the output SingleCellExperiment `.rds` file.
#'
#' @return The output path, returned invisibly on success. If suitable raw count data cannot be identified, the function may return `NULL` invisibly without writing an output object.
#' @export
convert_anndata_to_sce <- function(input, output) {
  if (!file.exists(input)) {
    stop("File does not exist. Please recheck the path.")
  }
  input <- normalizePath(input, mustWork = TRUE)
  message("Reading AnnData object...")
  reticulate::py_require("anndata>=0.10")
  ad <- reticulate::import("anndata", convert = FALSE)
  scipy_sparse <- reticulate::import("scipy.sparse", convert = FALSE)
  adata <- ad$read_h5ad(input)
  message("Inspecting AnnData structure...")
  n_cells <- reticulate::py_to_r(adata$n_obs)
  n_features <- reticulate::py_to_r(adata$n_vars)
  message("Number of cells: ", n_cells)
  message("Number of features: ", n_features)
  layers <- reticulate::iterate(
    adata$layers$keys(),
    as.character
  ) |>
    unlist(use.names = FALSE)
  layers <- layers[!is.na(layers) & layers != "" & layers != "None"]
  reductions <- reticulate::iterate(
    adata$obsm$keys(),
    as.character
  ) |>
    unlist(use.names = FALSE)
  reductions <- reductions[!is.na(reductions) & reductions != "" & reductions != "None"]
  cell_names <- reticulate::iterate(
    adata$obs_names,
    as.character
  ) |>
    unlist(use.names = FALSE)
  feature_names <- reticulate::iterate(
    adata$var_names,
    as.character
  ) |>
    unlist(use.names = FALSE)
  raw_counts <- NULL
  counts_source <- NULL
  if ("counts" %in% layers) {
    mat <- adata$layers$`__getitem__`("counts")
    if (reticulate::py_to_r(scipy_sparse$issparse(mat))) {
      values <- reticulate::py_to_r(mat$data$astype("float64"))
    } else {
      values <- as.numeric(
        reticulate::py_to_r(mat$astype("float64"))
      )
    }
    if (
      all(values >= 0) &&
      all(abs(values - round(values)) < 1e-8)
    ) {
      raw_counts <- reticulate::py_to_r(
        mat$astype("float64")
      )
      counts_source <- "layers['counts']"
    }
  }
  if (is.null(raw_counts)) {
    if (reticulate::py_to_r(scipy_sparse$issparse(adata$X))) {
      values <- reticulate::py_to_r(
        adata$X$data$astype("float64")
      )
    } else {
      values <- as.numeric(
        reticulate::py_to_r(
          adata$X$astype("float64")
        )
      )
    }
    if (
      all(values >= 0) &&
      all(abs(values - round(values)) < 1e-8)
    ) {
      raw_counts <- reticulate::py_to_r(
        adata$X$astype("float64")
      )
      counts_source <- "X"
    }
  }
  if (is.null(raw_counts)) {
    message(
      "Raw counts-like matrix not found in X or layers['counts'].\n",
      "Layers present: ",
      if (length(layers) > 0) paste(layers, collapse = ", ") else "None",
      "\nConversion stopped."
    )
    return(invisible(NULL))
  }
  message("Using ", counts_source, " as SCE counts.")
  raw_counts <- Matrix::t(raw_counts)
  rownames(raw_counts) <- feature_names
  colnames(raw_counts) <- cell_names
  logcounts <- NULL
  if ("logcounts" %in% layers) {
    mat <- adata$layers$`__getitem__`("logcounts")
    logcounts <- reticulate::py_to_r(
      mat$astype("float64")
    )
    logcounts <- Matrix::t(logcounts)
    rownames(logcounts) <- feature_names
    colnames(logcounts) <- cell_names
    message("Using layers['logcounts'] as SCE logcounts.")
  } else if (counts_source != "X") {
    logcounts <- reticulate::py_to_r(
      adata$X$astype("float64")
    )
    logcounts <- Matrix::t(logcounts)
    rownames(logcounts) <- feature_names
    colnames(logcounts) <- cell_names
    message("Using X as SCE logcounts.")
  } else {
    message("Logcounts not found. SCE will contain counts only.")
  }
  metadata <- reticulate::py_to_r(adata$obs) |>
    as.data.frame()
  feature_metadata <- reticulate::py_to_r(adata$var) |>
    as.data.frame()
  rownames(metadata) <- cell_names
  rownames(feature_metadata) <- feature_names
  stopifnot(
    identical(colnames(raw_counts), rownames(metadata)),
    identical(rownames(raw_counts), rownames(feature_metadata))
  )
  if (!is.null(logcounts)) {
    stopifnot(
      identical(dim(raw_counts), dim(logcounts)),
      identical(rownames(raw_counts), rownames(logcounts)),
      identical(colnames(raw_counts), colnames(logcounts))
    )
  }
  stopifnot(
    !anyDuplicated(cell_names),
    !anyDuplicated(feature_names)
  )
  message("Matrix alignment checks passed!")
  assays_list <- list(
    counts = raw_counts
  )
  if (!is.null(logcounts)) {
    assays_list$logcounts <- logcounts
  }
  sce <- SingleCellExperiment::SingleCellExperiment(
    assays = assays_list,
    colData = S4Vectors::DataFrame(metadata),
    rowData = S4Vectors::DataFrame(feature_metadata)
  )
  if (length(reductions) > 0) {
    for (red in reductions) {
      emb <- reticulate::py_to_r(
        adata$obsm$`__getitem__`(red)
      ) |>
        as.matrix()
      rownames(emb) <- cell_names
      sce_red_name <- sub("^X_", "", red)
      SingleCellExperiment::reducedDim(
        sce,
        sce_red_name
      ) <- emb
      message(
        sce_red_name, ": ",
        nrow(emb), " cells x ",
        ncol(emb), " dimensions"
      )
    }
  }
  message("SingleCellExperiment object created")
  message("Writing RDS output file...")
  output <- path.expand(output)
  if (!dir.exists(dirname(output))) {
    stop("Output directory does not exist: ", dirname(output))
  }
  saveRDS(sce, file = output)
  message("RDS created successfully.")
  invisible(output)
}
