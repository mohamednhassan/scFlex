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

  layers <- layers[
    !is.na(layers) &
      layers != "" &
      layers != "None"
  ]

  reductions <- reticulate::iterate(
    adata$obsm$keys(),
    as.character
  ) |>
    unlist(use.names = FALSE)

  reductions <- reductions[
    !is.na(reductions) &
      reductions != "" &
      reductions != "None"
  ]

  # Get axis names
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

  if (anyDuplicated(cell_names)) {
    stop(
      "Duplicated cell names detected in AnnData.",
      call. = FALSE
    )
  }

  if (anyDuplicated(feature_names)) {
    stop(
      "Duplicated feature names detected in AnnData.",
      call. = FALSE
    )
  }

  # Find raw counts
  raw_counts <- NULL
  counts_source <- NULL

  # Prefer an explicit counts layer.
  # Explicit counts are allowed to contain fractional values,
  # for example after denoising or background correction.
  if ("counts" %in% layers) {
    mat <- adata$layers$`__getitem__`("counts")

    if (reticulate::py_to_r(scipy_sparse$issparse(mat))) {
      values <- reticulate::py_to_r(
        mat$data$astype("float64")
      )
    } else {
      values <- as.numeric(
        reticulate::py_to_r(
          mat$astype("float64")
        )
      )
    }

    if (
      all(is.finite(values)) &&
      all(values >= 0)
    ) {
      raw_counts <- reticulate::py_to_r(
        mat$astype("float64")
      )

      counts_source <- "layers['counts']"

      message("Using layers['counts'] as SCE counts.")
    } else {
      stop(
        "AnnData layers['counts'] contains negative or non-finite values.",
        call. = FALSE
      )
    }
  }

  # If no explicit counts layer exists, X is only accepted as
  # counts when it is finite, non-negative and integer-like.
  if (is.null(raw_counts)) {
    mat <- adata$X

    if (!is.null(mat)) {
      if (reticulate::py_to_r(scipy_sparse$issparse(mat))) {
        values <- reticulate::py_to_r(
          mat$data$astype("float64")
        )
      } else {
        values <- as.numeric(
          reticulate::py_to_r(
            mat$astype("float64")
          )
        )
      }

      if (
        all(is.finite(values)) &&
        all(values >= 0) &&
        all(abs(values - round(values)) < 1e-8)
      ) {
        raw_counts <- reticulate::py_to_r(
          mat$astype("float64")
        )

        counts_source <- "X"

        message("Using X as SCE counts.")
      }
    }
  }

  # Conversion requires a valid counts matrix.
  # Do not silently infer counts from raw.X or arbitrary layers.
  if (is.null(raw_counts)) {
    message(
      "Raw counts-like matrix not found in X or layers['counts']."
    )

    message(
      "Layers present: ",
      collapse_or_none(layers)
    )

    stop(
      "Conversion cannot continue because no valid raw counts matrix was found.",
      call. = FALSE
    )
  }

  # Find normalized expression data
  logcounts <- NULL

  if ("logcounts" %in% layers) {
    logcounts <- reticulate::py_to_r(
      adata$layers$`__getitem__`("logcounts")$astype("float64")
    )

    message("Using layers['logcounts'] as SCE logcounts.")
  } else if (!identical(counts_source, "X") && !is.null(adata$X)) {
    logcounts <- reticulate::py_to_r(
      adata$X$astype("float64")
    )

    message("Using X as SCE logcounts.")
  } else {
    message(
      "Logcounts not found. SCE will contain counts only."
    )
  }

  # AnnData = cells x features
  # SCE     = features x cells
  raw_counts <- Matrix::t(raw_counts)

  if (!is.null(logcounts)) {
    logcounts <- Matrix::t(logcounts)
  }

  # Cell metadata
  metadata <- reticulate::py_to_r(
    adata$obs
  ) |>
    as.data.frame()

  # Feature metadata
  feature_metadata <- reticulate::py_to_r(
    adata$var
  ) |>
    as.data.frame()

  # Restore dimnames
  rownames(raw_counts) <- feature_names
  colnames(raw_counts) <- cell_names

  if (!is.null(logcounts)) {
    rownames(logcounts) <- feature_names
    colnames(logcounts) <- cell_names
  }

  rownames(metadata) <- cell_names
  rownames(feature_metadata) <- feature_names

  # Alignment checks
  if (!identical(
    colnames(raw_counts),
    rownames(metadata)
  )) {
    stop(
      "Cell metadata is not aligned with the counts matrix.",
      call. = FALSE
    )
  }

  if (!identical(
    rownames(raw_counts),
    rownames(feature_metadata)
  )) {
    stop(
      "Feature metadata is not aligned with the counts matrix.",
      call. = FALSE
    )
  }

  if (!is.null(logcounts)) {
    if (
      !identical(dim(raw_counts), dim(logcounts)) ||
      !identical(rownames(raw_counts), rownames(logcounts)) ||
      !identical(colnames(raw_counts), colnames(logcounts))
    ) {
      stop(
        "Logcounts are not aligned with the counts matrix.",
        call. = FALSE
      )
    }
  }

  message("Matrix alignment checks passed!")

  # Build assay list
  assays_list <- list(
    counts = raw_counts
  )

  if (!is.null(logcounts)) {
    assays_list$logcounts <- logcounts
  }

  # Construct SingleCellExperiment
  sce <- SingleCellExperiment::SingleCellExperiment(
    assays = assays_list,
    colData = S4Vectors::DataFrame(metadata),
    rowData = S4Vectors::DataFrame(feature_metadata)
  )

  # Restore AnnData reductions
  if (length(reductions) > 0) {
    for (red in reductions) {
      emb <- reticulate::py_to_r(
        adata$obsm$`__getitem__`(red)
      ) |>
        as.matrix()

      rownames(emb) <- cell_names

      if (!identical(
        rownames(emb),
        colnames(sce)
      )) {
        stop(
          "Cell names in embedding '",
          red,
          "' do not match the SCE cell order.",
          call. = FALSE
        )
      }

      sce_red_name <- sub(
        "^X_",
        "",
        red
      )

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
    stop(
      "Output directory does not exist: ",
      dirname(output),
      call. = FALSE
    )
  }

  saveRDS(
    sce,
    file = output
  )

  message("RDS created successfully.")

  invisible(output)
}
