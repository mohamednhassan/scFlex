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
  ## Access assays available
  assays <- SummarizedExperiment::assayNames(obj)
  ## Access reductions
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
#' Converts an AnnData H5AD file to a SingleCellExperiment object while
#' preserving counts, normalized expression when available, cell metadata,
#' feature metadata, and dimensional reductions.
#'
#' An explicit `layers["counts"]` matrix is preferred as the counts assay.
#' Fractional values are accepted in an explicit counts layer as long as all
#' values are finite and non-negative. If no explicit counts layer is present,
#' `X` is used as counts only when its values are finite, non-negative, and
#' integer-like.
#'
#' Normalized expression is optional. If `layers["logcounts"]` is present it
#' is stored as the SCE `logcounts` assay. Otherwise, `X` is used as
#' `logcounts` when it was not already used as the counts matrix.
#'
#' @param input Path to an AnnData H5AD file.
#' @param output Path where the converted SingleCellExperiment RDS file
#'   should be written.
#'
#' @return Invisibly returns the output file path.
#'
#' @export
convert_anndata_to_sce <- function(input, output) {
  ## Check that the input file exists
  if (!file.exists(input)) {
    stop(
      "File does not exist. Please recheck the path.",
      call. = FALSE
    )
  }

  ## Normalize input path
  input <- normalizePath(input, mustWork = TRUE)

  message("Reading AnnData object...")

  ## Import required Python modules
  reticulate::py_require("anndata>=0.10")

  ad <- reticulate::import(
    "anndata",
    convert = FALSE
  )

  scipy_sparse <- reticulate::import(
    "scipy.sparse",
    convert = FALSE
  )

  ## Read AnnData object
  adata <- ad$read_h5ad(input)

  message("Inspecting AnnData structure...")

  ## Get object dimensions
  n_cells <- reticulate::py_to_r(adata$n_obs)
  n_features <- reticulate::py_to_r(adata$n_vars)

  message("Number of cells: ", n_cells)
  message("Number of features: ", n_features)

  ## Get available layers
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

  ## Get available reductions
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

  ## Get cell names
  cell_names <- reticulate::py_to_r(adata$obs_names$to_list())
  cell_names <- unlist(cell_names, use.names = FALSE)
  cell_names <- as.character(cell_names)

  ## Get feature names
  feature_names <- reticulate::py_to_r(adata$var_names$to_list())
  feature_names <- unlist(feature_names, use.names = FALSE)
  feature_names <- as.character(feature_names)

  ## Check for duplicated names
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

  ## Find raw counts
  raw_counts <- NULL
  counts_source <- NULL

  ## Prefer an explicit counts layer
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

    ## Explicit counts may contain fractional values
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

  ## Fall back to X only when it looks like raw counts
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

  ## Stop if no valid counts matrix was found
  if (is.null(raw_counts)) {
    stop(
      "Raw counts-like matrix not found in X or layers['counts'].\n",
      "Layers present: ",
      collapse_or_none(layers),
      "\nConversion cannot continue because no valid raw counts matrix was found.",
      call. = FALSE
    )
  }

  ## Find normalized expression data
  logcounts <- NULL

  if ("logcounts" %in% layers) {
    logcounts <- reticulate::py_to_r(
      adata$layers$`__getitem__`("logcounts")$astype("float64")
    )

    message("Using layers['logcounts'] as SCE logcounts.")

  } else if (
    !identical(counts_source, "X") &&
    !is.null(adata$X)
  ) {
    logcounts <- reticulate::py_to_r(
      adata$X$astype("float64")
    )

    message("Using X as SCE logcounts.")

  } else {
    message(
      "Normalized expression not found. ",
      "SCE will contain counts only."
    )
  }

  ## AnnData stores cells x features
  ## SCE stores features x cells
  raw_counts <- Matrix::t(raw_counts)

  if (!is.null(logcounts)) {
    logcounts <- Matrix::t(logcounts)
  }

  ## Get cell metadata
  metadata <- reticulate::py_to_r(
    adata$obs
  ) |>
    as.data.frame()

  ## Get feature metadata
  feature_metadata <- reticulate::py_to_r(
    adata$var
  ) |>
    as.data.frame()

  ## Restore matrix and metadata names
  rownames(raw_counts) <- feature_names
  colnames(raw_counts) <- cell_names

  if (!is.null(logcounts)) {
    rownames(logcounts) <- feature_names
    colnames(logcounts) <- cell_names
  }

  rownames(metadata) <- cell_names
  rownames(feature_metadata) <- feature_names

  ## Check counts and cell metadata alignment
  if (!identical(
    colnames(raw_counts),
    rownames(metadata)
  )) {
    stop(
      "Cell metadata is not aligned with the counts matrix.",
      call. = FALSE
    )
  }

  ## Check counts and feature metadata alignment
  if (!identical(
    rownames(raw_counts),
    rownames(feature_metadata)
  )) {
    stop(
      "Feature metadata is not aligned with the counts matrix.",
      call. = FALSE
    )
  }

  ## Check normalized expression alignment when present
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

  ## Build assay list
  assays_list <- list(
    counts = raw_counts
  )

  if (!is.null(logcounts)) {
    assays_list$logcounts <- logcounts
  }

  ## Create SingleCellExperiment object
  sce <- SingleCellExperiment::SingleCellExperiment(
    assays = assays_list,
    colData = S4Vectors::DataFrame(metadata),
    rowData = S4Vectors::DataFrame(feature_metadata)
  )

  ## Restore AnnData reductions
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

      ## Remove the conventional AnnData X_ prefix
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

  ## Prepare output path
  output <- path.expand(output)

  if (!dir.exists(dirname(output))) {
    stop(
      "Output directory does not exist: ",
      dirname(output),
      call. = FALSE
    )
  }

  ## Write output
  saveRDS(
    sce,
    file = output
  )

  message("RDS created successfully.")

  invisible(output)
}
