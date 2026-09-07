#' Convert a Seurat Assay5 object to a classic Assay
#'
#' Converts the default Seurat Assay5 assay to a classic Seurat Assay, adds it to the object as `RNA_classic`, makes it the default assay, and saves the updated object as RDS.
#'
#' @param input An RDS file containing a Seurat object whose default assay is an Assay5 object.
#' @param output Path to the output RDS file.
#'
#' @return The output path, returned invisibly.
#' @export
convert_seu_v5_to_classic <- function(input, output) {
  # Checking input
  info <- suppressMessages(inspect_sc(path_to_file = input))
  if (info$rna_structure != "v5") {
    stop("RNA assay is not Assay5. No conversion performed.")
  }
  message(
    "\nSource RNA assay is Assay5. Starting conversion to classic Assay..."
  )
  # Reading object
  obj <- readRDS(input)
  
  # Casting an Assay5 to a classic Assay
  source_assay <- SeuratObject::DefaultAssay(obj)
  rna_classic <- methods::as(obj[[source_assay]], Class = "Assay")
  SeuratObject::Key(rna_classic) <- "classicRNA_"
  obj[["RNA_classic"]] <- rna_classic
  SeuratObject::DefaultAssay(obj) <- "RNA_classic"
  message("\nConversion completed. RNA_classic assay has been added to the Seurat object")
  # Saving RDS
  message("\nSaving object...\n")
  saveRDS(obj, output)
  message("Object has been saved!\n")
  
  invisible(output)
}

#' Convert a classic Seurat Assay to Assay5
#'
#' Converts the default classic Seurat Assay to Assay5, adds it to the object as `RNA_assay5`, makes it the default assay, and saves the updated object as RDS.
#'
#' @param input An RDS file containing a Seurat object whose default assay is a classic Assay object.
#' @param output Path to the output RDS file.
#'
#' @return The output path, returned invisibly.
#' @export
convert_seu_classic_to_v5 <- function(input, output) {
  # Checking input
  info <- suppressMessages(inspect_sc(path_to_file = input))
  if (info$rna_structure != "classic") {
    stop("RNA assay is not classic. No conversion performed.")
  }
  # Check installed Seurat version
  seurat_version <- utils::packageVersion("Seurat")
  if (utils::compareVersion(as.character(seurat_version), "5.0.0") < 0) {
    stop(
      "Converting a classic Assay to Assay5 requires Seurat >= 5.0.0. ",
      "Current Seurat version: ",
      seurat_version
    )
  }
  message("\nSource RNA assay is classic. Starting conversion to Assay5...")
  
  # Reading object
  obj <- readRDS(input)
  # Casting classic Assay -> Assay5
  source_assay <- SeuratObject::DefaultAssay(obj)
  rna_assay5 <- methods::as(
    obj[[source_assay]],
    Class = "Assay5"
  )
  # Assign unique key
  SeuratObject::Key(rna_assay5) <- "RNAassay5_"
  # Add new assay
  obj[["RNA_assay5"]] <- rna_assay5
  SeuratObject::DefaultAssay(obj) <- "RNA_assay5"
  message("\nConversion completed. RNA_assay5 has been added to the Seurat object")
  # Saving RDS
  message("\nSaving object...")
  saveRDS(obj, output)
  message("Object has been saved!\n")
  invisible(output)
}

#' Convert AnnData to Seurat
#'
#' Reads an H5AD AnnData file and creates a Seurat object while preserving supported expression data, metadata, feature metadata, and dimensional reductions when available.
#'
#' @details
#' Raw counts are selected conservatively from an explicit `counts` layer when available or from `X` when it is non-negative and integer-like. The function does not guess raw counts from arbitrary named layers.
#'
#' @param input Path to an input `.h5ad` file.
#' @param output Path to the output Seurat `.rds` file.
#'
#' @return The output path, returned invisibly on success. If suitable raw count data cannot be identified, the function may return `NULL` invisibly without writing an output object.
#' @export
convert_anndata_to_seurat <- function(input, output) {
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
  )

  reductions <- reticulate::iterate(
    adata$obsm$keys(),
    as.character
  )

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
    stop("Duplicated cell names detected in AnnData.", call. = FALSE)
  }

  if (anyDuplicated(feature_names)) {
    stop("Duplicated feature names detected in AnnData.", call. = FALSE)
  }

  # Find raw counts
  raw_counts <- NULL
  counts_source <- NULL

  # Prefer an explicit counts layer.
  # Explicit counts may contain fractional values, for example after
  # background correction or denoising, so integer values are not required.
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

      message("Using layers['counts'] as Seurat counts.")
    } else {
      stop(
        "AnnData layers['counts'] contains negative or non-finite values.",
        call. = FALSE
      )
    }
  }

  # If there is no explicit counts layer, only use X when it clearly
  # looks like raw counts: finite, non-negative and integer-like.
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

        message("Using X as Seurat counts.")
      }
    }
  }

  # Conversion to Seurat requires a valid counts matrix.
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
  data_mat <- NULL

  if ("logcounts" %in% layers) {
    data_mat <- reticulate::py_to_r(
      adata$layers$`__getitem__`("logcounts")$astype("float64")
    )

    message("Using layers['logcounts'] as Seurat data.")
  } else if (!identical(counts_source, "X") && !is.null(adata$X)) {
    data_mat <- reticulate::py_to_r(
      adata$X$astype("float64")
    )

    message("Using X as Seurat data.")
  } else {
    message(
      "Normalized data not found. Seurat will contain counts only."
    )
  }

  # AnnData = cells x features
  # Seurat  = features x cells
  raw_counts <- Matrix::t(raw_counts)

  if (!is.null(data_mat)) {
    data_mat <- Matrix::t(data_mat)
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

  message(
    "Raw counts: ",
    paste(dim(raw_counts), collapse = " x ")
  )

  if (!is.null(data_mat)) {
    message(
      "Data: ",
      paste(dim(data_mat), collapse = " x ")
    )
  } else {
    message("Data: None")
  }

  message(
    "Metadata: ",
    paste(dim(metadata), collapse = " x ")
  )

  message(
    "Feature metadata: ",
    paste(dim(feature_metadata), collapse = " x ")
  )

  # Restore R dimnames
  rownames(raw_counts) <- feature_names
  colnames(raw_counts) <- cell_names

  if (!is.null(data_mat)) {
    rownames(data_mat) <- feature_names
    colnames(data_mat) <- cell_names
  }

  rownames(metadata) <- cell_names
  rownames(feature_metadata) <- feature_names

  # Alignment checks
  if (!identical(
    colnames(raw_counts),
    rownames(metadata)
  )) {
    stop(
      "Cell metadata is not aligned with the expression matrix.",
      call. = FALSE
    )
  }

  if (!identical(
    rownames(raw_counts),
    rownames(feature_metadata)
  )) {
    stop(
      "Feature metadata is not aligned with the expression matrix.",
      call. = FALSE
    )
  }

  if (!is.null(data_mat)) {
    if (
      !identical(dim(raw_counts), dim(data_mat)) ||
      !identical(rownames(raw_counts), rownames(data_mat)) ||
      !identical(colnames(raw_counts), colnames(data_mat))
    ) {
      stop(
        "Normalized data is not aligned with the counts matrix.",
        call. = FALSE
      )
    }
  }

  message("Matrix alignment checks passed!")

  # Extract AnnData embeddings
  reduction_list <- list()

  if (length(reductions) > 0) {
    for (red in reductions) {
      emb <- reticulate::py_to_r(
        adata$obsm$`__getitem__`(red)
      ) |>
        as.matrix()

      rownames(emb) <- cell_names

      if (!identical(
        rownames(emb),
        colnames(raw_counts)
      )) {
        stop(
          "Cell names in embedding '",
          red,
          "' do not match the Seurat cell order.",
          call. = FALSE
        )
      }

      reduction_list[[red]] <- emb

      message(
        red, ": ",
        nrow(emb), " cells x ",
        ncol(emb), " dimensions"
      )
    }
  }

  # Construct Seurat object
  message("\nCreating Seurat object...")

  seuratObj <- SeuratObject::CreateSeuratObject(
    counts = raw_counts
  )

  # Restore source cell metadata after object construction.
  # CreateSeuratObject() calculates its own nCount/nFeature columns,
  # but preservation requires the source metadata to win.
  if (!identical(
    rownames(metadata),
    colnames(seuratObj)
  )) {
    stop(
      "Cell metadata is not aligned with the created Seurat object.",
      call. = FALSE
    )
  }

  seuratObj@meta.data <- metadata

  # Restore normalized expression
  if (!is.null(data_mat)) {
    SeuratObject::LayerData(
      seuratObj[["RNA"]],
      layer = "data"
    ) <- data_mat
  }

  # Restore feature metadata
  if (ncol(feature_metadata) > 0) {
    if (!identical(
      rownames(feature_metadata),
      rownames(seuratObj)
    )) {
      stop(
        "Feature metadata is not aligned with the created Seurat object.",
        call. = FALSE
      )
    }

    seuratObj[["RNA"]]@meta.data <- feature_metadata
  }

  message(
    "Counts have been added. Now adding reductions..\n"
  )

  # Restore dimensional reductions
  for (red in names(reduction_list)) {
    seurat_red_name <- sub("^X_", "", red)

    emb <- reduction_list[[red]]

    key <- paste0(
      toupper(
        gsub(
          "[^A-Za-z0-9]",
          "",
          seurat_red_name
        )
      ),
      "_"
    )

    colnames(emb) <- paste0(
      key,
      seq_len(ncol(emb))
    )

    seuratObj[[seurat_red_name]] <-
      SeuratObject::CreateDimReducObject(
        embeddings = emb,
        assay = "RNA",
        key = key
      )
  }

  print(seuratObj)

  message(
    "\nExporting Seurat as RDS file...\n"
  )

  output <- path.expand(output)

  if (!dir.exists(dirname(output))) {
    stop(
      "Output directory does not exist: ",
      dirname(output),
      call. = FALSE
    )
  }

  saveRDS(
    seuratObj,
    file = output
  )

  message("Object created successfully.")

  invisible(output)
}

#' Convert Seurat to AnnData
#'
#' Reads a Seurat RDS object and writes an H5AD AnnData representation while preserving supported expression matrices, metadata, feature metadata, and reductions.
#'
#' @param input Path to an input Seurat `.rds` file.
#' @param output Path to the output `.h5ad` file.
#'
#' @return The output path, returned invisibly.
#' @export
convert_seurat_to_anndata <- function(input, output) {
  if (!file.exists(input)) {
    stop("File does not exist. Please recheck the path.")
  }
  input <- normalizePath(input, mustWork = TRUE)
  message("Reading object...")
  obj <- readRDS(input)
  if (!inherits(obj, "Seurat")) {
    stop("Input RDS is not a Seurat object.")
  }
  message("Inspecting object structure...\n")
  assay_name <- SeuratObject::DefaultAssay(obj)
  assay_obj <- obj[[assay_name]]
  message("Default assay: ", assay_name)
  counts_mat <- NULL
  data_mat <- NULL
  if (inherits(assay_obj, "Assay")) {
    message("Assay structure: classic Assay")
    counts_mat <- SeuratObject::GetAssayData(
      obj,
      assay = assay_name,
      layer = "counts"
    )
    data_mat <- SeuratObject::GetAssayData(
      obj,
      assay = assay_name,
      layer = "data"
    )
    if (length(data_mat) == 0) data_mat <- NULL
  } else if (inherits(assay_obj, "Assay5")) {
    message("Assay structure: Assay5")
    layers <- SeuratObject::Layers(assay_obj)
    if (any(grepl("^counts\\.", layers)) || any(grepl("^data\\.", layers))) {
      assay_obj <- SeuratObject::JoinLayers(assay_obj)
      layers <- SeuratObject::Layers(assay_obj)
    }
    if ("counts" %in% layers) {
      counts_mat <- SeuratObject::LayerData(assay_obj, layer = "counts")
    }
    if ("data" %in% layers) {
      data_mat <- SeuratObject::LayerData(assay_obj, layer = "data")
    }
  } else {
    stop("Unsupported Seurat assay structure.")
  }
  if (is.null(counts_mat)) {
    stop("Counts layer not found in the selected Seurat assay.")
  }
  if (length(counts_mat) == 0) counts_mat <- NULL
  if (!is.null(data_mat) && length(data_mat) == 0) data_mat <- NULL
  if (is.null(counts_mat)) {
    stop("Counts layer is empty in the selected Seurat assay.")
  }
  message("Counts: ", nrow(counts_mat), " features x ", ncol(counts_mat), " cells")
  if (is.null(data_mat)) {
    message("Data: None")
  } else {
    message("Data: ", nrow(data_mat), " features x ", ncol(data_mat), " cells")
  }
  X_source <- if (!is.null(data_mat)) data_mat else counts_mat
  X <- Matrix::t(X_source)
  counts <- Matrix::t(counts_mat)
  obs <- sanitize_obs(obj@meta.data)
  feature_metadata <- as.data.frame(assay_obj[[]])
  var <- feature_metadata
  rownames(var) <- rownames(counts_mat)
  stopifnot(
    identical(rownames(obs), rownames(X)),
    identical(rownames(var), colnames(X)),
    identical(rownames(obs), rownames(counts)),
    identical(rownames(var), colnames(counts))
  )
  stopifnot(
    !anyDuplicated(rownames(obs)),
    !anyDuplicated(rownames(var))
  )
  message("Matrix dimensions alignment checks passed!")
  reductions <- SeuratObject::Reductions(obj)
  reticulate::py_require("anndata>=0.10")
  ad <- reticulate::import("anndata")
  adata <- ad$AnnData(
    X = X,
    obs = obs,
    var = var
  )
  if (length(reductions) > 0) {
    for (red in reductions) {
      key <- paste0("X_", red)
      reduction_matrix <- SeuratObject::Embeddings(obj[[red]])
      reduction_py <- reticulate::r_to_py(unname(reduction_matrix))
      adata$obsm$`__setitem__`(key, reduction_py)
      message(red, ": ", nrow(reduction_matrix), " cells x ", ncol(reduction_matrix), " dimensions")
    }
  }
  adata$layers["counts"] <- counts
  message("AnnData object created")
  message("Writing h5ad output file...\n")
  output <- path.expand(output)
  adata$write_h5ad(output, convert_strings_to_categoricals = FALSE)
  message("H5AD created successfully.\n")
  invisible(output)
}

#' Convert Seurat to SingleCellExperiment
#'
#' Reads a Seurat RDS object and creates a SingleCellExperiment while preserving supported assays, cell metadata, feature metadata, and reductions.
#'
#' @param input Path to an input Seurat `.rds` file.
#' @param output Path to the output SingleCellExperiment `.rds` file.
#'
#' @return The output path, returned invisibly.
#' @export
convert_seurat_to_sce <- function(input, output) {
  if (!file.exists(input)) {
    stop("File does not exist. Please recheck the path.")
  }
  input <- normalizePath(input, mustWork = TRUE)
  message("Reading Seurat object...")
  obj <- readRDS(input)
  if (!inherits(obj, "Seurat")) {
    stop("Input RDS is not a Seurat object.")
  }
  message("Inspecting object structure...")
  assay_name <- SeuratObject::DefaultAssay(obj)
  assay_obj <- obj[[assay_name]]
  message("Default assay: ", assay_name)
  counts_mat <- NULL
  data_mat <- NULL
  if (inherits(assay_obj, "Assay")) {
    message("Assay structure: classic Assay")
    counts_mat <- SeuratObject::GetAssayData(
      obj,
      assay = assay_name,
      layer = "counts"
    )
    data_mat <- SeuratObject::GetAssayData(
      obj,
      assay = assay_name,
      layer = "data"
    )
  } else if (inherits(assay_obj, "Assay5")) {
    message("Assay structure: Assay5")
    layers <- SeuratObject::Layers(assay_obj)
    if ("counts" %in% layers) {
      counts_mat <- SeuratObject::LayerData(
        assay_obj,
        layer = "counts"
      )
    } else if (any(grepl("^counts\\.", layers))) {
      tmp_name <- ".scFlex_counts"
      assay_tmp <- SeuratObject::JoinLayers(
        assay_obj,
        layers = "counts",
        new = tmp_name
      )
      counts_mat <- SeuratObject::LayerData(
        assay_tmp,
        layer = tmp_name
      )
    }
    if ("data" %in% layers) {
      data_mat <- SeuratObject::LayerData(
        assay_obj,
        layer = "data"
      )
    } else if (any(grepl("^data\\.", layers))) {
      tmp_name <- ".scFlex_data"
      assay_tmp <- SeuratObject::JoinLayers(
        assay_obj,
        layers = "data",
        new = tmp_name
      )
      data_mat <- SeuratObject::LayerData(
        assay_tmp,
        layer = tmp_name
      )
    }
  } else {
    stop("Unsupported Seurat assay structure.")
  }
  if (is.null(counts_mat)) {
    stop("Counts could not be found in the Seurat assay.")
  }
  message(
    "Counts: ",
    nrow(counts_mat), " features x ",
    ncol(counts_mat), " cells"
  )
  if (!is.null(data_mat) && nrow(data_mat) > 0 && ncol(data_mat) > 0) {
    message(
      "Data: ",
      nrow(data_mat), " features x ",
      ncol(data_mat), " cells"
    )
  } else {
    data_mat <- NULL
    message("Data: None")
  }
  cell_names <- colnames(counts_mat)
  feature_names <- rownames(counts_mat)
  metadata <- obj@meta.data
  metadata <- metadata[colnames(counts_mat), , drop = FALSE]
  feature_metadata <- as.data.frame(assay_obj[[]])
  if (ncol(feature_metadata) == 0) {
    feature_metadata <- data.frame(row.names = rownames(counts_mat))
  } else {
    feature_metadata <- feature_metadata[rownames(counts_mat), , drop = FALSE]
  }
  stopifnot(
    identical(colnames(counts_mat), rownames(metadata)),
    identical(rownames(counts_mat), rownames(feature_metadata))
  )
  if (!is.null(data_mat)) {
    stopifnot(
      identical(dim(counts_mat), dim(data_mat)),
      identical(rownames(counts_mat), rownames(data_mat)),
      identical(colnames(counts_mat), colnames(data_mat))
    )
  }
  stopifnot(
    !anyDuplicated(cell_names),
    !anyDuplicated(feature_names)
  )
  message("Matrix alignment checks passed!")
  assays_list <- list(
    counts = counts_mat
  )
  if (!is.null(data_mat)) {
    assays_list$logcounts <- data_mat
  }
  sce <- SingleCellExperiment::SingleCellExperiment(
    assays = assays_list,
    colData = S4Vectors::DataFrame(metadata),
    rowData = S4Vectors::DataFrame(feature_metadata)
  )
  reductions <- SeuratObject::Reductions(obj)
  if (length(reductions) > 0) {
    for (red in reductions) {
      emb <- SeuratObject::Embeddings(obj[[red]])
      emb <- emb[colnames(counts_mat), , drop = FALSE]
      SingleCellExperiment::reducedDim(
        sce,
        red
      ) <- emb
      message(
        red, ": ",
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

#' Convert SingleCellExperiment to Seurat
#'
#' Reads a SingleCellExperiment RDS object and creates a Seurat object while preserving supported expression assays, metadata, feature metadata, and reductions.
#'
#' @param input Path to an input SingleCellExperiment `.rds` file.
#' @param output Path to the output Seurat `.rds` file.
#'
#' @return The output path, returned invisibly.
#' @export
convert_sce_to_seurat <- function(input, output) {
  if (!file.exists(input)) {
    stop("File does not exist. Please recheck the path.")
  }

  input <- normalizePath(input, mustWork = TRUE)

  message("Reading SingleCellExperiment object...")

  obj <- suppressPackageStartupMessages(
    readRDS(input)
  )

  if (!inherits(obj, "SingleCellExperiment")) {
    stop("Input RDS is not a SingleCellExperiment object.")
  }

  message("Inspecting object structure...")

  assays <- SummarizedExperiment::assayNames(obj)
  reductions <- SingleCellExperiment::reducedDimNames(obj)

  if (!"counts" %in% assays) {
    stop("counts assay not found in SingleCellExperiment.")
  }

  counts_mat <- SummarizedExperiment::assay(
    obj,
    "counts"
  )

  data_mat <- if ("logcounts" %in% assays) {
    SummarizedExperiment::assay(
      obj,
      "logcounts"
    )
  } else {
    NULL
  }

  message(
    "Counts: ",
    nrow(counts_mat), " features x ",
    ncol(counts_mat), " cells"
  )

  if (!is.null(data_mat)) {
    message(
      "Logcounts: ",
      nrow(data_mat), " features x ",
      ncol(data_mat), " cells"
    )
  } else {
    message("Logcounts: None")
  }

  cell_names <- colnames(obj)
  feature_names <- rownames(obj)

  metadata <- as.data.frame(
    SummarizedExperiment::colData(obj)
  )

  feature_metadata <- as.data.frame(
    SummarizedExperiment::rowData(obj)
  )

  rownames(metadata) <- cell_names
  rownames(feature_metadata) <- feature_names

  stopifnot(
    identical(colnames(counts_mat), rownames(metadata)),
    identical(rownames(counts_mat), rownames(feature_metadata))
  )

  if (!is.null(data_mat)) {
    stopifnot(
      identical(dim(counts_mat), dim(data_mat)),
      identical(rownames(counts_mat), rownames(data_mat)),
      identical(colnames(counts_mat), colnames(data_mat))
    )
  }

  stopifnot(
    !anyDuplicated(cell_names),
    !anyDuplicated(feature_names)
  )

  message("Matrix alignment checks passed!")

  seuratObj <- SeuratObject::CreateSeuratObject(
    counts = counts_mat
  )

  assay_name <- SeuratObject::DefaultAssay(seuratObj)

  if (!is.null(data_mat)) {
    seuratObj <- SeuratObject::SetAssayData(
      seuratObj,
      assay = assay_name,
      layer = "data",
      new.data = data_mat
    )
  }

  if (ncol(feature_metadata) > 0) {
    seuratObj[[assay_name]][[]] <- feature_metadata
  }

  if (length(reductions) > 0) {
    for (red in reductions) {
      emb <- SingleCellExperiment::reducedDim(
        obj,
        red
      )

      emb <- as.matrix(emb)
      emb <- emb[cell_names, , drop = FALSE]

      dim_names <- colnames(emb)
      key <- NULL

      if (!is.null(dim_names) && length(dim_names) > 0) {
        prefixes <- sub("[0-9]+$", "", dim_names)

        if (length(unique(prefixes)) == 1 &&
            nzchar(prefixes[1]) &&
            all(grepl("^[A-Za-z][A-Za-z0-9]*_$", prefixes))) {
          key <- prefixes[1]
        }
      }

      if (is.null(key)) {
        key <- paste0(
          toupper(red),
          "_"
        )
      }

      seuratObj[[red]] <- SeuratObject::CreateDimReducObject(
        embeddings = emb,
        key = key,
        assay = assay_name
      )

      message(
        red, ": ",
        nrow(emb), " cells x ",
        ncol(emb), " dimensions"
      )
    }
  }

  if (!identical(rownames(metadata), colnames(seuratObj))) {
    stop(
      "Cell metadata is not aligned with the created Seurat object.",
      call. = FALSE
    )
  }

  seuratObj@meta.data <- metadata

  message("Seurat object created")
  message("Writing RDS output file...")

  output <- path.expand(output)

  if (!dir.exists(dirname(output))) {
    stop("Output directory does not exist: ", dirname(output))
  }

  saveRDS(
    seuratObj,
    file = output
  )

  message("RDS created successfully.")

  invisible(output)
}
