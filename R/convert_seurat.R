#' Convert a Seurat Assay5 assay to a classic Assay
#'
#' Converts a selected Seurat v5 `Assay5` assay to a classic Seurat `Assay`
#' while preserving the supported assay data and object structure.
#'
#' @param input Path to an input Seurat `.rds` file.
#' @param output Path to the output Seurat `.rds` file.
#' @param assay Name of the assay to convert. Defaults to `"RNA"`.
#'
#' @return The output path, returned invisibly.
#' @export
convert_seu_v5_to_classic <- function(input, output, assay = "RNA") {
  ## Read object
  obj <- readRDS(input)

  ## Check that the requested assay exists
  available_assays <- SeuratObject::Assays(obj)

  if (!assay %in% available_assays) {
    stop(
      "Assay '", assay, "' was not found in the Seurat object. ",
      "Available assays: ", collapse_or_none(available_assays),
      ". Use inspect_sc() to inspect the available assays."
    )
  }

  ## Access requested assay
  source_assay <- obj[[assay]]

  ## Check that the requested assay uses the Assay5 structure
  if (!inherits(source_assay, "Assay5")) {
    stop(
      "Assay '", assay, "' is not an Assay5. No conversion performed."
    )
  }

  message(
    "\nSource assay '", assay,
    "' is Assay5. Starting conversion to classic Assay..."
  )

  ## Cast Assay5 to classic Assay
  classic_assay <- methods::as(source_assay, Class = "Assay")

  ## Assign unique key
  SeuratObject::Key(classic_assay) <- "classicRNA_"

  ## Add converted assay to object
  obj[["RNA_classic"]] <- classic_assay
  SeuratObject::DefaultAssay(obj) <- "RNA_classic"

  message(
    "\nConversion completed. RNA_classic assay has been added to the Seurat object."
  )

  ## Save RDS
  message("\nSaving object...\n")
  saveRDS(obj, output)
  message("Object has been saved!\n")

  invisible(output)
}

#' Convert a classic Seurat Assay to Assay5
#'
#' Converts a selected classic Seurat `Assay` to a Seurat v5 `Assay5`
#' while preserving the supported assay data and object structure.
#'
#' @param input Path to an input Seurat `.rds` file.
#' @param output Path to the output Seurat `.rds` file.
#' @param assay Name of the assay to convert. Defaults to `"RNA"`.
#'
#' @return The output path, returned invisibly.
#' @export
convert_seu_classic_to_v5 <- function(input, output, assay = "RNA") {
  ## Check installed Seurat version
  seurat_version <- utils::packageVersion("Seurat")

  if (seurat_version < base::package_version("5.0.0")) {
    stop(
      "Converting a classic Assay to Assay5 requires Seurat >= 5.0.0. ",
      "Current Seurat version: ",
      seurat_version
    )
  }

  ## Read object
  obj <- readRDS(input)

  ## Check that the requested assay exists
  available_assays <- SeuratObject::Assays(obj)

  if (!assay %in% available_assays) {
    stop(
      "Assay '", assay, "' was not found in the Seurat object. ",
      "Available assays: ", collapse_or_none(available_assays),
      ". Use inspect_sc() to inspect the available assays."
    )
  }

  ## Access requested assay
  source_assay <- obj[[assay]]

  ## Check that the requested assay uses the classic Assay structure
  if (!inherits(source_assay, "Assay")) {
    stop(
      "Assay '", assay, "' is not a classic Assay. No conversion performed."
    )
  }

  message(
    "\nSource assay '", assay,
    "' is classic. Starting conversion to Assay5..."
  )

  ## Cast classic Assay to Assay5
  assay5 <- methods::as(
    source_assay,
    Class = "Assay5"
  )

  ## Assign unique key
  SeuratObject::Key(assay5) <- "RNAassay5_"

  ## Add converted assay to object
  obj[["RNA_assay5"]] <- assay5
  SeuratObject::DefaultAssay(obj) <- "RNA_assay5"

  message(
    "\nConversion completed. RNA_assay5 has been added to the Seurat object."
  )

  ## Save RDS
  message("\nSaving object...\n")
  saveRDS(obj, output)
  message("Object has been saved!\n")

  invisible(output)
}

#' Convert AnnData to Seurat
#'
#' Reads an AnnData H5AD file and writes a Seurat RDS object while preserving
#' supported counts, normalized expression, metadata, feature metadata, and
#' dimensional reductions.
#'
#' @param input Path to an input `.h5ad` file.
#' @param output Path to the output Seurat `.rds` file.
#'
#' @return The output path, returned invisibly.
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
  ## Reading anndata
  adata <- ad$read_h5ad(input)

  message("Inspecting AnnData structure...")

  n_cells <- reticulate::py_to_r(adata$n_obs)
  n_features <- reticulate::py_to_r(adata$n_vars)

  message("Number of cells: ", n_cells)
  message("Number of features: ", n_features)

  ## Accessing layers
  layers <- reticulate::iterate(
    adata$layers$keys(),
    as.character
  )
  ## Accessing reductions
  reductions <- reticulate::iterate(
    adata$obsm$keys(),
    as.character
  )

  ## Getting cell names
  cell_names <- reticulate::iterate(
    adata$obs_names,
    as.character
  ) |>
    unlist(use.names = FALSE)

  ## Getting feature names
  feature_names <- reticulate::iterate(
    adata$var_names,
    as.character
  ) |>
    unlist(use.names = FALSE)


  # Check layers["counts"] and X for raw counts
  raw_counts <- NULL

  if ("counts" %in% layers) {

    mat <- adata$layers$`__getitem__`("counts")
    ## check if the matrix is sparse
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
    ## check if the matrix has raw counts, and integers
    if (
      all(is.finite(values)) &&
      all(values >= 0)
    ) {

      raw_counts <- reticulate::py_to_r(
        mat$astype("float64")
      )

      message("Using layers['counts'] as Seurat counts.")
    } else {
      stop(
        "layers['counts'] contains negative or non-finite values."
      )
    }
  }

  if (is.null(raw_counts)) {

    mat <- adata$X

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

      message("Using X as Seurat counts.")
    }
  }


  if (is.null(raw_counts)) {

    message(
      "Raw counts-like matrix not found in X or layers['counts'].\n",
      "Layers present: ",
      if (length(layers) > 0) {
        paste(layers, collapse = ", ")
      } else {
        "None"
      },
      "\nConversion stopped."
    )

    stop(
      "Conversion cannot continue because no valid raw counts matrix was found."
    )
  }


  normalized_data <- NULL

  if ("logcounts" %in% layers) {
    normalized_data <- reticulate::py_to_r(
      adata$layers$`__getitem__`("logcounts")$astype("float64")
    )
  } else if ("counts" %in% layers && !is.null(adata$X)) {
    candidate_data <- reticulate::py_to_r(
      adata$X$astype("float64")
    )

    if (
      !identical(dim(candidate_data), dim(raw_counts)) ||
      Matrix::nnzero(candidate_data - raw_counts) > 0
    ) {
      normalized_data <- candidate_data
    }
  }


  # AnnData = cells x features
  # Seurat  = features x cells
  ## Transposing the raw counts matrix
  raw_counts <- Matrix::t(raw_counts)

  if (!is.null(normalized_data)) {
    normalized_data <- Matrix::t(normalized_data)
  }

  ## Metadata
  metadata <- reticulate::py_to_r(adata$obs) |>
    as.data.frame()

  feature_metadata <- reticulate::py_to_r(adata$var) |>
    as.data.frame()

  message(
    "Raw counts: ",
    paste(dim(raw_counts), collapse = " x ")
  )

  message(
    "Metadata: ",
    paste(dim(metadata), collapse = " x ")
  )


  # Restore Seurat dimnames
  rownames(raw_counts) <- feature_names
  colnames(raw_counts) <- cell_names

  if (!is.null(normalized_data)) {
    rownames(normalized_data) <- feature_names
    colnames(normalized_data) <- cell_names
  }

  rownames(metadata) <- cell_names
  rownames(feature_metadata) <- feature_names


  # Alignment checks
  stopifnot(
    identical(colnames(raw_counts), rownames(metadata)),
    identical(rownames(raw_counts), feature_names)
  )

  stopifnot(
    !anyDuplicated(cell_names),
    !anyDuplicated(feature_names)
  )

  message("Matrix alignment checks passed!")


  # Extract AnnData embeddings
  reduction_list <- list()

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
        "' do not match the Seurat cell order."
      )
    }

    reduction_list[[red]] <- emb

    message(
      red, ": ",
      nrow(emb), " cells x ",
      ncol(emb), " dimensions"
    )
  }


  # Construct Seurat object
  message("\nCreating Seurat object...")

  seuratObj <- SeuratObject::CreateSeuratObject(
    counts = raw_counts,
    meta.data = metadata
  )

  if (!is.null(normalized_data)) {
    SeuratObject::LayerData(
      seuratObj,
      assay = "RNA",
      layer = "data"
    ) <- normalized_data
  }

  if (ncol(feature_metadata) > 0) {
    seuratObj[["RNA"]] <- SeuratObject::AddMetaData(
      seuratObj[["RNA"]],
      metadata = feature_metadata
    )
  }

  message(
    "Counts have been added. Now adding reductions..\n"
  )


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

  saveRDS(
    seuratObj,
    file = output
  )

  message(
    "Object created successfully."
  )

  invisible(output)
}

#' Convert Seurat to AnnData
#'
#' Reads a Seurat RDS object and writes an H5AD AnnData file while preserving
#' supported counts, normalized expression, metadata, feature metadata, and
#' dimensional reductions.
#'
#' @param input Path to an input Seurat `.rds` file.
#' @param output Path to the output `.h5ad` file.
#' @param assay Name of the Seurat assay to convert. Defaults to `"RNA"`.
#'
#' @return The output path, returned invisibly.
#' @export
convert_seurat_to_anndata <- function(input, output, assay = "RNA") {
  message("Reading object...")

  ## Reading object
  obj <- readRDS(input)

  message("Inspecting object structure...\n")

  ## Check that the requested assay exists
  available_assays <- SeuratObject::Assays(obj)

  if (!assay %in% available_assays) {
    stop(
      "Assay '", assay, "' was not found in the Seurat object. ",
      "Available assays: ", collapse_or_none(available_assays), "."
    )
  }

  ## Access requested assay
  source_assay <- obj[[assay]]

  ## Extract expression matrices
  if (inherits(source_assay, "Assay5")) {
    message("Assay structure: Assay5\n")

    assay_layers <- SeuratObject::Layers(source_assay)

    counts_mat <- if ("counts" %in% assay_layers) {
      SeuratObject::LayerData(
        source_assay,
        layer = "counts"
      )
    } else if (any(grepl("^counts\\.", assay_layers))) {
      tmp_name <- ".scFlex_counts"
      assay_tmp <- SeuratObject::JoinLayers(
        source_assay,
        layers = "counts",
        new = tmp_name
      )
      SeuratObject::LayerData(
        assay_tmp,
        layer = tmp_name
      )
    } else {
      NULL
    }

    data_mat <- if ("data" %in% assay_layers) {
      SeuratObject::LayerData(
        source_assay,
        layer = "data"
      )
    } else if (any(grepl("^data\\.", assay_layers))) {
      tmp_name <- ".scFlex_data"
      assay_tmp <- SeuratObject::JoinLayers(
        source_assay,
        layers = "data",
        new = tmp_name
      )
      SeuratObject::LayerData(
        assay_tmp,
        layer = tmp_name
      )
    } else {
      NULL
    }

  } else if (inherits(source_assay, "Assay")) {
    message("Assay structure: classic Assay\n")

    counts_mat <- SeuratObject::GetAssayData(
      source_assay,
      layer = "counts"
    )

    data_mat <- SeuratObject::GetAssayData(
      source_assay,
      layer = "data"
    )

    ## Treat an empty data matrix as absent
    if (nrow(data_mat) == 0 || ncol(data_mat) == 0) {
      data_mat <- NULL
    }

  } else {
    stop(
      "Assay '", assay, "' has an unsupported assay structure."
    )
  }

  ## Counts are required for this conversion
  if (is.null(counts_mat)) {
    stop(
      "No counts matrix was found in assay '", assay, "'."
    )
  }

  ## Print matrix dimensions
  message("Counts: ",
          nrow(counts_mat), " features x ",
          ncol(counts_mat), " cells"
  )

  if (!is.null(data_mat)) {
    message(
      "Data: ",
      nrow(data_mat), " features x ",
      ncol(data_mat), " cells"
    )
  } else {
    message("Data: None")
  }

  ## Build AnnData counts matrix (transposed)
  counts <- Matrix::t(counts_mat)

  ## Use normalized data as X when available
  if (!is.null(data_mat)) {
    X <- Matrix::t(data_mat)
    message("Using normalized data as AnnData X.")
  } else {
    X <- counts
    message(
      "Normalized data not found. Using counts as AnnData X."
    )
  }

  ## Cell metadata
  obs <- obj@meta.data
  obs <- sanitize_obs(obs)

  ## Feature metadata
  var <- as.data.frame(source_assay[[]])

  if (ncol(var) == 0) {
    var <- data.frame(
      row.names = rownames(counts_mat)
    )
  } else {
    var <- var[
      rownames(counts_mat),
      ,
      drop = FALSE
    ]
  }

  message(
    "X: ",
    paste(dim(X), collapse = " x ")
  )

  message(
    "Raw counts: ",
    paste(dim(counts), collapse = " x ")
  )

  message(
    "obs: ",
    paste(dim(obs), collapse = " x ")
  )

  message(
    "var: ",
    paste(dim(var), collapse = " x "),
    "\n"
  )

  ## Check matrix alignment
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

  message("Matrix dimension and alignment checks passed!")

  ## Get reductions
  reductions <- SeuratObject::Reductions(obj)

  ## Import AnnData
  reticulate::py_require("anndata>=0.10")
  ad <- reticulate::import("anndata")

  ## Create AnnData object
  adata <- ad$AnnData(
    X = X,
    obs = obs,
    var = var
  )

  ## Add reductions
  if (length(reductions) > 0) {
    for (red in reductions) {
      embeddings <- SeuratObject::Embeddings(
        obj[[red]]
      )

      embeddings <- unname(embeddings)

      key <- paste0("X_", red)

      adata$obsm$`__setitem__`(
        key,
        reticulate::r_to_py(embeddings)
      )
    }
  }

  message("AnnData object created\n")

  ## Add raw counts layer
  adata$layers$`__setitem__`(
    "counts",
    reticulate::r_to_py(counts)
  )

  message("All matrices have been added\n")

  ## Write H5AD file
  message("Writing h5ad output file...\n")

  output <- path.expand(output)

  adata$write_h5ad(
    output,
    convert_strings_to_categoricals = FALSE
  )

  message("H5AD created successfully.\n")

  invisible(output)
}

#' Convert Seurat to SingleCellExperiment
#'
#' Reads a Seurat RDS object and writes a SingleCellExperiment RDS object while
#' preserving supported counts, normalized expression, metadata, feature
#' metadata, and dimensional reductions.
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
  ## Checking assay structure
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
      tmp_name <- ".scTransit_counts"
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
      tmp_name <- ".scTransit_data"
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
      SingleCellExperiment::reducedDim(sce, red) <- emb
      message(red, ": ", nrow(emb), " cells x ",
              ncol(emb), " dimensions")
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
#' Reads a SingleCellExperiment RDS object and writes a Seurat RDS object while
#' preserving supported counts, normalized expression, metadata, feature
#' metadata, and dimensional reductions.
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
    counts = counts_mat,
    meta.data = metadata
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
      key <- if (
        !is.null(colnames(emb)) &&
        length(colnames(emb)) > 0 &&
        grepl("^[A-Za-z][A-Za-z0-9]*_[0-9]+$", colnames(emb)[1])
      ) {
        sub("[0-9]+$", "", colnames(emb)[1])
      } else {
        paste0(
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