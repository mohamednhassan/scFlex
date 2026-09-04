#' Convert Loom to AnnData
#'
#' Reads a Loom file through AnnData and writes the resulting object as H5AD while preserving components supported by the Loom-to-AnnData reader.
#'
#' @details
#' Loom is a more limited interchange format than Seurat, SingleCellExperiment, or AnnData, so only components representable by the current mapping are preserved.
#'
#' @param input Path to an input `.loom` file.
#' @param output Path to the output `.h5ad` file.
#'
#' @return The output path, returned invisibly.
#' @export
convert_loom_to_anndata <- function(input, output) {
  if (!file.exists(input)) {
    stop("File does not exist. Please recheck the path.")
  }
  input <- normalizePath(input, mustWork = TRUE)
  message("Reading Loom object...")
  reticulate::py_require("anndata>=0.10")
  reticulate::py_require("loompy>=3.0")
  ad_io <- reticulate::import("anndata.io", convert = FALSE)
  adata <- ad_io$read_loom(
    input,
    sparse = TRUE
  )
  message("Inspecting Loom structure...")
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
  message(
    "Layers [", length(layers), "]: ",
    if (length(layers) > 0) paste(layers, collapse = ", ") else "None"
  )
  message(
    "Reductions [", length(reductions), "]: ",
    if (length(reductions) > 0) paste(reductions, collapse = ", ") else "None"
  )
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

#' Convert AnnData to Loom
#'
#' Reads an AnnData H5AD file and writes a Loom file using loompy, preserving supported matrices, cell and feature attributes, and reductions where representable.
#'
#' @details
#' Loom is a more limited interchange format than Seurat, SingleCellExperiment, or AnnData, so only components representable by the current mapping are preserved.
#'
#' @param input Path to an input `.h5ad` file.
#' @param output Path to the output `.loom` file.
#'
#' @return The output path, returned invisibly.
#' @export
convert_anndata_to_loom <- function(input, output) {
  if (!file.exists(input)) {
    stop("File does not exist. Please recheck the path.")
  }
  input <- normalizePath(input, mustWork = TRUE)
  message("Reading AnnData object...")
  reticulate::py_require("anndata>=0.10")
  reticulate::py_require("loompy>=3.0")
  ad <- reticulate::import("anndata", convert = FALSE)
  loompy <- reticulate::import("loompy", convert = FALSE)
  scipy_sparse <- reticulate::import("scipy.sparse", convert = FALSE)
  numpy <- reticulate::import("numpy", convert = FALSE)
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
  message(
    "Layers [", length(layers), "]: ",
    if (length(layers) > 0) paste(layers, collapse = ", ") else "None"
  )
  cell_names <- reticulate::iterate(
    adata$obs_names,
    as.character
  ) |>
    unlist(use.names = FALSE)
  feature_names_original <- reticulate::iterate(
    adata$var_names,
    as.character
  ) |>
    unlist(use.names = FALSE)
  stopifnot(
    length(cell_names) == n_cells,
    length(feature_names_original) == n_features
  )
  if (anyDuplicated(cell_names)) {
    stop("Duplicated cell names found in AnnData object.", call. = FALSE)
  }
  feature_names <- feature_names_original
  if (anyDuplicated(feature_names)) {
    n_duplicates <- sum(duplicated(feature_names))
    message(
      "Duplicated feature names detected: ",
      n_duplicates,
      ". Making Loom feature names unique while preserving original names."
    )
    feature_names <- make.unique(feature_names)
  }
  row_attrs <- reticulate::dict()
  col_attrs <- reticulate::dict()
  row_attrs$`__setitem__`(
    "Gene",
    numpy$array(feature_names)
  )
  col_attrs$`__setitem__`(
    "CellID",
    numpy$array(cell_names)
  )
  var <- reticulate::py_to_r(adata$var) |>
    as.data.frame()
  obs <- reticulate::py_to_r(adata$obs) |>
    as.data.frame()
  if (!identical(feature_names, feature_names_original)) {
    var$original_feature_name <- feature_names_original
  }
  if (ncol(var) > 0) {
    for (nm in colnames(var)) {
      values <- var[[nm]]
      if (is.factor(values)) {
        values <- as.character(values)
      }
      if (is.logical(values)) {
        values <- as.integer(values)
      }
      if (is.character(values)) {
        values[is.na(values)] <- ""
      }
      row_attrs$`__setitem__`(
        nm,
        numpy$array(values)
      )
    }
  }
  if (ncol(obs) > 0) {
    for (nm in colnames(obs)) {
      values <- obs[[nm]]
      if (is.factor(values)) {
        values <- as.character(values)
      }
      if (is.logical(values)) {
        values <- as.integer(values)
      }
      if (is.character(values)) {
        values[is.na(values)] <- ""
      }
      col_attrs$`__setitem__`(
        nm,
        numpy$array(values)
      )
    }
  }
  loom_layers <- reticulate::dict()
  X <- adata$X
  if (reticulate::py_to_r(scipy_sparse$issparse(X))) {
    X <- X$T$tocsc()
  } else {
    X <- X$T
  }
  loom_layers$`__setitem__`(
    "",
    X
  )
  if (length(layers) > 0) {
    for (layer_name in layers) {
      mat <- adata$layers$`__getitem__`(layer_name)
      if (reticulate::py_to_r(scipy_sparse$issparse(mat))) {
        mat <- mat$T$tocsc()
      } else {
        mat <- mat$T
      }
      loom_layers$`__setitem__`(
        layer_name,
        mat
      )
    }
  }
  message("Matrix dimensions alignment checks passed!")
  message("Writing Loom output file...")
  output <- path.expand(output)
  if (!dir.exists(dirname(output))) {
    stop("Output directory does not exist: ", dirname(output))
  }
  loompy$create(
    output,
    loom_layers,
    row_attrs = row_attrs,
    col_attrs = col_attrs
  )
  message("Loom created successfully.")
  invisible(output)
}

#' Convert Loom to Seurat
#'
#' Reads a Loom file through AnnData and creates a Seurat object while preserving supported count data, metadata, feature metadata, and reductions.
#'
#' @details
#' Raw counts are selected conservatively from an explicit `counts` layer when available or from `X` when it is non-negative and integer-like. The function does not guess raw counts from arbitrary named layers.
#' Loom is a more limited interchange format than Seurat, SingleCellExperiment, or AnnData, so only components representable by the current mapping are preserved.
#'
#' @param input Path to an input `.loom` file.
#' @param output Path to the output Seurat `.rds` file.
#'
#' @return The output path, returned invisibly on success. If suitable raw count data cannot be identified, the function may return `NULL` invisibly without writing an output object.
#' @export
convert_loom_to_seurat <- function(input, output) {
  if (!file.exists(input)) {
    stop("File does not exist. Please recheck the path.")
  }
  input <- normalizePath(input, mustWork = TRUE)
  message("Reading Loom object...")
  reticulate::py_require("anndata>=0.10")
  reticulate::py_require("loompy>=3.0")
  ad_io <- reticulate::import("anndata.io", convert = FALSE)
  scipy_sparse <- reticulate::import("scipy.sparse", convert = FALSE)
  adata <- ad_io$read_loom(input, sparse = TRUE)
  message("Inspecting Loom structure...")
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
  
  cell_names <- reticulate::iterate(
    adata$obs_names,
    as.character
  ) |>
    unlist(use.names = FALSE)
  
  feature_names_original <- reticulate::iterate(
    adata$var_names,
    as.character
  ) |>
    unlist(use.names = FALSE)
  
  if (anyDuplicated(cell_names)) {
    stop("Duplicated cell names found in Loom file.")
  }
  
  feature_names <- feature_names_original
  
  if (anyDuplicated(feature_names)) {
    n_duplicates <- sum(duplicated(feature_names))
    message(
      "Duplicated feature names detected: ",
      n_duplicates,
      ". Making Seurat feature names unique while preserving original names."
    )
    feature_names <- make.unique(feature_names)
  }
  
  is_count_like <- function(mat) {
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
    
    if (length(values) == 0) {
      return(TRUE)
    }
    
    all(
      is.finite(values) &
        values >= 0 &
        abs(values - round(values)) < 1e-8
    )
  }
  
  counts_source <- NULL
  counts_py <- NULL
  
  if ("counts" %in% layers) {
    candidate <- adata$layers$`__getitem__`("counts")
    
    if (is_count_like(candidate)) {
      counts_py <- candidate
      counts_source <- "counts"
      message("Using layers['counts'] as Seurat counts.")
    }
  }
  
  if (is.null(counts_py) && is_count_like(adata$X)) {
    counts_py <- adata$X
    counts_source <- "X"
    message("Using X as Seurat counts.")
  }
  
  if (is.null(counts_py)) {
    message("Raw counts-like matrix not found in X or layers['counts'].")
    message(
      "Layers present: ",
      if (length(layers) > 0) paste(layers, collapse = ", ") else "None"
    )
    message("Conversion stopped.")
    return(invisible(NULL))
  }
  
  data_py <- NULL
  
  if ("logcounts" %in% layers) {
    data_py <- adata$layers$`__getitem__`("logcounts")
    message("Using layers['logcounts'] as Seurat normalized data.")
  } else if (counts_source != "X") {
    data_py <- adata$X
    message("Using X as Seurat normalized data.")
  } else {
    message("Normalized data not found. Seurat object will contain counts only.")
  }
  
  if (reticulate::py_to_r(scipy_sparse$issparse(counts_py))) {
    counts_mat <- reticulate::py_to_r(
      counts_py$astype("float64")$T$tocsc()
    )
  } else {
    counts_mat <- Matrix::Matrix(
      t(
        reticulate::py_to_r(
          counts_py$astype("float64")
        )
      ),
      sparse = TRUE
    )
  }
  
  rownames(counts_mat) <- feature_names
  colnames(counts_mat) <- cell_names
  
  data_mat <- NULL
  
  if (!is.null(data_py)) {
    if (reticulate::py_to_r(scipy_sparse$issparse(data_py))) {
      data_mat <- reticulate::py_to_r(
        data_py$astype("float64")$T$tocsc()
      )
    } else {
      data_mat <- Matrix::Matrix(
        t(
          reticulate::py_to_r(
            data_py$astype("float64")
          )
        ),
        sparse = TRUE
      )
    }
    
    rownames(data_mat) <- feature_names
    colnames(data_mat) <- cell_names
  }
  
  metadata <- reticulate::py_to_r(adata$obs) |>
    as.data.frame()
  
  rownames(metadata) <- cell_names
  
  feature_metadata <- reticulate::py_to_r(adata$var) |>
    as.data.frame()
  
  if (nrow(feature_metadata) != length(feature_names)) {
    stop("Feature metadata dimensions do not match the expression matrix.")
  }
  
  feature_metadata$original_feature_name <- feature_names_original
  rownames(feature_metadata) <- feature_names
  
  stopifnot(
    identical(colnames(counts_mat), rownames(metadata)),
    identical(rownames(counts_mat), rownames(feature_metadata))
  )
  
  if (!is.null(data_mat)) {
    stopifnot(
      identical(rownames(counts_mat), rownames(data_mat)),
      identical(colnames(counts_mat), colnames(data_mat))
    )
  }
  
  message("Matrix dimensions alignment checks passed!")
  
  seuratObj <- SeuratObject::CreateSeuratObject(
    counts = counts_mat,
    meta.data = metadata,
    assay = "RNA"
  )
  
  if (!is.null(data_mat)) {
    seuratObj <- SeuratObject::SetAssayData(
      seuratObj,
      assay = "RNA",
      layer = "data",
      new.data = data_mat
    )
  }
  
  if (ncol(feature_metadata) > 0) {
    seuratObj[["RNA"]][[]] <- feature_metadata
  }
  
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
  
  if (length(reductions) > 0) {
    for (red in reductions) {
      emb <- reticulate::py_to_r(
        adata$obsm$`__getitem__`(red)
      )
      
      emb <- as.matrix(emb)
      rownames(emb) <- cell_names
      
      red_name <- sub("^X_", "", red)
      
      seuratObj[[red_name]] <- SeuratObject::CreateDimReducObject(
        embeddings = emb,
        key = paste0(toupper(red_name), "_"),
        assay = "RNA"
      )
      
      message(
        red_name, ": ",
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
  
  saveRDS(seuratObj, output)
  
  message("RDS created successfully.")
  invisible(output)
}

#' Convert Seurat to Loom
#'
#' Reads a Seurat RDS object and writes a Loom file using loompy, preserving supported expression matrices, attributes, and reductions where representable.
#'
#' @details
#' Loom is a more limited interchange format than Seurat, SingleCellExperiment, or AnnData, so only components representable by the current mapping are preserved.
#'
#' @param input Path to an input Seurat `.rds` file.
#' @param output Path to the output `.loom` file.
#'
#' @return The output path, returned invisibly.
#' @export
convert_seurat_to_loom <- function(input, output) {
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
    
    if (nrow(data_mat) == 0 || ncol(data_mat) == 0) {
      data_mat <- NULL
    }
  } else if (inherits(assay_obj, "Assay5")) {
    message("Assay structure: Assay5")
    
    layer_names <- SeuratObject::Layers(assay_obj)
    
    count_layers <- grep(
      "^counts($|\\.)",
      layer_names,
      value = TRUE
    )
    
    data_layers <- grep(
      "^data($|\\.)",
      layer_names,
      value = TRUE
    )
    
    if (length(count_layers) == 0) {
      stop("Counts layer not found in Seurat assay.")
    }
    
    if ("counts" %in% layer_names) {
      counts_mat <- SeuratObject::LayerData(
        assay_obj,
        layer = "counts"
      )
    } else {
      tmp_name <- paste0(
        "sctransit_counts_",
        as.integer(Sys.time())
      )
      
      tmp_assay <- SeuratObject::JoinLayers(
        assay_obj,
        layers = "counts",
        new = tmp_name
      )
      
      counts_mat <- SeuratObject::LayerData(
        tmp_assay,
        layer = tmp_name
      )
    }
    
    if (length(data_layers) > 0) {
      if ("data" %in% layer_names) {
        data_mat <- SeuratObject::LayerData(
          assay_obj,
          layer = "data"
        )
      } else {
        tmp_name <- paste0(
          "sctransit_data_",
          as.integer(Sys.time())
        )
        
        tmp_assay <- SeuratObject::JoinLayers(
          assay_obj,
          layers = "data",
          new = tmp_name
        )
        
        data_mat <- SeuratObject::LayerData(
          tmp_assay,
          layer = tmp_name
        )
      }
    }
  } else {
    stop(
      "Unsupported Seurat assay class: ",
      paste(class(assay_obj), collapse = ", ")
    )
  }
  
  message(
    "Counts: ",
    nrow(counts_mat),
    " features x ",
    ncol(counts_mat),
    " cells"
  )
  
  if (!is.null(data_mat)) {
    message(
      "Data: ",
      nrow(data_mat),
      " features x ",
      ncol(data_mat),
      " cells"
    )
  } else {
    message("Data: None")
  }
  
  cell_names <- colnames(counts_mat)
  feature_names <- rownames(counts_mat)
  
  metadata <- obj@meta.data
  metadata <- metadata[cell_names, , drop = FALSE]
  
  feature_metadata <- as.data.frame(
    assay_obj[[]]
  )
  
  if (ncol(feature_metadata) == 0) {
    feature_metadata <- data.frame(
      row.names = feature_names
    )
  } else {
    feature_metadata <- feature_metadata[
      feature_names,
      ,
      drop = FALSE
    ]
  }
  
  stopifnot(
    identical(colnames(counts_mat), rownames(metadata)),
    identical(rownames(counts_mat), rownames(feature_metadata))
  )
  
  if (!is.null(data_mat)) {
    stopifnot(
      identical(rownames(counts_mat), rownames(data_mat)),
      identical(colnames(counts_mat), colnames(data_mat))
    )
  }
  
  stopifnot(
    !anyDuplicated(cell_names),
    !anyDuplicated(feature_names)
  )
  
  message("Matrix dimensions alignment checks passed!")
  
  reticulate::py_require("loompy>=3.0")
  loompy <- reticulate::import("loompy", convert = FALSE)
  scipy_sparse <- reticulate::import("scipy.sparse", convert = FALSE)
  numpy <- reticulate::import("numpy", convert = FALSE)
  
  row_attrs <- reticulate::dict()
  col_attrs <- reticulate::dict()
  
  row_attrs$`__setitem__`(
    "Gene",
    numpy$array(feature_names)
  )
  
  col_attrs$`__setitem__`(
    "CellID",
    numpy$array(cell_names)
  )
  
  if (ncol(feature_metadata) > 0) {
    for (nm in colnames(feature_metadata)) {
      if (nm == "Gene") {
        next
      }
      
      values <- feature_metadata[[nm]]
      
      if (is.factor(values)) {
        values <- as.character(values)
      }
      
      if (is.logical(values)) {
        values <- as.integer(values)
      }
      
      if (is.character(values)) {
        values[is.na(values)] <- ""
      }
      
      if (!is.atomic(values)) {
        values <- as.character(values)
        values[is.na(values)] <- ""
      }
      
      row_attrs$`__setitem__`(
        nm,
        numpy$array(values)
      )
    }
  }
  
  if (ncol(metadata) > 0) {
    for (nm in colnames(metadata)) {
      if (nm == "CellID") {
        next
      }
      
      values <- metadata[[nm]]
      
      if (is.factor(values)) {
        values <- as.character(values)
      }
      
      if (is.logical(values)) {
        values <- as.integer(values)
      }
      
      if (is.character(values)) {
        values[is.na(values)] <- ""
      }
      
      if (!is.atomic(values)) {
        values <- as.character(values)
        values[is.na(values)] <- ""
      }
      
      col_attrs$`__setitem__`(
        nm,
        numpy$array(values)
      )
    }
  }
  
  reductions <- names(obj@reductions)
  
  if (length(reductions) > 0) {
    for (red in reductions) {
      emb <- SeuratObject::Embeddings(
        obj[[red]]
      )
      
      emb <- emb[cell_names, , drop = FALSE]
      
      col_attrs$`__setitem__`(
        paste0("X_", red),
        numpy$array(unname(as.matrix(emb)))
      )
      
      message(
        red, ": ",
        nrow(emb),
        " cells x ",
        ncol(emb),
        " dimensions"
      )
    }
  }
  
  counts_py <- reticulate::r_to_py(counts_mat)
  
  if (reticulate::py_to_r(scipy_sparse$issparse(counts_py))) {
    counts_py <- counts_py$tocsc()
  }
  
  data_py <- NULL
  
  if (!is.null(data_mat)) {
    data_py <- reticulate::r_to_py(data_mat)
    
    if (reticulate::py_to_r(scipy_sparse$issparse(data_py))) {
      data_py <- data_py$tocsc()
    }
  }
  
  loom_layers <- reticulate::dict()
  
  if (!is.null(data_py)) {
    loom_layers$`__setitem__`(
      "",
      data_py
    )
    
    loom_layers$`__setitem__`(
      "counts",
      counts_py
    )
    
    loom_layers$`__setitem__`(
      "logcounts",
      data_py
    )
  } else {
    loom_layers$`__setitem__`(
      "",
      counts_py
    )
    
    loom_layers$`__setitem__`(
      "counts",
      counts_py
    )
  }
  
  message("Writing Loom output file...")
  
  output <- path.expand(output)
  
  if (!dir.exists(dirname(output))) {
    stop("Output directory does not exist: ", dirname(output))
  }
  
  loompy$create(
    output,
    loom_layers,
    row_attrs = row_attrs,
    col_attrs = col_attrs
  )
  
  message("Loom created successfully.")
  invisible(output)
}

#' Convert SingleCellExperiment to Loom
#'
#' Reads a SingleCellExperiment RDS object and writes a Loom file using loompy, preserving supported assays, attributes, and reductions where representable.
#'
#' @details
#' Loom is a more limited interchange format than Seurat, SingleCellExperiment, or AnnData, so only components representable by the current mapping are preserved.
#'
#' @param input Path to an input SingleCellExperiment `.rds` file.
#' @param output Path to the output `.loom` file.
#'
#' @return The output path, returned invisibly.
#' @export
convert_sce_to_loom <- function(input, output) {
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
    nrow(counts_mat),
    " features x ",
    ncol(counts_mat),
    " cells"
  )
  
  if (!is.null(data_mat)) {
    message(
      "Logcounts: ",
      nrow(data_mat),
      " features x ",
      ncol(data_mat),
      " cells"
    )
  } else {
    message("Logcounts: None")
  }
  
  cell_names <- colnames(obj)
  feature_names_original <- rownames(obj)

  if (anyDuplicated(cell_names)) {
    stop(
      "Duplicated cell names found in SingleCellExperiment object.",
      call. = FALSE
    )
  }

  feature_names <- feature_names_original
  if (anyDuplicated(feature_names)) {
    n_duplicates <- sum(duplicated(feature_names))
    message(
      "Duplicated feature names detected: ",
      n_duplicates,
      ". Making Loom feature names unique while preserving original names."
    )
    feature_names <- make.unique(feature_names)
  }

  metadata <- as.data.frame(
    SummarizedExperiment::colData(obj)
  )

  feature_metadata <- as.data.frame(
    SummarizedExperiment::rowData(obj)
  )

  rownames(metadata) <- cell_names

  if (!identical(feature_names, feature_names_original)) {
    feature_metadata$original_feature_name <- feature_names_original
  }

  stopifnot(
    identical(colnames(counts_mat), cell_names),
    identical(rownames(counts_mat), feature_names_original),
    nrow(metadata) == length(cell_names),
    nrow(feature_metadata) == length(feature_names_original)
  )

  if (!is.null(data_mat)) {
    stopifnot(
      identical(rownames(counts_mat), rownames(data_mat)),
      identical(colnames(counts_mat), colnames(data_mat))
    )
  }

  message("Matrix dimensions alignment checks passed!")
  
  reticulate::py_require("loompy>=3.0")
  loompy <- reticulate::import("loompy", convert = FALSE)
  scipy_sparse <- reticulate::import("scipy.sparse", convert = FALSE)
  numpy <- reticulate::import("numpy", convert = FALSE)
  
  row_attrs <- reticulate::dict()
  col_attrs <- reticulate::dict()
  
  row_attrs$`__setitem__`(
    "Gene",
    numpy$array(feature_names)
  )
  
  col_attrs$`__setitem__`(
    "CellID",
    numpy$array(cell_names)
  )
  
  if (ncol(feature_metadata) > 0) {
    for (nm in colnames(feature_metadata)) {
      if (nm == "Gene") {
        next
      }
      
      values <- feature_metadata[[nm]]
      
      if (is.factor(values)) {
        values <- as.character(values)
      }
      
      if (is.logical(values)) {
        values <- as.integer(values)
      }
      
      if (is.character(values)) {
        values[is.na(values)] <- ""
      }
      
      if (!is.atomic(values)) {
        values <- as.character(values)
        values[is.na(values)] <- ""
      }
      
      row_attrs$`__setitem__`(
        nm,
        numpy$array(values)
      )
    }
  }
  
  if (ncol(metadata) > 0) {
    for (nm in colnames(metadata)) {
      if (nm == "CellID") {
        next
      }
      
      values <- metadata[[nm]]
      
      if (is.factor(values)) {
        values <- as.character(values)
      }
      
      if (is.logical(values)) {
        values <- as.integer(values)
      }
      
      if (is.character(values)) {
        values[is.na(values)] <- ""
      }
      
      if (!is.atomic(values)) {
        values <- as.character(values)
        values[is.na(values)] <- ""
      }
      
      col_attrs$`__setitem__`(
        nm,
        numpy$array(values)
      )
    }
  }
  
  if (length(reductions) > 0) {
    for (red in reductions) {
      emb <- SingleCellExperiment::reducedDim(
        obj,
        red
      )
      
      col_attrs$`__setitem__`(
        paste0("X_", red),
        numpy$array(unname(as.matrix(emb)))
      )
      
      message(
        red, ": ",
        nrow(emb),
        " cells x ",
        ncol(emb),
        " dimensions"
      )
    }
  }
  
  counts_py <- reticulate::r_to_py(counts_mat)
  
  if (reticulate::py_to_r(scipy_sparse$issparse(counts_py))) {
    counts_py <- counts_py$tocsc()
  }
  
  data_py <- NULL
  
  if (!is.null(data_mat)) {
    data_py <- reticulate::r_to_py(data_mat)
    
    if (reticulate::py_to_r(scipy_sparse$issparse(data_py))) {
      data_py <- data_py$tocsc()
    }
  }
  
  loom_layers <- reticulate::dict()
  
  if (!is.null(data_py)) {
    loom_layers$`__setitem__`(
      "",
      data_py
    )
    
    loom_layers$`__setitem__`(
      "counts",
      counts_py
    )
    
    loom_layers$`__setitem__`(
      "logcounts",
      data_py
    )
  } else {
    loom_layers$`__setitem__`(
      "",
      counts_py
    )
    
    loom_layers$`__setitem__`(
      "counts",
      counts_py
    )
  }
  
  message("Writing Loom output file...")
  
  output <- path.expand(output)
  
  if (!dir.exists(dirname(output))) {
    stop("Output directory does not exist: ", dirname(output))
  }
  
  loompy$create(
    output,
    loom_layers,
    row_attrs = row_attrs,
    col_attrs = col_attrs
  )
  
  message("Loom created successfully.")
  invisible(output)
}

#' Convert Loom to SingleCellExperiment
#'
#' Reads a Loom file through AnnData and creates a SingleCellExperiment while preserving supported count data, normalized data, metadata, feature metadata, and reductions.
#'
#' @details
#' Raw counts are selected conservatively from an explicit `counts` layer when available or from `X` when it is non-negative and integer-like. The function does not guess raw counts from arbitrary named layers.
#' Loom is a more limited interchange format than Seurat, SingleCellExperiment, or AnnData, so only components representable by the current mapping are preserved.
#'
#' @param input Path to an input `.loom` file.
#' @param output Path to the output SingleCellExperiment `.rds` file.
#'
#' @return The output path, returned invisibly on success. If suitable raw count data cannot be identified, the function may return `NULL` invisibly without writing an output object.
#' @export
convert_loom_to_sce <- function(input, output) {
  if (!file.exists(input)) {
    stop("File does not exist. Please recheck the path.")
  }
  input <- normalizePath(input, mustWork = TRUE)
  message("Reading Loom object...")
  reticulate::py_require("anndata>=0.10")
  reticulate::py_require("loompy>=3.0")
  ad_io <- reticulate::import("anndata.io", convert = FALSE)
  scipy_sparse <- reticulate::import("scipy.sparse", convert = FALSE)
  adata <- ad_io$read_loom(input, sparse = TRUE)
  
  message("Inspecting Loom structure...")
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
  
  cell_names <- reticulate::iterate(
    adata$obs_names,
    as.character
  ) |>
    unlist(use.names = FALSE)
  
  feature_names_original <- reticulate::iterate(
    adata$var_names,
    as.character
  ) |>
    unlist(use.names = FALSE)
  
  if (anyDuplicated(cell_names)) {
    stop("Duplicated cell names found in Loom file.")
  }
  
  feature_names <- feature_names_original
  
  if (anyDuplicated(feature_names)) {
    n_duplicates <- sum(duplicated(feature_names))
    message(
      "Duplicated feature names detected: ",
      n_duplicates,
      ". Making SCE feature names unique while preserving original names."
    )
    feature_names <- make.unique(feature_names)
  }
  
  is_count_like <- function(mat) {
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
    
    if (length(values) == 0) {
      return(TRUE)
    }
    
    all(
      is.finite(values) &
        values >= 0 &
        abs(values - round(values)) < 1e-8
    )
  }
  
  counts_source <- NULL
  counts_py <- NULL
  
  if ("counts" %in% layers) {
    candidate <- adata$layers$`__getitem__`("counts")
    
    if (is_count_like(candidate)) {
      counts_py <- candidate
      counts_source <- "counts"
      message("Using layers['counts'] as SCE counts.")
    }
  }
  
  if (is.null(counts_py) && is_count_like(adata$X)) {
    counts_py <- adata$X
    counts_source <- "X"
    message("Using X as SCE counts.")
  }
  
  if (is.null(counts_py)) {
    message("Raw counts-like matrix not found in X or layers['counts'].")
    message(
      "Layers present: ",
      if (length(layers) > 0) paste(layers, collapse = ", ") else "None"
    )
    message("Conversion stopped.")
    return(invisible(NULL))
  }
  
  data_py <- NULL
  
  if ("logcounts" %in% layers) {
    data_py <- adata$layers$`__getitem__`("logcounts")
    message("Using layers['logcounts'] as SCE logcounts.")
  } else if (counts_source != "X") {
    data_py <- adata$X
    message("Using X as SCE logcounts.")
  } else {
    message("Normalized data not found. SCE object will contain counts only.")
  }
  
  if (reticulate::py_to_r(scipy_sparse$issparse(counts_py))) {
    counts_mat <- reticulate::py_to_r(
      counts_py$astype("float64")$T$tocsc()
    )
  } else {
    counts_mat <- Matrix::Matrix(
      t(
        reticulate::py_to_r(
          counts_py$astype("float64")
        )
      ),
      sparse = TRUE
    )
  }
  
  rownames(counts_mat) <- feature_names
  colnames(counts_mat) <- cell_names
  
  data_mat <- NULL
  
  if (!is.null(data_py)) {
    if (reticulate::py_to_r(scipy_sparse$issparse(data_py))) {
      data_mat <- reticulate::py_to_r(
        data_py$astype("float64")$T$tocsc()
      )
    } else {
      data_mat <- Matrix::Matrix(
        t(
          reticulate::py_to_r(
            data_py$astype("float64")
          )
        ),
        sparse = TRUE
      )
    }
    
    rownames(data_mat) <- feature_names
    colnames(data_mat) <- cell_names
  }
  
  metadata <- suppressWarnings(
    reticulate::py_to_r(adata$obs)
  ) |>
    as.data.frame()
  
  rownames(metadata) <- cell_names
  
  feature_metadata <- suppressWarnings(
    reticulate::py_to_r(adata$var)
  ) |>
    as.data.frame()
  
  if (nrow(feature_metadata) != length(feature_names)) {
    stop("Feature metadata dimensions do not match the expression matrix.")
  }
  
  feature_metadata$original_feature_name <- feature_names_original
  rownames(feature_metadata) <- feature_names
  
  stopifnot(
    identical(colnames(counts_mat), rownames(metadata)),
    identical(rownames(counts_mat), rownames(feature_metadata))
  )
  
  if (!is.null(data_mat)) {
    stopifnot(
      identical(rownames(counts_mat), rownames(data_mat)),
      identical(colnames(counts_mat), colnames(data_mat))
    )
  }
  
  message("Matrix dimensions alignment checks passed!")
  
  sce_assays <- list(
    counts = counts_mat
  )
  
  if (!is.null(data_mat)) {
    sce_assays$logcounts <- data_mat
  }
  
  sce <- SingleCellExperiment::SingleCellExperiment(
    assays = sce_assays,
    colData = S4Vectors::DataFrame(metadata),
    rowData = S4Vectors::DataFrame(feature_metadata)
  )
  
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
  
  if (length(reductions) > 0) {
    for (red in reductions) {
      emb <- reticulate::py_to_r(
        adata$obsm$`__getitem__`(red)
      )
      
      emb <- as.matrix(emb)
      rownames(emb) <- cell_names
      
      red_name <- sub("^X_", "", red)
      
      SingleCellExperiment::reducedDim(
        sce,
        red_name
      ) <- emb
      
      message(
        red_name, ": ",
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
  
  saveRDS(sce, output)
  
  message("RDS created successfully.")
  invisible(output)
}
