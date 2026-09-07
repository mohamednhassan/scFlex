#' Inspect an RDS single-cell object
#'
#' Dispatches an already loaded R object to the Seurat or SingleCellExperiment inspector according to its class.
#'
#' @param obj An R object read from an RDS file.
#'
#' @return An inspection information list returned invisibly by the class-specific inspector.
#' @keywords internal
inspect_rds_object <- function(obj) {
  
  if (inherits(obj, "Seurat")) {
    message("Seurat object detected.\n")
    return(inspect_seurat(obj))
  } else if (inherits(obj, "SingleCellExperiment")) {
    message("SingleCellExperiment object detected.\n")
    return(inspect_sce(obj))
  } else {
    stop("RDS file contains an unsupported object of class: ",
         paste(class(obj), collapse = ", "), call. = FALSE)
  }
}

#' Inspect a single-cell data file
#'
#' Inspects a Seurat, SingleCellExperiment, AnnData, or Loom file and reports its major structural components without converting the file.
#'
#' @details
#' The file extension determines how the file is opened. RDS files are then classified as Seurat or SingleCellExperiment objects.
#' The function prints a concise summary and returns the same information invisibly for programmatic use.
#'
#' @param path_to_file Character string giving the path to an input `.rds`, `.h5ad`, or `.loom` file.
#'
#' @return A named list containing format-specific structural information, returned invisibly.
#' @export
inspect_sc <- function(path_to_file) {
  
  if (!file.exists(path_to_file)) {
    stop("File does not exist. Please recheck the path.", call. = FALSE)
  }
  
  path_to_file <- normalizePath(path_to_file,mustWork = TRUE)
  file_ext <- tools::file_ext(path_to_file) |> tolower()
  if (file_ext == "rds") {
    message("RDS file detected.\n")
    obj <- suppressPackageStartupMessages(readRDS(path_to_file))
    return(inspect_rds_object(obj))
  } else if (file_ext == "h5ad") {
    message("H5AD file detected.\n")
    return(inspect_anndata(path_to_file))
  } else if (file_ext == "loom") {
    message("Loom file detected.\n")
    return(inspect_loom(path_to_file))
  } else {
    stop("Unsupported file extension: .", file_ext, call. = FALSE)
  }
  
}

#' Inspect a Seurat object
#'
#' Reports dimensions, assays, the default assay structure, reductions, graphs, neighbors, metadata columns, and a sample of cell names from a Seurat object.
#'
#' @param obj A Seurat object.
#'
#' @return A named list of structural information, returned invisibly.
#' @keywords internal
inspect_seurat <- function(obj) {
  
  info_list <- list()
  version <- as.character(obj@version)
  
  info_list$num_cells <- ncol(obj)
  info_list$num_features <- nrow(obj)
  info_list$assays <- SeuratObject::Assays(obj)
  info_list$graphs <- SeuratObject::Graphs(obj)
  info_list$neighbors <- SeuratObject::Neighbors(obj)
  info_list$reductions <- SeuratObject::Reductions(obj)
  info_list$metadata_cols <- colnames(obj@meta.data)
  
  info_list$cellnames_sample <- sample(
    colnames(obj),
    size = min(6, ncol(obj)),
    replace = FALSE
  )
  message("Version of Seurat: ", version, "\n")
  message("Number of cells: ", info_list$num_cells)
  message("Number of features: ", info_list$num_features,"\n")
  message("Assays [", 
          length(info_list$assays), "]: ",
          collapse_or_none(info_list$assays))
  
  # RNA assay
  default_assay <- SeuratObject::DefaultAssay(obj)
  info_list$default_assay <- default_assay
  message("Default assay: ", default_assay)
  rna <- obj[[default_assay]]
  if (inherits(rna, "Assay5")) {
    info_list$rna_structure <- "v5"
    info_list$rna_class <- class(rna)[1]
    info_list$rna_layers <- SeuratObject::Layers(rna)
    message("RNA assay structure: Assay5")
    message("RNA layers [", length(info_list$rna_layers), "]: ",
            collapse_or_none(info_list$rna_layers),
            "\n")
  } else if (inherits(rna, "Assay")) {
    info_list$rna_structure <- "classic"
    info_list$rna_class <- class(rna)[1]
    info_list$rna_slots <- methods::slotNames(rna)
    message("RNA assay structure: classic Assay")
    message("RNA slots [", length(info_list$rna_slots), "]: ",
            collapse_or_none(info_list$rna_slots),
            "\n")
  } else {
    info_list$rna_structure <- "other"
    info_list$rna_class <- class(rna)[1]
    message("RNA assay structure: unsupported/other (",
            info_list$rna_class,
            ")\n")
  }
  message("Reductions [", length(info_list$reductions), "]: ",
          collapse_or_none(info_list$reductions))
  message("Graphs [", length(info_list$graphs), "]: ",
          collapse_or_none(info_list$graphs))
  message("Neighbors [", length(info_list$neighbors), "]: ",
          collapse_or_none(info_list$neighbors),
          "\n")
  message("Metadata columns [", 
          length(info_list$metadata_cols), "]: ",
          collapse_or_none(info_list$metadata_cols),
          "\n")
  message("Sample of cell names: ",
          paste(info_list$cellnames_sample, collapse = ", "))
  invisible(info_list)
}

#' Inspect a SingleCellExperiment object
#'
#' Reports dimensions, assays, reduced dimensions, alternative experiments, cell metadata, feature metadata, and a sample of cell names.
#'
#' @param obj A SingleCellExperiment object.
#'
#' @return A named list of structural information, returned invisibly.
#' @keywords internal
inspect_sce <- function(obj) {
  
  info_list <- list()
  # Basic information
  info_list$num_cells <- ncol(obj)
  info_list$num_features <- nrow(obj)
  info_list$assays <- SummarizedExperiment::assayNames(obj)
  info_list$reductions <-
    SingleCellExperiment::reducedDimNames(obj)
  info_list$metadata_cols <-
    colnames(SummarizedExperiment::colData(obj))
  info_list$feature_metadata_cols <-
    colnames(SummarizedExperiment::rowData(obj))
  info_list$alt_experiments <-
    SingleCellExperiment::altExpNames(obj)
  # Sample cell names
  cell_names <- colnames(obj)
  if (!is.null(cell_names) && length(cell_names) > 0) {
    info_list$cellnames_sample <- sample(cell_names,
                                         size = min(6, length(cell_names)),
                                         replace = FALSE)
  } else {
    info_list$cellnames_sample <- character(0)
  }
  # Printing information
  message("Input of class SingleCellExperiment\n")
  message("Number of cells: ", info_list$num_cells)
  message("Number of features: ", info_list$num_features, "\n")
  message("Assays [", 
          length(info_list$assays), "]: ",
          collapse_or_none(info_list$assays))
  message("Reductions [", 
          length(info_list$reductions), "]: ",
          collapse_or_none(info_list$reductions))
  message("Alternative experiments [",
          length(info_list$alt_experiments),
          "]: ",
          collapse_or_none(info_list$alt_experiments),
          "\n")
  message("Metadata columns [",
          length(info_list$metadata_cols),
          "]: ",
          collapse_or_none(info_list$metadata_cols),
          "\n")
  message("Feature metadata columns [",
          length(info_list$feature_metadata_cols),
          "]: ",
          collapse_or_none(info_list$feature_metadata_cols),
          "\n")
  message("Sample of cell names: ",
          collapse_or_none(info_list$cellnames_sample))
  invisible(info_list)
}

#' Inspect an AnnData file
#'
#' Reads an H5AD file with Python AnnData and reports its dimensions, X matrix, layers, reductions, graphs, feature matrices, unstructured metadata, raw slot, metadata columns, and cell names.
#'
#' @param path_to_file Character string giving the path to an `.h5ad` file.
#'
#' @return A named list of AnnData structural information, returned invisibly.
#' @keywords internal
inspect_anndata <- function(path_to_file) {
  message("Reading AnnData object...\n")
  reticulate::py_require("anndata>=0.10")
  ad <- reticulate::import("anndata")
  adata <- ad$read_h5ad(normalizePath(path_to_file, mustWork = TRUE))
  message("Input of class AnnData\n")
  n_cells <- adata$n_obs
  n_features <- adata$n_vars
  message("Number of cells: ", n_cells)
  message("Number of features: ", n_features, "\n")
  x_present <- !is.null(adata$X)
  message("X present: ", if (x_present) "Yes" else "No")
  layers <- adata$layers$keys()$`__iter__`() |>
    reticulate::iterate()
  layers <- unlist(layers, use.names = FALSE)
  layers <- as.character(layers)
  layers <- layers[!is.na(layers) & layers != "" & layers != "None"]
  message(
    "Layers [", length(layers), "]: ",
    if (length(layers) > 0) paste(layers, collapse = ", ") else "None",
    "\n"
  )
  reductions <- adata$obsm$keys()$`__iter__`() |>
    reticulate::iterate()
  reductions <- unlist(reductions, use.names = FALSE)
  reductions <- as.character(reductions)
  reductions <- reductions[!is.na(reductions) & reductions != "" & reductions != "None"]
  message(
    "Reductions (obsm) [", length(reductions), "]: ",
    if (length(reductions) > 0) paste(reductions, collapse = ", ") else "None"
  )
  graphs <- adata$obsp$keys()$`__iter__`() |>
    reticulate::iterate()
  graphs <- unlist(graphs, use.names = FALSE)
  graphs <- as.character(graphs)
  graphs <- graphs[!is.na(graphs) & graphs != "" & graphs != "None"]
  message(
    "Graphs (obsp) [", length(graphs), "]: ",
    if (length(graphs) > 0) paste(graphs, collapse = ", ") else "None"
  )
  feature_matrices <- adata$varm$keys()$`__iter__`() |>
    reticulate::iterate()
  feature_matrices <- unlist(feature_matrices, use.names = FALSE)
  feature_matrices <- as.character(feature_matrices)
  feature_matrices <- feature_matrices[
    !is.na(feature_matrices) &
      feature_matrices != "" &
      feature_matrices != "None"
  ]
  message(
    "Feature matrices (varm) [", length(feature_matrices), "]: ",
    if (length(feature_matrices) > 0) paste(feature_matrices, collapse = ", ") else "None"
  )
  uns <- names(adata$uns)
  if (is.null(uns)) {
    uns <- character(0)
  }
  uns <- as.character(uns)
  uns <- uns[!is.na(uns) & uns != "" & uns != "None"]
  message(
    "Unstructured metadata (uns) [", length(uns), "]: ",
    if (length(uns) > 0) paste(uns, collapse = ", ") else "None"
  )
  raw_present <- !is.null(adata$raw)
  message("Raw present: ", if (raw_present) "Yes" else "No", "\n")
  metadata_columns <- colnames(as.data.frame(adata$obs))
  message(
    "Metadata columns (obs) [", length(metadata_columns), "]: ",
    if (length(metadata_columns) > 0) paste(metadata_columns, collapse = ", ") else "None",
    "\n"
  )
  feature_metadata_columns <- colnames(as.data.frame(adata$var))
  message(
    "Feature metadata columns (var) [", length(feature_metadata_columns), "]: ",
    if (length(feature_metadata_columns) > 0) paste(feature_metadata_columns, collapse = ", ") else "None",
    "\n"
  )
  cell_names <- adata$obs_names$to_list()
  cell_names <- unlist(cell_names, use.names = FALSE)
  if (length(cell_names) > 0) {
    cellnames_sample <- sample(
      cell_names,
      size = min(6, length(cell_names)),
      replace = FALSE
    )
  } else {
    cellnames_sample <- character(0)
  }
  message(
    "Sample of cell names: ",
    if (length(cellnames_sample) > 0) paste(cellnames_sample, collapse = ", ") else "None"
  )
  info_list <- list(
    class = "AnnData",
    n_cells = n_cells,
    n_features = n_features,
    x_present = x_present,
    layers = layers,
    reductions = reductions,
    graphs = graphs,
    feature_matrices = feature_matrices,
    uns = uns,
    raw_present = raw_present,
    metadata_columns = metadata_columns,
    feature_metadata_columns = feature_metadata_columns,
    cellnames_sample = cellnames_sample
  )
  invisible(info_list)
}

#' Inspect a Loom file
#'
#' Inspects the HDF5 structure of a Loom file and reports matrix dimensions, layers, row and column attributes, graphs, and a sample of cell names when available.
#'
#' @param obj Character string giving the path to a `.loom` file.
#'
#' @return A named list of Loom structural information, returned invisibly.
#' @keywords internal
inspect_loom <- function(obj) {
  if (!requireNamespace("hdf5r", quietly = TRUE)) {
    stop(
      "The 'hdf5r' package is required to inspect Loom files.",
      call. = FALSE
    )
  }
  message("Reading Loom object...\n")
  loom <- hdf5r::H5File$new(obj, mode = "r")
  on.exit(loom$close_all())
  info_list <- list()
  root_names <- names(loom)
  if ("matrix" %in% root_names) {
    matrix_dims <- loom[["matrix"]]$dims
    info_list$num_cells <- matrix_dims[1]
    info_list$num_features <- matrix_dims[2]
  } else {
    info_list$num_features <- NA_integer_
    info_list$num_cells <- NA_integer_
  }
  
  if ("layers" %in% root_names) {
    info_list$layers <- names(loom[["layers"]])
  } else {
    info_list$layers <- character(0)
  }
  
  if ("col_attrs" %in% root_names) {
    info_list$metadata_cols <- names(loom[["col_attrs"]])
  } else {
    info_list$metadata_cols <- character(0)
  }
  
  if ("row_attrs" %in% root_names) {
    info_list$feature_metadata_cols <- names(loom[["row_attrs"]])
  } else {
    info_list$feature_metadata_cols <- character(0)
  }
  
  if ("col_graphs" %in% root_names) {
    info_list$col_graphs <- names(loom[["col_graphs"]])
  } else {
    info_list$col_graphs <- character(0)
  }
  
  if ("row_graphs" %in% root_names) {
    info_list$row_graphs <- names(loom[["row_graphs"]])
  } else {
    info_list$row_graphs <- character(0)
  }
  
  cell_names <- character(0)
  
  if ("col_attrs" %in% root_names) {
    col_attrs <- loom[["col_attrs"]]
    
    possible_cell_names <- c(
      "CellID",
      "cell_names",
      "CellName",
      "Barcode",
      "barcodes"
    )
    
    cell_name_field <- possible_cell_names[
      possible_cell_names %in% names(col_attrs)
    ]
    
    if (length(cell_name_field) > 0) {
      cell_names <- as.character(
        col_attrs[[cell_name_field[1]]][]
      )
    }
  }
  
  if (length(cell_names) > 0) {
    info_list$cellnames_sample <- sample(
      cell_names,
      size = min(6, length(cell_names)),
      replace = FALSE
    )
  } else {
    info_list$cellnames_sample <- character(0)
  }
  
  message("Input of class Loom\n")
  message("Number of cells: ",
          info_list$num_cells)
  message("Number of features: ",
          info_list$num_features,
          "\n")
  message("Layers [",
          length(info_list$layers),
          "]: ",
          collapse_or_none(info_list$layers),
          "\n")
  message("Metadata columns / col_attrs [",
          length(info_list$metadata_cols),
          "]: ",
          collapse_or_none(info_list$metadata_cols),
          "\n")
  message("Feature metadata columns / row_attrs [",
          length(info_list$feature_metadata_cols),
          "]: ",
          collapse_or_none(info_list$feature_metadata_cols),
          "\n")
  message("Column graphs / col_graphs [",
          length(info_list$col_graphs),
          "]: ",
          collapse_or_none(info_list$col_graphs))
  message("Row graphs / row_graphs [",
          length(info_list$row_graphs),
          "]: ",
          collapse_or_none(info_list$row_graphs),
          "\n")
  message("Sample of cell names: ",
          collapse_or_none(info_list$cellnames_sample))
  
  invisible(info_list)
}
