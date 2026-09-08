#' Inspect an R single-cell object
#'
#' Internal helper that identifies whether an object read from an RDS file is
#' a Seurat or SingleCellExperiment object and dispatches it to the
#' corresponding inspection function.
#'
#' @param obj An R object read from an RDS file.
#'
#' @return Invisibly returns a list containing information about the object.
#' @keywords internal
inspect_rds_object <- function(obj) {
  ## If seurat
  if (inherits(obj, "Seurat")) {
    message("Seurat object detected.\n")
    return(inspect_seurat(obj))
    ## If SCE
  } else if (inherits(obj, "SingleCellExperiment")) {
    message("SingleCellExperiment object detected.\n")
    return(inspect_sce(obj))
  } else {
    stop("RDS file contains an unsupported object of class: ",
         paste(class(obj), collapse = ", "), call. = FALSE)
  }
}
#' Inspect a single-cell object
#'
#' Inspects the structure and contents of a supported single-cell data file
#' without performing a conversion. Supported inputs include Seurat and
#' SingleCellExperiment objects stored as RDS files, AnnData H5AD files,
#' and Loom files.
#'
#' The function reports basic information such as the number of cells and
#' features, available assays or layers, dimensional reductions, metadata,
#' and a sample of cell names. Additional format-specific information is
#' reported when available.
#'
#' @param path_to_file Path to an RDS, H5AD, or Loom file.
#'
#' @return Invisibly returns a list containing information about the
#'   inspected object. The contents of the list depend on the input format.
#'
#' @examples
#' \dontrun{
#' inspect_sc("example.rds")
#' inspect_sc("example.h5ad")
#' inspect_sc("example.loom")
#' }
#'
#' @export
inspect_sc <- function(path_to_file) {
  ## check if file exists
  if (!file.exists(path_to_file)) {
    stop("File does not exist. Please recheck the path.", call. = FALSE)
  }
  ## normalize path to remove any special characters (~, etc.) and keep it in lower case
  path_to_file <- normalizePath(path_to_file,mustWork = TRUE)
  file_ext <- tools::file_ext(path_to_file) |> tolower()
  ## Check if the file has rds extenstion
  if (file_ext == "rds") {
    message("RDS file detected.\n")
    obj <- suppressPackageStartupMessages(readRDS(path_to_file))
    return(inspect_rds_object(obj))
    ## Check if the file has h5ad extenstion
  } else if (file_ext == "h5ad") {
    message("H5AD file detected.\n")
    return(inspect_anndata(path_to_file))
    ## Check if the file has loom extenstion
  } else if (file_ext == "loom") {
    message("Loom file detected.\n")
    return(inspect_loom(path_to_file))
  } else {
    stop("Unsupported file extension: .", file_ext, call. = FALSE)
  }

}
#' Inspect a Seurat object
#'
#' Internal helper that reports the structure of a Seurat object, including
#' its dimensions, assays, default assay structure, reductions, graphs,
#' neighbors, metadata columns, and a sample of cell names.
#'
#' @param obj A Seurat object.
#'
#' @return Invisibly returns a list containing information about the
#'   Seurat object.
#' @keywords internal
inspect_seurat <- function(obj) {
  ## Initialize a list of information
  info_list <- list()

  version <- as.character(obj@version)
  ## Basic object information
  info_list$num_cells <- ncol(obj)
  info_list$num_features <- nrow(obj)
  info_list$assays <- SeuratObject::Assays(obj)
  info_list$graphs <- SeuratObject::Graphs(obj)
  info_list$neighbors <- SeuratObject::Neighbors(obj)
  info_list$reductions <- SeuratObject::Reductions(obj)
  info_list$metadata_cols <- colnames(obj@meta.data)

  ## Random sample of cell names
  info_list$cellnames_sample <- sample(
    colnames(obj),
    size = min(6, ncol(obj)),
    replace = FALSE
  )

  message("Version of Seurat: ", version, "\n")
  message("Number of cells: ", info_list$num_cells)
  message("Number of features: ", info_list$num_features, "\n")

  ## Available assays
  message(
    "Assays [", length(info_list$assays), "]: ",
    collapse_or_none(info_list$assays)
  )

  ## Default assay
  default_assay <- SeuratObject::DefaultAssay(obj)
  info_list$default_assay <- default_assay
  message("Default assay: ", default_assay)
  assay <- obj[[default_assay]]

  ## Check if the default assay uses the Seurat v5 Assay5 structure
  if (inherits(assay, "Assay5")) {
    info_list$assay_structure <- "v5"
    info_list$assay_class <- class(assay)[1]
    info_list$assay_layers <- SeuratObject::Layers(assay)

    message("Default assay structure: Assay5")
    message(
      "Default assay layers [", length(info_list$assay_layers), "]: ",
      collapse_or_none(info_list$assay_layers),
      "\n"
    )

    ## Check if the default assay uses the classic Seurat Assay structure
  } else if (inherits(assay, "Assay")) {
    info_list$assay_structure <- "classic"
    info_list$assay_class <- class(assay)[1]
    info_list$assay_slots <- methods::slotNames(assay)

    message("Default assay structure: classic Assay")
    message(
      "Default assay slots [", length(info_list$assay_slots), "]: ",
      collapse_or_none(info_list$assay_slots),
      "\n"
    )

    ## Catch unsupported or unexpected assay classes
  } else {
    info_list$assay_structure <- "other"
    info_list$assay_class <- class(assay)[1]

    message(
      "Default assay structure: unsupported/other (",
      info_list$assay_class,
      ")\n"
    )
  }

  ## Reductions, graphs and neighbors
  message(
    "Reductions [", length(info_list$reductions), "]: ",
    collapse_or_none(info_list$reductions)
  )

  message(
    "Graphs [", length(info_list$graphs), "]: ",
    collapse_or_none(info_list$graphs)
  )

  message(
    "Neighbors [", length(info_list$neighbors), "]: ",
    collapse_or_none(info_list$neighbors),
    "\n"
  )

  ## Cell metadata
  message(
    "Metadata columns [", length(info_list$metadata_cols), "]: ",
    collapse_or_none(info_list$metadata_cols),
    "\n"
  )

  ## Example cell names
  message(
    "Sample of cell names: ",
    paste(info_list$cellnames_sample, collapse = ", ")
  )

  invisible(info_list)
}
#' Inspect a SingleCellExperiment object
#'
#' Internal helper that reports the structure of a SingleCellExperiment
#' object, including its dimensions, assays, reduced dimensions, alternative
#' experiments, cell metadata, feature metadata, and a sample of cell names.
#'
#' @param obj A SingleCellExperiment object.
#'
#' @return Invisibly returns a list containing information about the
#'   SingleCellExperiment object.
#' @keywords internal
inspect_sce <- function(obj) {
  ## Initialize information list
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
  message("Alteassaytive experiments [",
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
#' Inspect an AnnData object
#'
#' Internal helper that reads an H5AD file and reports its dimensions,
#' primary matrix, layers, cell-level matrices, graphs, feature-level
#' matrices, unstructured metadata, raw data, metadata columns, and a
#' sample of cell names.
#'
#' @param path_to_file Path to an H5AD file.
#'
#' @return Invisibly returns a list containing information about the
#'   AnnData object.
#' @keywords internal
inspect_anndata <- function(path_to_file) {
  message("Reading AnnData object...\n")

  reticulate::py_require("anndata>=0.10")
  ad <- reticulate::import("anndata")

  ## Read the h5ad file
  adata <- ad$read_h5ad(normalizePath(path_to_file, mustWork = TRUE))
  message("Input of class AnnData\n")

  ## Number of cells and genes present
  n_cells <- adata$n_obs
  n_features <- adata$n_vars

  message("Number of cells: ", n_cells)
  message("Number of features: ", n_features, "\n")

  ## Check if .X matrix is present
  x_present <- !is.null(adata$X)
  message("X present: ", if (x_present) "Yes" else "No")

  ## Iterate over layers
  layers <- adata$layers$keys()$`__iter__`() |>
    reticulate::iterate()

  layers <- unlist(layers, use.names = FALSE)
  layers <- as.character(layers)
  layers <- layers[
    !is.na(layers) &
      layers != "" &
      layers != "None"
  ]

  message(
    "Layers [", length(layers), "]: ",
    collapse_or_none(layers),
    "\n"
  )

  ## Iterate over reductions
  reductions <- adata$obsm$keys()$`__iter__`() |>
    reticulate::iterate()

  reductions <- unlist(reductions, use.names = FALSE)
  reductions <- as.character(reductions)
  reductions <- reductions[
    !is.na(reductions) &
      reductions != "" &
      reductions != "None"
  ]

  message(
    "Reductions (obsm) [", length(reductions), "]: ",
    collapse_or_none(reductions)
  )

  ## Iterate over graphs
  graphs <- adata$obsp$keys()$`__iter__`() |>
    reticulate::iterate()

  graphs <- unlist(graphs, use.names = FALSE)
  graphs <- as.character(graphs)
  graphs <- graphs[
    !is.na(graphs) &
      graphs != "" &
      graphs != "None"
  ]

  message(
    "Graphs (obsp) [", length(graphs), "]: ",
    collapse_or_none(graphs)
  )

  ## Check gene loadings
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
    collapse_or_none(feature_matrices)
  )

  ## Check unstructured metadata
  uns <- names(adata$uns)

  if (is.null(uns)) {
    uns <- character(0)
  }

  uns <- as.character(uns)
  uns <- uns[
    !is.na(uns) &
      uns != "" &
      uns != "None"
  ]

  message(
    "Unstructured metadata (uns) [", length(uns), "]: ",
    collapse_or_none(uns)
  )

  ## Check if raw matrix is present
  raw_present <- !is.null(adata$raw)
  message("Raw present: ", if (raw_present) "Yes" else "No", "\n")

  ## Cell metadata columns
  metadata_columns <- colnames(as.data.frame(adata$obs))

  message(
    "Metadata columns (obs) [", length(metadata_columns), "]: ",
    collapse_or_none(metadata_columns),
    "\n"
  )

  ## Feature metadata columns
  feature_metadata_columns <- colnames(as.data.frame(adata$var))

  message(
    "Feature metadata columns (var) [", length(feature_metadata_columns), "]: ",
    collapse_or_none(feature_metadata_columns),
    "\n"
  )

  ## Capture cell names
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
    collapse_or_none(cellnames_sample)
  )

  ## Update the information list
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
#' Internal helper that reads the HDF5 structure of a Loom file and reports
#' its dimensions, layers, cell and feature attributes, graphs, and a sample
#' of cell names when available.
#'
#' @param obj Path to a Loom file.
#'
#' @return Invisibly returns a list containing information about the
#'   Loom file.
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
  ## Initialize information list
  info_list <- list()
  ## Get top-level loom entries
  root_names <- names(loom)
  if ("matrix" %in% root_names) {
    ## Get matrix dimensions
    matrix_dims <- loom[["matrix"]]$dims
    ## Number of cells
    info_list$num_cells <- matrix_dims[1]
    ## Number of features
    info_list$num_features <- matrix_dims[2]
  } else {
    info_list$num_features <- NA_integer_
    info_list$num_cells <- NA_integer_
  }
  ## Access layers
  if ("layers" %in% root_names) {
    info_list$layers <- names(loom[["layers"]])
  } else {
    info_list$layers <- character(0)
  }
  ## Access cell metadata
  if ("col_attrs" %in% root_names) {
    info_list$metadata_cols <- names(loom[["col_attrs"]])
  } else {
    info_list$metadata_cols <- character(0)
  }
  ## Access feature metadata
  if ("row_attrs" %in% root_names) {
    info_list$feature_metadata_cols <- names(loom[["row_attrs"]])
  } else {
    info_list$feature_metadata_cols <- character(0)
  }
  ## Access cell graphs
  if ("col_graphs" %in% root_names) {
    info_list$col_graphs <- names(loom[["col_graphs"]])
  } else {
    info_list$col_graphs <- character(0)
  }
  ## Access feature graphs
  if ("row_graphs" %in% root_names) {
    info_list$row_graphs <- names(loom[["row_graphs"]])
  } else {
    info_list$row_graphs <- character(0)
  }

  ## Capture cell names if available
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
  message("Number of cells: ", info_list$num_cells)
  message("Number of features: ", info_list$num_features, "\n")
  message("Layers [", length(info_list$layers), "]: ",
          collapse_or_none(info_list$layers),"\n")
  message("Metadata columns / col_attrs [", length(info_list$metadata_cols),"]: ",
          collapse_or_none(info_list$metadata_cols),"\n")
  message("Feature metadata columns / row_attrs [",
          length(info_list$feature_metadata_cols),
          "]: ", collapse_or_none(info_list$feature_metadata_cols), "\n")
  message("Column graphs / col_graphs [", length(info_list$col_graphs), "]: ",
          collapse_or_none(info_list$col_graphs))
  message("Row graphs / row_graphs [", length(info_list$row_graphs), "]: ",
          collapse_or_none(info_list$row_graphs), "\n")
  message("Sample of cell names: ", collapse_or_none(info_list$cellnames_sample))
  invisible(info_list)
}