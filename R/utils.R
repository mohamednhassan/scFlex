.collapse_or_none <- function(x) {
  if (length(x) == 0L) return("None")
  paste(x, collapse = ", ")
}

.validate_input_file <- function(path) {
  if (!is.character(path) || length(path) != 1L || is.na(path) || !nzchar(path)) {
    stop("`input` must be a single non-empty file path.", call. = FALSE)
  }
  if (!file.exists(path)) {
    stop("Input file does not exist: ", path, call. = FALSE)
  }
  invisible(TRUE)
}

.validate_output_file <- function(path, overwrite = FALSE) {
  if (!is.character(path) || length(path) != 1L || is.na(path) || !nzchar(path)) {
    stop("`output` must be a single non-empty file path.", call. = FALSE)
  }
  out_dir <- dirname(path)
  if (!dir.exists(out_dir)) {
    stop("Output directory does not exist: ", out_dir, call. = FALSE)
  }
  if (file.exists(path) && !isTRUE(overwrite)) {
    stop(
      "Output file already exists: ", path,
      "\nSet `overwrite = TRUE` to replace it.",
      call. = FALSE
    )
  }
  invisible(TRUE)
}

.require_namespace <- function(pkg, reason = NULL) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    msg <- paste0("Package '", pkg, "' is required")
    if (!is.null(reason)) msg <- paste0(msg, " ", reason)
    stop(msg, ".", call. = FALSE)
  }
  invisible(TRUE)
}

.import_anndata <- function(require_loom = FALSE) {
  .require_namespace("reticulate", "for AnnData/Loom conversion")
  reticulate::py_require("anndata>=0.10")
  if (isTRUE(require_loom)) {
    reticulate::py_require("loompy>=3.0")
  }
  tryCatch(
    reticulate::import("anndata", convert = TRUE),
    error = function(e) {
      stop(
        "Python package 'anndata' could not be imported.\n",
        conditionMessage(e),
        call. = FALSE
      )
    }
  )
}

.python_keys <- function(x) {
  if (is.null(x)) return(character())
  out <- tryCatch(
    x$keys()$`__iter__`() |> reticulate::iterate(),
    error = function(e) character()
  )
  as.character(out)
}

.sanitize_dataframe <- function(x, axis = "metadata") {
  x <- as.data.frame(x)
  ids <- rownames(x)

  x[] <- lapply(names(x), function(column) {
    value <- x[[column]]

    if (is.factor(value)) return(value)
    if (is.atomic(value) && is.null(dim(value))) return(value)
    if (is.matrix(value) && ncol(value) == 1L) return(as.vector(value[, 1]))
    if (is.data.frame(value) && ncol(value) == 1L) return(value[[1]])

    stop(
      axis, " column '", column,
      "' is not one-dimensional and cannot be represented safely. ",
      "Class: ", paste(class(value), collapse = ", "),
      if (!is.null(dim(value))) paste0("; dimensions: ", paste(dim(value), collapse = " x ")) else "",
      call. = FALSE
    )
  })

  rownames(x) <- ids
  x
}

.make_reduction_key <- function(name) {
  clean <- toupper(gsub("[^A-Za-z0-9]", "", name))
  if (!nzchar(clean)) clean <- "REDUCTION"
  paste0(clean, "_")
}

.detect_format <- function(input, format = NULL) {
  valid <- c("seurat", "sce", "anndata", "loom")

  if (!is.null(format)) {
    format <- tolower(format)
    if (!format %in% valid) {
      stop("Unsupported format '", format, "'. Choose from: ", paste(valid, collapse = ", "), ".", call. = FALSE)
    }
    return(format)
  }

  ext <- tolower(tools::file_ext(input))
  if (ext == "h5ad") return("anndata")
  if (ext == "loom") return("loom")

  if (ext %in% c("rds", "rda", "rdata")) {
    if (ext != "rds") {
      stop("Automatic R-object detection currently requires an `.rds` file. Supply `from` explicitly for other R formats.", call. = FALSE)
    }
    obj <- readRDS(input)
    if (inherits(obj, "Seurat")) return("seurat")
    if (inherits(obj, "SingleCellExperiment")) return("sce")
  }

  stop(
    "Could not determine the input format from '", input,
    "'. Supply `from =` explicitly.",
    call. = FALSE
  )
}

.target_from_extension <- function(output, to = NULL) {
  valid <- c("seurat", "sce", "anndata", "loom")
  if (!is.null(to)) {
    to <- tolower(to)
    if (!to %in% valid) stop("Unsupported target format '", to, "'.", call. = FALSE)
    return(to)
  }

  ext <- tolower(tools::file_ext(output))
  if (ext == "h5ad") return("anndata")
  if (ext == "loom") return("loom")
  if (ext == "rds") {
    stop("`.rds` can contain either Seurat or SingleCellExperiment. Please supply `to = 'seurat'` or `to = 'sce'`.", call. = FALSE)
  }
  stop("Could not determine target format from output extension. Supply `to =` explicitly.", call. = FALSE)
}


.assay5_layer_matches <- function(layers, layer) {
  layers[layers == layer | startsWith(layers, paste0(layer, "."))]
}

.read_assay5_layer <- function(assay_obj, layer, cell_order = NULL, feature_order = NULL) {
  layers <- SeuratObject::Layers(assay_obj)
  matches <- .assay5_layer_matches(layers, layer)

  if (length(matches) == 0L) {
    return(list(matrix = NULL, layers = character(), joined = FALSE))
  }

  # A canonical layer takes precedence when it already exists.
  if (layer %in% matches) {
    mat <- SeuratObject::LayerData(assay_obj, layer = layer)
    used <- layer
    joined <- FALSE
  } else if (length(matches) == 1L) {
    mat <- SeuratObject::LayerData(assay_obj, layer = matches[[1L]])
    used <- matches
    joined <- FALSE
  } else {
    tmp_name <- paste0("sctransit_", gsub("[^A-Za-z0-9_.]", "_", layer))
    while (tmp_name %in% SeuratObject::Layers(assay_obj)) {
      tmp_name <- paste0(tmp_name, "_")
    }

    joined_assay <- SeuratObject::JoinLayers(
      assay_obj,
      layers = matches,
      new = tmp_name
    )
    mat <- SeuratObject::LayerData(joined_assay, layer = tmp_name)
    used <- matches
    joined <- TRUE
  }

  if (!is.null(feature_order) && setequal(rownames(mat), feature_order)) {
    mat <- mat[feature_order, , drop = FALSE]
  }
  if (!is.null(cell_order) && setequal(colnames(mat), cell_order)) {
    mat <- mat[, cell_order, drop = FALSE]
  }

  list(matrix = mat, layers = used, joined = joined)
}

.new_bundle <- function(
  counts = NULL,
  data = NULL,
  cell_metadata = NULL,
  feature_metadata = NULL,
  reductions = list(),
  source_format = NULL,
  source_assay = NULL,
  notes = character()
) {
  structure(
    list(
      counts = counts,
      data = data,
      cell_metadata = cell_metadata,
      feature_metadata = feature_metadata,
      reductions = reductions,
      source_format = source_format,
      source_assay = source_assay,
      notes = notes
    ),
    class = "scTransit_bundle"
  )
}

.validate_bundle <- function(x) {
  if (!inherits(x, "scTransit_bundle")) {
    stop("Internal conversion bundle is invalid.", call. = FALSE)
  }

  mats <- Filter(Negate(is.null), list(counts = x$counts, data = x$data))
  if (length(mats) == 0L) {
    stop("No expression matrix was found in the source object.", call. = FALSE)
  }

  reference <- mats[[1L]]
  if (is.null(rownames(reference)) || is.null(colnames(reference))) {
    stop("Expression matrices require feature and cell names.", call. = FALSE)
  }
  if (anyDuplicated(rownames(reference))) stop("Duplicated feature names detected.", call. = FALSE)
  if (anyDuplicated(colnames(reference))) stop("Duplicated cell names detected.", call. = FALSE)

  for (nm in names(mats)) {
    mat <- mats[[nm]]
    if (!identical(dim(mat), dim(reference)) ||
        !identical(rownames(mat), rownames(reference)) ||
        !identical(colnames(mat), colnames(reference))) {
      stop("Expression matrix '", nm, "' is not aligned with the other matrices.", call. = FALSE)
    }
  }

  if (!is.null(x$cell_metadata) && !identical(rownames(x$cell_metadata), colnames(reference))) {
    stop("Cell metadata are not aligned with expression-matrix columns.", call. = FALSE)
  }
  if (!is.null(x$feature_metadata) && !identical(rownames(x$feature_metadata), rownames(reference))) {
    stop("Feature metadata are not aligned with expression-matrix rows.", call. = FALSE)
  }

  for (nm in names(x$reductions)) {
    emb <- x$reductions[[nm]]
    if (nrow(emb) != ncol(reference)) {
      stop("Reduction '", nm, "' has ", nrow(emb), " rows but the object has ", ncol(reference), " cells.", call. = FALSE)
    }
  }

  invisible(TRUE)
}

.conversion_result <- function(output, from, to, notes = character()) {
  structure(
    list(output = output, from = from, to = to, notes = notes),
    class = "scTransit_result"
  )
}
