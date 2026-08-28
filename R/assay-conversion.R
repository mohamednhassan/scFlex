#' Convert a Seurat Assay5 to a classic Assay
#'
#' @param input Path to a Seurat `.rds` file.
#' @param output Output `.rds` path.
#' @param assay Name of the Assay5 to convert.
#' @param new_assay Name for the converted classic assay.
#' @param overwrite Replace an existing output file.
#' @return Invisibly returns the output path.
#' @export
convert_seu_v5_to_classic <- function(input, output, assay = "RNA", new_assay = "RNA_classic", overwrite = FALSE) {
  .validate_input_file(input)
  .validate_output_file(output, overwrite)
  obj <- readRDS(input)
  if (!inherits(obj, "Seurat")) stop("Input is not a Seurat object.", call. = FALSE)
  if (!assay %in% SeuratObject::Assays(obj)) stop("Assay '", assay, "' was not found.", call. = FALSE)
  if (!inherits(obj[[assay]], "Assay5")) stop("Assay '", assay, "' is not an Assay5.", call. = FALSE)
  if (new_assay %in% SeuratObject::Assays(obj)) stop("Assay '", new_assay, "' already exists.", call. = FALSE)

  converted <- methods::as(obj[[assay]], "Assay")
  SeuratObject::Key(converted) <- paste0(gsub("[^A-Za-z0-9]", "", new_assay), "_")
  obj[[new_assay]] <- converted
  saveRDS(obj, output)
  message("Converted assay added as '", new_assay, "'. Saved: ", output)
  invisible(output)
}

#' Convert a classic Seurat Assay to Assay5
#'
#' @param input Path to a Seurat `.rds` file.
#' @param output Output `.rds` path.
#' @param assay Name of the classic assay to convert.
#' @param new_assay Name for the resulting Assay5.
#' @param overwrite Replace an existing output file.
#' @return Invisibly returns the output path.
#' @export
convert_seu_classic_to_v5 <- function(input, output, assay = "RNA", new_assay = "RNA_assay5", overwrite = FALSE) {
  .validate_input_file(input)
  .validate_output_file(output, overwrite)
  .require_namespace("Seurat")
  if (utils::packageVersion("Seurat") < base::numeric_version("5.0.0")) {
    stop("Classic Assay -> Assay5 conversion requires Seurat >= 5.0.0.", call. = FALSE)
  }
  obj <- readRDS(input)
  if (!inherits(obj, "Seurat")) stop("Input is not a Seurat object.", call. = FALSE)
  if (!assay %in% SeuratObject::Assays(obj)) stop("Assay '", assay, "' was not found.", call. = FALSE)
  if (!inherits(obj[[assay]], "Assay") || inherits(obj[[assay]], "Assay5")) stop("Assay '", assay, "' is not a classic Assay.", call. = FALSE)
  if (new_assay %in% SeuratObject::Assays(obj)) stop("Assay '", new_assay, "' already exists.", call. = FALSE)

  converted <- methods::as(obj[[assay]], "Assay5")
  SeuratObject::Key(converted) <- paste0(gsub("[^A-Za-z0-9]", "", new_assay), "_")
  obj[[new_assay]] <- converted
  saveRDS(obj, output)
  message("Converted assay added as '", new_assay, "'. Saved: ", output)
  invisible(output)
}
