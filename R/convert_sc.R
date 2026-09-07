#' Convert between supported single-cell formats
#'
#' Dispatches a file conversion between Seurat, SingleCellExperiment, AnnData, and Loom using the corresponding scFlex conversion function.
#'
#' @details
#' Source and destination must differ. Unsupported format names or unimplemented conversion paths produce an error.
#'
#' @param input Character string giving the path to the input file.
#' @param source Character string naming the source format. Supported values are `"seurat"`, `"sce"`, `"anndata"`, and `"loom"`. Matching is case-insensitive.
#' @param destination Character string naming the destination format. Supported values are `"seurat"`, `"sce"`, `"anndata"`, and `"loom"`. Matching is case-insensitive.
#' @param output Character string giving the output file path.
#'
#' @return The result returned by the selected conversion function, invisibly. Current converters normally return the output path invisibly; some conversions may return `NULL` invisibly when required count data cannot be identified.
#' @export
convert_sc <- function(input, source, destination, output) {
  valid_formats <- c(
    "seurat",
    "sce",
    "anndata",
    "loom"
  )
  
  source <- tolower(source)
  destination <- tolower(destination)
  
  if (!source %in% valid_formats) {
    stop(
      "Unsupported source format: ",
      source,
      "\nSupported formats: ",
      paste(valid_formats, collapse = ", "),
      call. = FALSE
    )
  }
  
  if (!destination %in% valid_formats) {
    stop(
      "Unsupported destination format: ",
      destination,
      "\nSupported formats: ",
      paste(valid_formats, collapse = ", "),
      call. = FALSE
    )
  }
  
  if (source == destination) {
    stop(
      "Source and destination formats are the same: ",
      source,
      call. = FALSE
    )
  }
  
  conversion_key <- paste(
    source,
    destination,
    sep = "_to_"
  )
  
  message("Converting ", source, " -> ", destination, "...")
  
  result <- switch(
    conversion_key,
    
    "seurat_to_anndata" =
      convert_seurat_to_anndata(
        input = input,
        output = output
      ),
    
    "anndata_to_seurat" =
      convert_anndata_to_seurat(
        input = input,
        output = output
      ),
    
    "seurat_to_sce" =
      convert_seurat_to_sce(
        input = input,
        output = output
      ),
    
    "sce_to_seurat" =
      convert_sce_to_seurat(
        input = input,
        output = output
      ),
    
    "sce_to_anndata" =
      convert_sce_to_anndata(
        input = input,
        output = output
      ),
    
    "anndata_to_sce" =
      convert_anndata_to_sce(
        input = input,
        output = output
      ),
    
    "seurat_to_loom" =
      convert_seurat_to_loom(
        input = input,
        output = output
      ),
    
    "loom_to_seurat" =
      convert_loom_to_seurat(
        input = input,
        output = output
      ),
    
    "sce_to_loom" =
      convert_sce_to_loom(
        input = input,
        output = output
      ),
    
    "loom_to_sce" =
      convert_loom_to_sce(
        input = input,
        output = output
      ),
    
    "anndata_to_loom" =
      convert_anndata_to_loom(
        input = input,
        output = output
      ),
    
    "loom_to_anndata" =
      convert_loom_to_anndata(
        input = input,
        output = output
      ),
    
    stop(
      "Conversion path not implemented: ",
      source,
      " -> ",
      destination,
      call. = FALSE
    )
  )
  
  invisible(result)
}
