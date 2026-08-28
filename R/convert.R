#' Convert among Seurat, SingleCellExperiment, AnnData, and Loom
#'
#' `convert_sc()` is the main scTransit interface. It extracts a source object
#' into a common internal representation and constructs the requested target.
#' Unsupported optional components are skipped with informative messages;
#' required expression data trigger explicit errors rather than silent guessing.
#'
#' @param input Input file path. Seurat and SCE inputs should be `.rds`; AnnData
#'   inputs `.h5ad`; Loom inputs `.loom`.
#' @param output Output file path.
#' @param from Optional source format. One of `"seurat"`, `"sce"`, `"anndata"`, `"loom"`.
#' @param to Optional target format. Required for `.rds` outputs because both
#'   Seurat and SCE use RDS.
#' @param assay Seurat assay name. Default `"RNA"`.
#' @param counts_layer Raw-count layer/assay name. Default `"counts"`.
#' @param data_layer Normalized-data layer/assay name. Default `"data"` for
#'   Seurat sources; SCE normalized assay is controlled by `sce_data_assay`.
#' @param sce_counts_assay SCE raw-count assay name. Default `"counts"`.
#' @param sce_data_assay SCE normalized assay name. Default `"logcounts"`.
#' @param include_reductions Preserve compatible dimensional reductions.
#' @param strict If `TRUE`, missing matrices that are required to construct the
#'   target object stop conversion. AnnData can validly use raw counts in `X`;
#'   when normalized data are absent, counts are therefore written to `X`
#'   without normalization and this is reported in the conversion notes.
#' @param overwrite Replace an existing output file.
#' @return Invisibly returns an `scTransit_result` containing the output path,
#'   source and target formats, and conversion notes.
#' @export
convert_sc <- function(
  input,
  output,
  from = NULL,
  to = NULL,
  assay = "RNA",
  counts_layer = "counts",
  data_layer = "data",
  sce_counts_assay = "counts",
  sce_data_assay = "logcounts",
  include_reductions = TRUE,
  strict = TRUE,
  overwrite = FALSE
) {
  .validate_input_file(input)
  .validate_output_file(output, overwrite)
  from <- .detect_format(input, from)
  to <- .target_from_extension(output, to)

  if (identical(from, to)) {
    stop("Source and target formats are both '", from, "'. No conversion is required.", call. = FALSE)
  }

  message("scTransit: ", from, " -> ", to)

  bundle <- switch(
    from,
    seurat = .extract_seurat(
      readRDS(input), assay = assay, counts_layer = counts_layer,
      data_layer = data_layer, include_reductions = include_reductions
    ),
    sce = .extract_sce(
      readRDS(input), counts_assay = sce_counts_assay,
      data_assay = sce_data_assay, include_reductions = include_reductions
    ),
    anndata = .extract_anndata(
      .read_anndata_file(input), counts_layer = counts_layer,
      include_reductions = include_reductions, source_format = "anndata"
    ),
    loom = .extract_anndata(
      .read_loom_file(input), counts_layer = counts_layer,
      include_reductions = include_reductions, source_format = "loom"
    )
  )

  notes <- bundle$notes
  if (from == "loom" || to == "loom") {
    notes <- c(notes, "Loom supports a smaller data model than Seurat, SCE, and AnnData; some components may not be representable.")
  }

  if (to == "seurat") {
    out_obj <- .build_seurat(bundle, assay = assay, strict = strict)
    saveRDS(out_obj, output)
  } else if (to == "sce") {
    out_obj <- .build_sce(bundle, counts_assay = sce_counts_assay, data_assay = sce_data_assay)
    saveRDS(out_obj, output)
  } else if (to == "anndata") {
    built <- .build_anndata(bundle, counts_layer = counts_layer, strict = strict, include_reductions = include_reductions)
    notes <- c(notes, built$notes)
    .write_anndata_file(built$adata, output)
  } else if (to == "loom") {
    built <- .build_anndata(bundle, counts_layer = counts_layer, strict = strict, include_reductions = include_reductions)
    notes <- c(notes, built$notes)
    .write_loom_file(built$adata, output, include_reductions = include_reductions)
  }

  notes <- unique(notes[nzchar(notes)])
  message("Conversion completed: ", output)
  if (length(notes)) {
    message("Notes:")
    for (note in notes) message("  - ", note)
  }
  invisible(.conversion_result(output, from, to, notes))
}

#' @rdname convert_sc
#' @export
convert_seurat_to_anndata <- function(input, output, ...) convert_sc(input, output, from = "seurat", to = "anndata", ...)
#' @rdname convert_sc
#' @export
convert_anndata_to_seurat <- function(input, output, ...) convert_sc(input, output, from = "anndata", to = "seurat", ...)
#' @rdname convert_sc
#' @export
convert_seurat_to_sce <- function(input, output, ...) convert_sc(input, output, from = "seurat", to = "sce", ...)
#' @rdname convert_sc
#' @export
convert_sce_to_seurat <- function(input, output, ...) convert_sc(input, output, from = "sce", to = "seurat", ...)
#' @rdname convert_sc
#' @export
convert_sce_to_anndata <- function(input, output, ...) convert_sc(input, output, from = "sce", to = "anndata", ...)
#' @rdname convert_sc
#' @export
convert_anndata_to_sce <- function(input, output, ...) convert_sc(input, output, from = "anndata", to = "sce", ...)
#' @rdname convert_sc
#' @export
convert_seurat_to_loom <- function(input, output, ...) convert_sc(input, output, from = "seurat", to = "loom", ...)
#' @rdname convert_sc
#' @export
convert_loom_to_seurat <- function(input, output, ...) convert_sc(input, output, from = "loom", to = "seurat", ...)
#' @rdname convert_sc
#' @export
convert_sce_to_loom <- function(input, output, ...) convert_sc(input, output, from = "sce", to = "loom", ...)
#' @rdname convert_sc
#' @export
convert_loom_to_sce <- function(input, output, ...) convert_sc(input, output, from = "loom", to = "sce", ...)
#' @rdname convert_sc
#' @export
convert_anndata_to_loom <- function(input, output, ...) convert_sc(input, output, from = "anndata", to = "loom", ...)
#' @rdname convert_sc
#' @export
convert_loom_to_anndata <- function(input, output, ...) convert_sc(input, output, from = "loom", to = "anndata", ...)
