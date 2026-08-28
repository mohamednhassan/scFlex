\name{convert_sc}
\alias{convert_sc}
\alias{convert_seurat_to_anndata}
\alias{convert_anndata_to_seurat}
\alias{convert_seurat_to_sce}
\alias{convert_sce_to_seurat}
\alias{convert_sce_to_anndata}
\alias{convert_anndata_to_sce}
\alias{convert_seurat_to_loom}
\alias{convert_loom_to_seurat}
\alias{convert_sce_to_loom}
\alias{convert_loom_to_sce}
\alias{convert_anndata_to_loom}
\alias{convert_loom_to_anndata}
\title{Convert among supported single-cell formats}
\usage{
convert_sc(
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
)
}
\arguments{
  \item{input}{Input single-cell object or file path.}
  \item{output}{Output file path.}
  \item{from}{Source format. If \code{NULL}, scTransit attempts to detect it automatically.}
  \item{to}{Target format. Required when it cannot be inferred from the output extension.}
  \item{assay}{Seurat assay to use. Default is \code{"RNA"}.}
  \item{counts_layer}{Base name of the raw-count layer. Split Seurat v5 \code{counts.*} layers are detected and joined automatically.}
  \item{data_layer}{Base name of the normalized-data layer. Split Seurat v5 \code{data.*} layers are detected and joined automatically.}
  \item{sce_counts_assay}{Name of the SingleCellExperiment assay containing raw counts.}
  \item{sce_data_assay}{Name of the SingleCellExperiment assay containing normalized expression values.}
  \item{include_reductions}{Logical. Whether compatible dimensionality reductions should be transferred.}
  \item{strict}{Logical. Whether matrices required to construct the target object must be present. AnnData can use counts in \code{X}; counts-only sources are not normalized by scTransit.}
  \item{overwrite}{Logical. Whether an existing output file may be overwritten.}
}
\description{Converts among Seurat, SingleCellExperiment, AnnData, and Loom using a shared internal representation and alignment checks.}
\value{Invisibly returns an scTransit conversion result.}
