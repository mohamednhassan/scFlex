\name{inspect_sc}
\alias{inspect_sc}
\alias{check_input_info}
\title{Inspect a supported single-cell object}
\usage{
inspect_sc(
  input,
  format = NULL,
  assay = "RNA",
  counts_layer = "counts",
  data_layer = "data"
)

check_input_info(path_to_file)
}
\arguments{
  \item{input}{Input single-cell object file path.}
  \item{format}{Optional source format. If \code{NULL}, scTransit attempts to detect it automatically.}
  \item{assay}{Seurat assay to inspect. Default is \code{"RNA"}.}
  \item{counts_layer}{Base name of the layer containing raw counts. Split \code{counts.*} layers are also detected.}
  \item{data_layer}{Base name of the layer containing normalized expression values. Split \code{data.*} layers are also detected.}
  \item{path_to_file}{Path to a Seurat object saved as an RDS file, used by \code{check_input_info()}.}
}
\description{Reports major components of Seurat, SingleCellExperiment, AnnData, or Loom input without modifying the source. For Seurat v5 Assay5 inputs, canonical and split count/data layers are reported separately.}
\value{Invisibly returns a named list.}
