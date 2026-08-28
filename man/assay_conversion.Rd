\name{convert_seu_v5_to_classic}
\alias{convert_seu_v5_to_classic}
\alias{convert_seu_classic_to_v5}
\title{Convert Seurat assay structure}
\usage{
convert_seu_v5_to_classic(
  input,
  output,
  assay = "RNA",
  new_assay = "RNA_classic",
  overwrite = FALSE
)

convert_seu_classic_to_v5(
  input,
  output,
  assay = "RNA",
  new_assay = "RNA_assay5",
  overwrite = FALSE
)
}
\arguments{
  \item{input}{Path to the input Seurat object saved as an RDS file.}
  \item{output}{Path where the converted Seurat object should be saved.}
  \item{assay}{Name of the Seurat assay to convert. Default is \code{"RNA"}.}
  \item{new_assay}{Name to assign to the converted assay in the output object.}
  \item{overwrite}{Logical. Whether an existing output file may be overwritten.}
}
\description{Converts an individual Seurat assay between classic \code{Assay} and Seurat v5 \code{Assay5} representations.}
\value{Invisibly returns the output path.}
