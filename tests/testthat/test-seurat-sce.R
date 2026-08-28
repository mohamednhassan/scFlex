test_that("Seurat -> SCE preserves matrices and metadata", {
  skip_if_not_installed("SeuratObject")
  skip_if_not_installed("SingleCellExperiment")
  skip_if_not_installed("SummarizedExperiment")

  counts <- Matrix::Matrix(matrix(c(1,0,2,3,0,1), nrow = 2), sparse = TRUE)
  rownames(counts) <- c("g1", "g2")
  colnames(counts) <- c("c1", "c2", "c3")

  seu <- SeuratObject::CreateSeuratObject(counts = counts)
  if (inherits(seu[["RNA"]], "Assay5")) {
    SeuratObject::LayerData(seu[["RNA"]], "data") <- log1p(counts)
  } else {
    seu <- SeuratObject::SetAssayData(seu, assay = "RNA", slot = "data", new.data = log1p(counts))
  }

  bundle <- scTransit:::.extract_seurat(seu)
  sce <- scTransit:::.build_sce(bundle)

  expect_equal(dim(sce), dim(seu))
  expect_equal(
    as.matrix(SummarizedExperiment::assay(sce, "counts")),
    as.matrix(counts)
  )
})
