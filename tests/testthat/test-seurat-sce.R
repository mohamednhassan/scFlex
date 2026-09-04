test_that("Seurat to SCE preserves counts and cell metadata", {
  skip_if_not_installed("SeuratObject")
  skip_if_not_installed("SingleCellExperiment")
  skip_if_not_installed("SummarizedExperiment")

  counts <- Matrix::Matrix(
    matrix(c(1, 0, 2, 3, 0, 1), nrow = 2),
    sparse = TRUE
  )
  rownames(counts) <- c("g1", "g2")
  colnames(counts) <- c("c1", "c2", "c3")

  seu <- SeuratObject::CreateSeuratObject(counts = counts)
  seu$group <- c("A", "A", "B")

  input <- tempfile(fileext = ".rds")
  output <- tempfile(fileext = ".rds")
  saveRDS(seu, input)

  expect_invisible(convert_seurat_to_sce(input, output))
  expect_true(file.exists(output))

  sce <- readRDS(output)
  expect_s4_class(sce, "SingleCellExperiment")
  expect_equal(
    as.matrix(SummarizedExperiment::assay(sce, "counts")),
    as.matrix(counts)
  )
  expect_identical(
    as.character(SummarizedExperiment::colData(sce)$group),
    c("A", "A", "B")
  )
})
