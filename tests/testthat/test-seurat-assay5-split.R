test_that("split Seurat Assay5 count layers convert to SCE", {
  skip_if_not_installed("SeuratObject")
  skip_if_not_installed("SingleCellExperiment")
  skip_if_not_installed("SummarizedExperiment")

  counts <- Matrix::Matrix(
    matrix(c(1, 0, 2, 3, 0, 1, 4, 0, 2, 1, 0, 5), nrow = 3),
    sparse = TRUE
  )
  rownames(counts) <- paste0("g", 1:3)
  colnames(counts) <- paste0("c", 1:4)

  seu <- SeuratObject::CreateSeuratObject(counts = counts)
  skip_if_not(inherits(seu[["RNA"]], "Assay5"), "Requires SeuratObject Assay5")

  seu$sample <- c("sample1", "sample1", "sample2", "sample2")
  seu[["RNA"]] <- split(seu[["RNA"]], f = seu$sample)

  expect_true(
    all(c("counts.sample1", "counts.sample2") %in% SeuratObject::Layers(seu[["RNA"]]))
  )

  input <- tempfile(fileext = ".rds")
  output <- tempfile(fileext = ".rds")
  saveRDS(seu, input)

  expect_invisible(convert_seurat_to_sce(input, output))
  sce <- readRDS(output)

  expect_equal(
    as.matrix(SummarizedExperiment::assay(sce, "counts")),
    as.matrix(counts)
  )
})
