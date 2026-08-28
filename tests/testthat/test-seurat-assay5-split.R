test_that("split Seurat v5 count layers are joined during extraction", {
  skip_if_not_installed("SeuratObject")

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

  expect_true(all(c("counts.sample1", "counts.sample2") %in% SeuratObject::Layers(seu[["RNA"]])))

  bundle <- scTransit:::.extract_seurat(seu)

  expect_equal(dim(bundle$counts), dim(counts))
  expect_identical(rownames(bundle$counts), rownames(counts))
  expect_identical(colnames(bundle$counts), colnames(counts))
  expect_equal(as.matrix(bundle$counts), as.matrix(counts))
  expect_true(any(grepl("Joined split Seurat count layers", bundle$notes, fixed = TRUE)))
})

test_that("counts-only bundles can populate AnnData X without normalization", {
  skip_if_not_installed("reticulate")

  counts <- Matrix::Matrix(matrix(c(1, 0, 2, 3), nrow = 2), sparse = TRUE)
  rownames(counts) <- c("g1", "g2")
  colnames(counts) <- c("c1", "c2")

  bundle <- scTransit:::.new_bundle(
    counts = counts,
    cell_metadata = data.frame(row.names = colnames(counts)),
    feature_metadata = data.frame(row.names = rownames(counts))
  )

  # This test is skipped if a usable Python AnnData environment is unavailable.
  ad_available <- tryCatch({
    scTransit:::.import_anndata()
    TRUE
  }, error = function(e) FALSE)
  skip_if_not(ad_available, "Python AnnData is unavailable")

  built <- scTransit:::.build_anndata(bundle, strict = TRUE)
  expect_true(any(grepl("raw counts were written to AnnData X", built$notes, fixed = TRUE)))
})
