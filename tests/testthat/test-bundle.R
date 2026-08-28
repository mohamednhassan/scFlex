test_that("bundle validation catches misalignment", {
  counts <- Matrix::Matrix(matrix(1:6, nrow = 2), sparse = TRUE)
  rownames(counts) <- c("g1", "g2")
  colnames(counts) <- c("c1", "c2", "c3")

  data <- counts
  colnames(data) <- c("c3", "c2", "c1")

  x <- scTransit:::.new_bundle(counts = counts, data = data)
  expect_error(scTransit:::.validate_bundle(x), "not aligned")
})
