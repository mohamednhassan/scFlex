test_that("output validation prevents accidental overwrite", {
  x <- tempfile(fileext = ".rds")
  saveRDS(1, x)
  expect_error(scTransit:::.validate_output_file(x, overwrite = FALSE), "already exists")
  expect_invisible(scTransit:::.validate_output_file(x, overwrite = TRUE))
})

test_that("format detection recognizes file extensions", {
  x <- tempfile(fileext = ".h5ad")
  file.create(x)
  expect_equal(scTransit:::.detect_format(x), "anndata")

  y <- tempfile(fileext = ".loom")
  file.create(y)
  expect_equal(scTransit:::.detect_format(y), "loom")
})
