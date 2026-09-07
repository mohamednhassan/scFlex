test_that("convert_sc validates format arguments before dispatch", {
  expect_error(
    convert_sc("input", "unsupported", "sce", "output"),
    "Unsupported source format"
  )
  expect_error(
    convert_sc("input", "seurat", "unsupported", "output"),
    "Unsupported destination format"
  )
  expect_error(
    convert_sc("input", "sce", "sce", "output"),
    "Source and destination formats are the same"
  )
})
