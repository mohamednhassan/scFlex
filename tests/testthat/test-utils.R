test_that("internal formatting helper handles empty and populated input", {
  expect_identical(scFlex:::collapse_or_none(character(0)), "None")
  expect_identical(scFlex:::collapse_or_none(c("a", "b")), "a, b")
})

test_that("inspect_sc rejects unsupported file extensions", {
  x <- tempfile(fileext = ".txt")
  file.create(x)
  expect_error(inspect_sc(x), "Unsupported file extension")
})
