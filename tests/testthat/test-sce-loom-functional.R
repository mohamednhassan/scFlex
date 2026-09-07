# Functional validation for SingleCellExperiment <-> Loom conversions

testthat::test_that("SCE -> Loom preserves core matrix, names, metadata, and reductions where supported", {
  testthat::skip_if_not_installed("SingleCellExperiment")
  testthat::skip_if_not_installed("SummarizedExperiment")
  testthat::skip_if_not_installed("S4Vectors")
  testthat::skip_if_not_installed("reticulate")

  reticulate::py_require("loompy>=3.0")
  loompy <- reticulate::import("loompy", convert = FALSE)

  counts <- Matrix::Matrix(
    matrix(
      c(
        1, 0, 2, 0,
        0, 3, 0, 1,
        4, 0, 2, 0,
        0, 1, 0, 5,
        2, 2, 0, 1
      ),
      nrow = 5,
      byrow = TRUE
    ),
    sparse = TRUE
  )

  rownames(counts) <- paste0("gene", 1:5)
  colnames(counts) <- paste0("cell", 1:4)

  logcounts <- log1p(counts)

  sce <- SingleCellExperiment::SingleCellExperiment(
    assays = list(
      counts = counts,
      logcounts = logcounts
    ),
    colData = S4Vectors::DataFrame(
      donor = c("D1", "D1", "D2", "D2"),
      score = c(1.1, 2.2, 3.3, 4.4),
      row.names = colnames(counts)
    ),
    rowData = S4Vectors::DataFrame(
      symbol = paste0("SYM", 1:5),
      gene_type = c("pc", "pc", "lnc", "pc", "lnc"),
      row.names = rownames(counts)
    )
  )

  SingleCellExperiment::reducedDim(sce, "PCA") <- matrix(
    c(
      0.1, 0.2,
      0.3, 0.4,
      0.5, 0.6,
      0.7, 0.8
    ),
    nrow = 4,
    byrow = TRUE,
    dimnames = list(colnames(counts), c("PC1", "PC2"))
  )

  input_file <- tempfile(fileext = ".rds")
  output_file <- tempfile(fileext = ".loom")
  saveRDS(sce, input_file)

  result <- scFlex::convert_sce_to_loom(
    input = input_file,
    output = output_file
  )

  testthat::expect_true(file.exists(output_file))
  if (!is.null(result)) {
    testthat::expect_identical(
      normalizePath(result, mustWork = TRUE),
      normalizePath(output_file, mustWork = TRUE)
    )
  }

  ds <- loompy$connect(output_file, mode = "r")
  on.exit(ds$close(), add = TRUE)

  shape <- as.integer(reticulate::py_to_r(ds$shape))
  testthat::expect_identical(shape, c(5L, 4L))

  full_slice <- reticulate::py_eval("slice(None)", convert = FALSE)
  loom_matrix <- reticulate::py_to_r(
    ds$`__getitem__`(
      reticulate::tuple(full_slice, full_slice)
    )
  )

  testthat::expect_equal(
    unname(as.matrix(loom_matrix)),
    unname(as.matrix(logcounts)),
    tolerance = 1e-7
  )

  ra_keys <- unlist(
    reticulate::iterate(ds$ra$keys(), as.character),
    use.names = FALSE
  )
  ca_keys <- unlist(
    reticulate::iterate(ds$ca$keys(), as.character),
    use.names = FALSE
  )

  testthat::expect_true("Gene" %in% ra_keys)
  testthat::expect_true("CellID" %in% ca_keys)

  gene_ids <- as.character(
    reticulate::py_to_r(ds$ra$`__getitem__`("Gene"))
  )
  cell_ids <- as.character(
    reticulate::py_to_r(ds$ca$`__getitem__`("CellID"))
  )

  testthat::expect_identical(gene_ids, rownames(counts))
  testthat::expect_identical(cell_ids, colnames(counts))

  testthat::expect_true("symbol" %in% ra_keys)
  testthat::expect_true("gene_type" %in% ra_keys)
  testthat::expect_true("donor" %in% ca_keys)
  testthat::expect_true("score" %in% ca_keys)

  testthat::expect_identical(
    as.character(reticulate::py_to_r(ds$ra$`__getitem__`("symbol"))),
    paste0("SYM", 1:5)
  )

  testthat::expect_identical(
    as.character(reticulate::py_to_r(ds$ca$`__getitem__`("donor"))),
    c("D1", "D1", "D2", "D2")
  )

  layer_keys <- unlist(
    reticulate::iterate(ds$layers$keys(), as.character),
    use.names = FALSE
  )
  testthat::expect_true("counts" %in% layer_keys)

  counts_loom <- reticulate::py_to_r(
    ds$layers$`__getitem__`("counts")$`__getitem__`(
      reticulate::tuple(full_slice, full_slice)
    )
  )

  testthat::expect_equal(
    unname(as.matrix(counts_loom)),
    unname(as.matrix(counts)),
    tolerance = 0
  )

  pca_keys <- ca_keys[grepl("PCA|pca", ca_keys)]
  testthat::expect_true(length(pca_keys) >= 1)
})


testthat::test_that("Loom -> SCE preserves primary matrix, names, and metadata", {
  testthat::skip_if_not_installed("SingleCellExperiment")
  testthat::skip_if_not_installed("SummarizedExperiment")
  testthat::skip_if_not_installed("reticulate")

  reticulate::py_require("loompy>=3.0")
  np <- reticulate::import("numpy", convert = FALSE)
  loompy <- reticulate::import("loompy", convert = FALSE)

  genes <- paste0("gene", 1:4)
  cells <- paste0("cell", 1:3)

  mat <- matrix(
    c(
      1, 0, 2,
      0, 3, 1,
      4, 0, 0,
      2, 1, 5
    ),
    nrow = 4,
    byrow = TRUE
  )

  row_attrs <- reticulate::dict(
    Gene = np$array(genes),
    symbol = np$array(paste0("SYM", 1:4))
  )

  col_attrs <- reticulate::dict(
    CellID = np$array(cells),
    donor = np$array(c("A", "A", "B"))
  )

  input_file <- tempfile(fileext = ".loom")
  output_file <- tempfile(fileext = ".rds")

  loompy$create(
    input_file,
    np$array(mat)$astype("float64"),
    row_attrs = row_attrs,
    col_attrs = col_attrs
  )

  result <- scFlex::convert_loom_to_sce(
    input = input_file,
    output = output_file
  )

  testthat::expect_true(file.exists(output_file))
  if (!is.null(result)) {
    testthat::expect_identical(
      normalizePath(result, mustWork = TRUE),
      normalizePath(output_file, mustWork = TRUE)
    )
  }

  sce <- readRDS(output_file)

  testthat::expect_s4_class(sce, "SingleCellExperiment")
  testthat::expect_equal(dim(sce), c(4, 3))
  testthat::expect_identical(colnames(sce), cells)
  testthat::expect_identical(rownames(sce), genes)

  assay_names <- SummarizedExperiment::assayNames(sce)
  testthat::expect_true("counts" %in% assay_names)

  counts_rt <- SummarizedExperiment::assay(sce, "counts")

  testthat::expect_equal(
    unname(as.matrix(counts_rt)),
    unname(mat),
    tolerance = 0
  )

  cd <- as.data.frame(SummarizedExperiment::colData(sce))
  rd <- as.data.frame(SummarizedExperiment::rowData(sce))

  testthat::expect_identical(
    as.character(cd$donor),
    c("A", "A", "B")
  )

  if ("symbol" %in% colnames(rd)) {
    testthat::expect_identical(
      as.character(rd$symbol),
      paste0("SYM", 1:4)
    )
  }
})


testthat::test_that("SCE -> Loom -> SCE round trip preserves supported core content", {
  testthat::skip_if_not_installed("SingleCellExperiment")
  testthat::skip_if_not_installed("SummarizedExperiment")
  testthat::skip_if_not_installed("S4Vectors")
  testthat::skip_if_not_installed("reticulate")

  reticulate::py_require("loompy>=3.0")

  counts <- Matrix::Matrix(
    matrix(
      c(
        1, 0, 2,
        0, 3, 1,
        4, 0, 0,
        2, 1, 5
      ),
      nrow = 4,
      byrow = TRUE
    ),
    sparse = TRUE
  )
  rownames(counts) <- paste0("g", 1:4)
  colnames(counts) <- paste0("c", 1:3)

  sce <- SingleCellExperiment::SingleCellExperiment(
    assays = list(counts = counts),
    colData = S4Vectors::DataFrame(
      sample = c("S1", "S1", "S2"),
      row.names = colnames(counts)
    )
  )

  input_file <- tempfile(fileext = ".rds")
  loom_file <- tempfile(fileext = ".loom")
  output_file <- tempfile(fileext = ".rds")

  saveRDS(sce, input_file)

  scFlex::convert_sce_to_loom(
    input = input_file,
    output = loom_file
  )
  scFlex::convert_loom_to_sce(
    input = loom_file,
    output = output_file
  )

  rt <- readRDS(output_file)

  testthat::expect_identical(colnames(rt), colnames(sce))
  testthat::expect_identical(rownames(rt), rownames(sce))

  counts_rt <- SummarizedExperiment::assay(rt, "counts")

  testthat::expect_equal(
    unname(as.matrix(counts_rt)),
    unname(as.matrix(counts)),
    tolerance = 0
  )

  cd <- as.data.frame(SummarizedExperiment::colData(rt))
  testthat::expect_identical(
    as.character(cd$sample),
    c("S1", "S1", "S2")
  )
})


testthat::test_that("counts-only SCE converts to Loom without inventing normalized data", {
  testthat::skip_if_not_installed("SingleCellExperiment")
  testthat::skip_if_not_installed("reticulate")

  reticulate::py_require("loompy>=3.0")
  loompy <- reticulate::import("loompy", convert = FALSE)

  counts <- Matrix::Matrix(
    matrix(
      c(
        1, 0,
        2, 3,
        0, 4
      ),
      nrow = 3,
      byrow = TRUE
    ),
    sparse = TRUE
  )
  rownames(counts) <- paste0("g", 1:3)
  colnames(counts) <- paste0("c", 1:2)

  sce <- SingleCellExperiment::SingleCellExperiment(
    assays = list(counts = counts)
  )

  input_file <- tempfile(fileext = ".rds")
  output_file <- tempfile(fileext = ".loom")
  saveRDS(sce, input_file)

  scFlex::convert_sce_to_loom(
    input = input_file,
    output = output_file
  )

  testthat::expect_true(file.exists(output_file))

  ds <- loompy$connect(output_file, mode = "r")
  on.exit(ds$close(), add = TRUE)

  full_slice <- reticulate::py_eval("slice(None)", convert = FALSE)
  primary <- reticulate::py_to_r(
    ds$`__getitem__`(
      reticulate::tuple(full_slice, full_slice)
    )
  )

  testthat::expect_equal(
    unname(as.matrix(primary)),
    unname(as.matrix(counts)),
    tolerance = 0
  )
})


testthat::test_that("SCE duplicate feature names are made unique at Loom boundary", {
  testthat::skip_if_not_installed("SingleCellExperiment")
  testthat::skip_if_not_installed("reticulate")

  reticulate::py_require("loompy>=3.0")
  loompy <- reticulate::import("loompy", convert = FALSE)

  counts <- Matrix::Matrix(
    matrix(
      c(
        1, 0,
        2, 3,
        0, 4
      ),
      nrow = 3,
      byrow = TRUE
    ),
    sparse = TRUE
  )
  rownames(counts) <- c("geneA", "geneA", "geneB")
  colnames(counts) <- c("cell1", "cell2")

  sce <- SingleCellExperiment::SingleCellExperiment(
    assays = list(counts = counts)
  )

  input_file <- tempfile(fileext = ".rds")
  output_file <- tempfile(fileext = ".loom")
  saveRDS(sce, input_file)

  scFlex::convert_sce_to_loom(
    input = input_file,
    output = output_file
  )

  ds <- loompy$connect(output_file, mode = "r")
  on.exit(ds$close(), add = TRUE)

  gene_ids <- as.character(
    reticulate::py_to_r(ds$ra$`__getitem__`("Gene"))
  )

  testthat::expect_identical(length(gene_ids), 3L)
  testthat::expect_false(anyDuplicated(gene_ids) > 0)

  ra_keys <- unlist(
    reticulate::iterate(ds$ra$keys(), as.character),
    use.names = FALSE
  )

  if ("original_feature_name" %in% ra_keys) {
    original_names <- as.character(
      reticulate::py_to_r(
        ds$ra$`__getitem__`("original_feature_name")
      )
    )

    testthat::expect_identical(
      original_names,
      c("geneA", "geneA", "geneB")
    )
  }
})
