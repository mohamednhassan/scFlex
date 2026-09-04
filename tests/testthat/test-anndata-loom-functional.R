# Functional validation for AnnData <-> Loom conversions

testthat::test_that("AnnData -> Loom preserves core matrix and annotations", {
  testthat::skip_if_not_installed("reticulate")

  reticulate::py_require(c("anndata>=0.10", "loompy>=3.0"))
  ad <- reticulate::import("anndata", convert = FALSE)
  np <- reticulate::import("numpy", convert = FALSE)
  pd <- reticulate::import("pandas", convert = FALSE)
  loompy <- reticulate::import("loompy", convert = FALSE)

  cells <- paste0("cell", 1:4)
  genes <- paste0("gene", 1:5)

  x <- matrix(
    c(
      0.1, 0.2, 0.3, 0.4, 0.5,
      1.1, 1.2, 1.3, 1.4, 1.5,
      2.1, 2.2, 2.3, 2.4, 2.5,
      3.1, 3.2, 3.3, 3.4, 3.5
    ),
    nrow = length(cells),
    byrow = TRUE
  )

  counts <- matrix(
    c(
      1, 0, 2, 0, 3,
      0, 4, 1, 0, 2,
      2, 1, 0, 3, 0,
      5, 0, 1, 2, 1
    ),
    nrow = length(cells),
    byrow = TRUE
  )

  obs <- pd$DataFrame(
    reticulate::dict(
      donor = c("D1", "D1", "D2", "D2"),
      score = c(1.1, 2.2, 3.3, 4.4)
    ),
    index = cells
  )

  var <- pd$DataFrame(
    reticulate::dict(
      gene_type = c("pc", "pc", "lnc", "pc", "lnc"),
      gc = c(0.41, 0.52, 0.63, 0.44, 0.55)
    ),
    index = genes
  )

  adata <- ad$AnnData(
    X = np$array(x),
    obs = obs,
    var = var
  )

  adata$layers$`__setitem__`(
    "counts",
    np$array(counts)$astype("float64")
  )

  pca <- matrix(
    c(
      0.1, 0.2,
      0.3, 0.4,
      0.5, 0.6,
      0.7, 0.8
    ),
    nrow = length(cells),
    byrow = TRUE
  )
  adata$obsm$`__setitem__`("X_pca", np$array(pca))

  input_file <- tempfile(fileext = ".h5ad")
  output_file <- tempfile(fileext = ".loom")
  adata$write_h5ad(input_file)

  result <- scTransit::convert_anndata_to_loom(
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
  testthat::expect_identical(shape, c(length(genes), length(cells)))

  full_slice <- reticulate::py_eval("slice(None)", convert = FALSE)
  loom_matrix <- reticulate::py_to_r(
    ds$`__getitem__`(
      reticulate::tuple(
        full_slice,
        full_slice
      )
    )
  )

  testthat::expect_equal(
    unname(as.matrix(loom_matrix)),
    unname(t(x)),
    tolerance = 1e-12
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

  gene_ids <- reticulate::py_to_r(ds$ra$`__getitem__`("Gene"))
  cell_ids <- reticulate::py_to_r(ds$ca$`__getitem__`("CellID"))

  testthat::expect_identical(as.character(gene_ids), genes)
  testthat::expect_identical(as.character(cell_ids), cells)

  testthat::expect_true("gene_type" %in% ra_keys)
  testthat::expect_true("gc" %in% ra_keys)
  testthat::expect_true("donor" %in% ca_keys)
  testthat::expect_true("score" %in% ca_keys)

  testthat::expect_identical(
    as.character(reticulate::py_to_r(ds$ra$`__getitem__`("gene_type"))),
    c("pc", "pc", "lnc", "pc", "lnc")
  )

  testthat::expect_equal(
    as.numeric(reticulate::py_to_r(ds$ra$`__getitem__`("gc"))),
    c(0.41, 0.52, 0.63, 0.44, 0.55)
  )

  testthat::expect_identical(
    as.character(reticulate::py_to_r(ds$ca$`__getitem__`("donor"))),
    c("D1", "D1", "D2", "D2")
  )

  testthat::expect_equal(
    as.numeric(reticulate::py_to_r(ds$ca$`__getitem__`("score"))),
    c(1.1, 2.2, 3.3, 4.4)
  )
})


testthat::test_that("Loom -> AnnData preserves primary matrix and axis names", {
  testthat::skip_if_not_installed("reticulate")

  reticulate::py_require(c("anndata>=0.10", "loompy>=3.0"))
  ad <- reticulate::import("anndata", convert = FALSE)
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
    nrow = length(genes),
    byrow = TRUE
  )

  row_attrs <- reticulate::dict(
    Gene = np$array(genes),
    gene_type = np$array(c("pc", "lnc", "pc", "pc"))
  )

  col_attrs <- reticulate::dict(
    CellID = np$array(cells),
    donor = np$array(c("A", "A", "B"))
  )

  input_file <- tempfile(fileext = ".loom")
  output_file <- tempfile(fileext = ".h5ad")

  loompy$create(
    input_file,
    np$array(mat)$astype("float64"),
    row_attrs = row_attrs,
    col_attrs = col_attrs
  )

  result <- scTransit::convert_loom_to_anndata(
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

  adata <- ad$read_h5ad(output_file)

  obs_names <- unlist(
    reticulate::iterate(adata$obs_names, as.character),
    use.names = FALSE
  )
  var_names <- unlist(
    reticulate::iterate(adata$var_names, as.character),
    use.names = FALSE
  )

  testthat::expect_identical(obs_names, cells)
  testthat::expect_identical(var_names, genes)

  x <- reticulate::py_to_r(adata$X)

  testthat::expect_equal(
    unname(as.matrix(x)),
    unname(t(mat)),
    tolerance = 0
  )

  obs <- reticulate::py_to_r(adata$obs)
  var <- reticulate::py_to_r(adata$var)

  testthat::expect_identical(
    as.character(obs$donor),
    c("A", "A", "B")
  )

  testthat::expect_identical(
    as.character(var$gene_type),
    c("pc", "lnc", "pc", "pc")
  )
})


testthat::test_that("AnnData -> Loom -> AnnData round trip preserves supported core content", {
  testthat::skip_if_not_installed("reticulate")

  reticulate::py_require(c("anndata>=0.10", "loompy>=3.0"))
  ad <- reticulate::import("anndata", convert = FALSE)
  np <- reticulate::import("numpy", convert = FALSE)
  pd <- reticulate::import("pandas", convert = FALSE)

  cells <- paste0("c", 1:4)
  genes <- paste0("g", 1:5)

  x <- matrix(
    seq_len(length(cells) * length(genes)) / 10,
    nrow = length(cells),
    ncol = length(genes)
  )

  obs <- pd$DataFrame(
    reticulate::dict(
      sample = c("S1", "S1", "S2", "S2")
    ),
    index = cells
  )

  var <- pd$DataFrame(
    reticulate::dict(
      symbol = paste0("SYM", 1:5)
    ),
    index = genes
  )

  adata <- ad$AnnData(
    X = np$array(x),
    obs = obs,
    var = var
  )

  input_file <- tempfile(fileext = ".h5ad")
  loom_file <- tempfile(fileext = ".loom")
  output_file <- tempfile(fileext = ".h5ad")

  adata$write_h5ad(input_file)

  scTransit::convert_anndata_to_loom(
    input = input_file,
    output = loom_file
  )

  scTransit::convert_loom_to_anndata(
    input = loom_file,
    output = output_file
  )

  roundtrip <- ad$read_h5ad(output_file)

  obs_names <- unlist(
    reticulate::iterate(roundtrip$obs_names, as.character),
    use.names = FALSE
  )
  var_names <- unlist(
    reticulate::iterate(roundtrip$var_names, as.character),
    use.names = FALSE
  )

  testthat::expect_identical(obs_names, cells)
  testthat::expect_identical(var_names, genes)

  x_roundtrip <- reticulate::py_to_r(roundtrip$X)

  testthat::expect_equal(
    unname(as.matrix(x_roundtrip)),
    unname(x),
    tolerance = 1e-7
  )

  obs_rt <- reticulate::py_to_r(roundtrip$obs)
  var_rt <- reticulate::py_to_r(roundtrip$var)

  testthat::expect_identical(
    as.character(obs_rt$sample),
    c("S1", "S1", "S2", "S2")
  )

  testthat::expect_identical(
    as.character(var_rt$symbol),
    paste0("SYM", 1:5)
  )
})


testthat::test_that("AnnData -> Loom rejects duplicate cell names", {
  testthat::skip_if_not_installed("reticulate")

  reticulate::py_require(c("anndata>=0.10", "loompy>=3.0"))
  ad <- reticulate::import("anndata", convert = FALSE)
  np <- reticulate::import("numpy", convert = FALSE)
  pd <- reticulate::import("pandas", convert = FALSE)

  cells <- c("cell1", "cell1", "cell2")
  genes <- c("gene1", "gene2")

  adata <- ad$AnnData(
    X = np$array(
      matrix(
        c(
          1, 2,
          3, 4,
          5, 6
        ),
        nrow = 3,
        byrow = TRUE
      )
    ),
    obs = pd$DataFrame(index = cells),
    var = pd$DataFrame(index = genes)
  )

  input_file <- tempfile(fileext = ".h5ad")
  output_file <- tempfile(fileext = ".loom")
  adata$write_h5ad(input_file)

  testthat::expect_error(
    scTransit::convert_anndata_to_loom(
      input = input_file,
      output = output_file
    ),
    regexp = "duplicate|Duplicate|unique|cell",
    ignore.case = TRUE
  )

  testthat::expect_false(file.exists(output_file))
})


testthat::test_that("duplicate feature names are made unique at Loom boundary", {
  testthat::skip_if_not_installed("reticulate")

  reticulate::py_require(c("anndata>=0.10", "loompy>=3.0"))
  ad <- reticulate::import("anndata", convert = FALSE)
  np <- reticulate::import("numpy", convert = FALSE)
  pd <- reticulate::import("pandas", convert = FALSE)
  loompy <- reticulate::import("loompy", convert = FALSE)

  cells <- c("cell1", "cell2")
  genes <- c("geneA", "geneA", "geneB")

  adata <- ad$AnnData(
    X = np$array(
      matrix(
        c(
          1, 2, 3,
          4, 5, 6
        ),
        nrow = 2,
        byrow = TRUE
      )
    ),
    obs = pd$DataFrame(index = cells),
    var = pd$DataFrame(index = genes)
  )

  input_file <- tempfile(fileext = ".h5ad")
  output_file <- tempfile(fileext = ".loom")
  adata$write_h5ad(input_file)

  scTransit::convert_anndata_to_loom(
    input = input_file,
    output = output_file
  )

  ds <- loompy$connect(output_file, mode = "r")
  on.exit(ds$close(), add = TRUE)

  gene_ids <- as.character(
    reticulate::py_to_r(ds$ra$`__getitem__`("Gene"))
  )

  testthat::expect_identical(length(gene_ids), length(genes))
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

    testthat::expect_identical(original_names, genes)
  }
})
