# Functional validation for SingleCellExperiment <-> AnnData conversions

testthat::test_that("SCE -> AnnData preserves core data components", {
  testthat::skip_on_cran()
  testthat::skip_if_not_installed("reticulate")
  testthat::skip_if_not_installed("SingleCellExperiment")
  testthat::skip_if_not_installed("SummarizedExperiment")
  testthat::skip_if_not_installed("S4Vectors")
  testthat::skip_if_not_installed("Matrix")

  reticulate::py_require("anndata>=0.10")
  ad <- reticulate::import("anndata", convert = FALSE)

  genes <- paste0("gene", 1:6)
  cells <- paste0("cell", 1:5)

  counts_dense <- matrix(
    c(
      1, 0, 3, 0, 2,
      0, 4, 0, 1, 0,
      5, 0, 2, 0, 1,
      0, 2, 0, 3, 0,
      1, 1, 0, 0, 4,
      0, 0, 2, 5, 1
    ),
    nrow = length(genes),
    byrow = TRUE,
    dimnames = list(genes, cells)
  )

  counts <- Matrix::Matrix(counts_dense, sparse = TRUE)
  logcounts <- Matrix::Matrix(log1p(counts_dense), sparse = TRUE)

  metadata <- data.frame(
    donor = c("D1", "D1", "D2", "D2", "D3"),
    condition = factor(c("Lean", "Lean", "Obese", "Obese", "Lean")),
    score = c(1.2, 2.4, 3.6, 4.8, 6.0),
    row.names = cells
  )

  feature_metadata <- data.frame(
    gene_type = c("protein_coding", "protein_coding", "lncRNA",
                  "protein_coding", "lncRNA", "protein_coding"),
    gc = c(0.41, 0.52, 0.63, 0.44, 0.55, 0.66),
    row.names = genes
  )

  pca <- matrix(
    c(
      0.1, 0.5,
      0.2, 0.4,
      0.3, 0.3,
      0.4, 0.2,
      0.5, 0.1
    ),
    nrow = length(cells),
    byrow = TRUE,
    dimnames = list(cells, c("PC_1", "PC_2"))
  )

  sce <- SingleCellExperiment::SingleCellExperiment(
    assays = list(
      counts = counts,
      logcounts = logcounts
    ),
    colData = S4Vectors::DataFrame(metadata),
    rowData = S4Vectors::DataFrame(feature_metadata)
  )
  SingleCellExperiment::reducedDim(sce, "pca") <- pca

  input_file <- tempfile(fileext = ".rds")
  output_file <- tempfile(fileext = ".h5ad")
  saveRDS(sce, input_file)

  result <- scFlex::convert_sce_to_anndata(
    input = input_file,
    output = output_file
  )

  testthat::expect_true(file.exists(output_file))
  testthat::expect_identical(
    normalizePath(result, mustWork = TRUE),
    normalizePath(output_file, mustWork = TRUE)
  )

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
  testthat::expect_equal(
    as.integer(reticulate::py_to_r(adata$n_obs)),
    length(cells)
  )
  testthat::expect_equal(
    as.integer(reticulate::py_to_r(adata$n_vars)),
    length(genes)
  )

  x <- reticulate::py_to_r(adata$X)
  counts_layer <- reticulate::py_to_r(
    adata$layers$`__getitem__`("counts")
  )
  logcounts_layer <- reticulate::py_to_r(
    adata$layers$`__getitem__`("logcounts")
  )

  testthat::expect_equal(
    unname(as.matrix(x)),
    unname(t(as.matrix(logcounts))),
    tolerance = 1e-12
  )
  testthat::expect_equal(
    unname(as.matrix(counts_layer)),
    unname(t(as.matrix(counts))),
    tolerance = 0
  )
  testthat::expect_equal(
    unname(as.matrix(logcounts_layer)),
    unname(t(as.matrix(logcounts))),
    tolerance = 1e-12
  )

  obs <- reticulate::py_to_r(adata$obs)
  testthat::expect_identical(rownames(obs), cells)
  testthat::expect_identical(as.character(obs$donor), metadata$donor)
  testthat::expect_identical(
    as.character(obs$condition),
    as.character(metadata$condition)
  )
  testthat::expect_equal(as.numeric(obs$score), metadata$score)

  var <- reticulate::py_to_r(adata$var)
  testthat::expect_identical(rownames(var), genes)
  testthat::expect_identical(
    as.character(var$gene_type),
    feature_metadata$gene_type
  )
  testthat::expect_equal(as.numeric(var$gc), feature_metadata$gc)

  pca_out <- reticulate::py_to_r(
    adata$obsm$`__getitem__`("X_pca")
  )
  testthat::expect_equal(
    as.matrix(pca_out),
    unname(pca),
    tolerance = 0
  )
})


testthat::test_that("AnnData -> SCE prefers layers['counts'] and preserves core content", {
  testthat::skip_on_cran()
  testthat::skip_if_not_installed("reticulate")
  testthat::skip_if_not_installed("SingleCellExperiment")
  testthat::skip_if_not_installed("SummarizedExperiment")
  testthat::skip_if_not_installed("Matrix")

  reticulate::py_require("anndata>=0.10")
  ad <- reticulate::import("anndata", convert = FALSE)
  np <- reticulate::import("numpy", convert = FALSE)
  pd <- reticulate::import("pandas", convert = FALSE)

  genes <- paste0("gene", 1:5)
  cells <- paste0("cell", 1:4)

  counts <- matrix(
    c(
      1, 0, 2, 3, 0,
      0, 4, 1, 0, 2,
      2, 2, 0, 1, 1,
      5, 0, 1, 2, 0
    ),
    nrow = length(cells),
    byrow = TRUE
  )
  logcounts <- log1p(counts) + 0.125

  obs <- pd$DataFrame(
    reticulate::dict(
      donor = c("A", "A", "B", "B"),
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
    X = np$array(logcounts),
    obs = obs,
    var = var
  )

  adata$layers$`__setitem__`(
    "counts",
    np$array(counts)$astype("float64")
  )
  adata$layers$`__setitem__`(
    "logcounts",
    np$array(logcounts)
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
  output_file <- tempfile(fileext = ".rds")
  adata$write_h5ad(input_file)

  result <- scFlex::convert_anndata_to_sce(
    input = input_file,
    output = output_file
  )

  testthat::expect_true(file.exists(output_file))
  testthat::expect_identical(
    normalizePath(result, mustWork = TRUE),
    normalizePath(output_file, mustWork = TRUE)
  )

  sce <- readRDS(output_file)

  testthat::expect_s4_class(sce, "SingleCellExperiment")
  testthat::expect_identical(rownames(sce), genes)
  testthat::expect_identical(colnames(sce), cells)

  sce_counts <- SummarizedExperiment::assay(sce, "counts")
  sce_logcounts <- SummarizedExperiment::assay(sce, "logcounts")

  testthat::expect_equal(
    unname(as.matrix(sce_counts)),
    unname(t(counts)),
    tolerance = 0
  )
  testthat::expect_equal(
    unname(as.matrix(sce_logcounts)),
    unname(t(logcounts)),
    tolerance = 1e-12
  )

  cell_meta <- as.data.frame(SummarizedExperiment::colData(sce))
  testthat::expect_identical(rownames(cell_meta), cells)
  testthat::expect_identical(as.character(cell_meta$donor), c("A", "A", "B", "B"))
  testthat::expect_equal(as.numeric(cell_meta$score), c(1.1, 2.2, 3.3, 4.4))

  feature_meta <- as.data.frame(SummarizedExperiment::rowData(sce))
  testthat::expect_identical(rownames(feature_meta), genes)
  testthat::expect_identical(
    as.character(feature_meta$gene_type),
    c("pc", "pc", "lnc", "pc", "lnc")
  )
  testthat::expect_equal(
    as.numeric(feature_meta$gc),
    c(0.41, 0.52, 0.63, 0.44, 0.55)
  )

  testthat::expect_true(
    "X_pca" %in% SingleCellExperiment::reducedDimNames(sce) ||
      "pca" %in% SingleCellExperiment::reducedDimNames(sce)
  )

  reduction_name <- if (
    "X_pca" %in% SingleCellExperiment::reducedDimNames(sce)
  ) {
    "X_pca"
  } else {
    "pca"
  }

  testthat::expect_equal(
    unname(as.matrix(SingleCellExperiment::reducedDim(sce, reduction_name))),
    unname(pca),
    tolerance = 0
  )
})


testthat::test_that("SCE -> AnnData -> SCE round trip preserves supported content", {
  testthat::skip_on_cran()
  testthat::skip_if_not_installed("reticulate")
  testthat::skip_if_not_installed("SingleCellExperiment")
  testthat::skip_if_not_installed("SummarizedExperiment")
  testthat::skip_if_not_installed("S4Vectors")
  testthat::skip_if_not_installed("Matrix")

  reticulate::py_require("anndata>=0.10")

  genes <- paste0("g", 1:5)
  cells <- paste0("c", 1:4)

  counts_dense <- matrix(
    c(
      1, 0, 2, 3,
      0, 4, 1, 0,
      2, 2, 0, 1,
      5, 0, 1, 2,
      0, 3, 2, 1
    ),
    nrow = length(genes),
    byrow = TRUE,
    dimnames = list(genes, cells)
  )

  counts <- Matrix::Matrix(counts_dense, sparse = TRUE)
  logcounts <- Matrix::Matrix(log1p(counts_dense), sparse = TRUE)

  metadata <- data.frame(
    sample = c("S1", "S1", "S2", "S2"),
    condition = c("A", "A", "B", "B"),
    row.names = cells
  )

  feature_metadata <- data.frame(
    symbol = paste0("SYM", 1:5),
    row.names = genes
  )

  pca <- matrix(
    seq_len(length(cells) * 2) / 10,
    nrow = length(cells),
    ncol = 2,
    dimnames = list(cells, c("PC_1", "PC_2"))
  )

  original <- SingleCellExperiment::SingleCellExperiment(
    assays = list(
      counts = counts,
      logcounts = logcounts
    ),
    colData = S4Vectors::DataFrame(metadata),
    rowData = S4Vectors::DataFrame(feature_metadata)
  )
  SingleCellExperiment::reducedDim(original, "pca") <- pca

  sce_file <- tempfile(fileext = ".rds")
  h5ad_file <- tempfile(fileext = ".h5ad")
  roundtrip_file <- tempfile(fileext = ".rds")

  saveRDS(original, sce_file)

  scFlex::convert_sce_to_anndata(
    input = sce_file,
    output = h5ad_file
  )
  scFlex::convert_anndata_to_sce(
    input = h5ad_file,
    output = roundtrip_file
  )

  roundtrip <- readRDS(roundtrip_file)

  testthat::expect_identical(rownames(roundtrip), genes)
  testthat::expect_identical(colnames(roundtrip), cells)

  testthat::expect_equal(
    as.matrix(SummarizedExperiment::assay(roundtrip, "counts")),
    as.matrix(SummarizedExperiment::assay(original, "counts")),
    tolerance = 0
  )

  testthat::expect_equal(
    as.matrix(SummarizedExperiment::assay(roundtrip, "logcounts")),
    as.matrix(SummarizedExperiment::assay(original, "logcounts")),
    tolerance = 1e-12
  )

  roundtrip_meta <- as.data.frame(SummarizedExperiment::colData(roundtrip))
  original_meta <- as.data.frame(SummarizedExperiment::colData(original))

  testthat::expect_identical(
    as.character(roundtrip_meta$sample),
    as.character(original_meta$sample)
  )
  testthat::expect_identical(
    as.character(roundtrip_meta$condition),
    as.character(original_meta$condition)
  )

  roundtrip_feature_meta <- as.data.frame(
    SummarizedExperiment::rowData(roundtrip)
  )
  original_feature_meta <- as.data.frame(
    SummarizedExperiment::rowData(original)
  )

  testthat::expect_identical(
    as.character(roundtrip_feature_meta$symbol),
    as.character(original_feature_meta$symbol)
  )

  reduction_names <- SingleCellExperiment::reducedDimNames(roundtrip)
  testthat::expect_true(
    any(reduction_names %in% c("pca", "X_pca"))
  )

  reduction_name <- reduction_names[
    reduction_names %in% c("pca", "X_pca")
  ][1]

  testthat::expect_equal(
    unname(as.matrix(SingleCellExperiment::reducedDim(roundtrip, reduction_name))),
    unname(pca),
    tolerance = 0
  )
})


testthat::test_that("AnnData -> SCE uses integer-like X as counts when counts layer is absent", {
  testthat::skip_on_cran()
  testthat::skip_if_not_installed("reticulate")
  testthat::skip_if_not_installed("SingleCellExperiment")
  testthat::skip_if_not_installed("SummarizedExperiment")

  reticulate::py_require("anndata>=0.10")
  ad <- reticulate::import("anndata", convert = FALSE)
  np <- reticulate::import("numpy", convert = FALSE)
  pd <- reticulate::import("pandas", convert = FALSE)

  cells <- c("c1", "c2", "c3")
  genes <- c("g1", "g2", "g3", "g4")

  counts <- matrix(
    c(
      1, 0, 2, 4,
      0, 3, 1, 0,
      2, 1, 0, 5
    ),
    nrow = 3,
    byrow = TRUE
  )

  adata <- ad$AnnData(
    X = np$array(counts)$astype("float64"),
    obs = pd$DataFrame(index = cells),
    var = pd$DataFrame(index = genes)
  )

  input_file <- tempfile(fileext = ".h5ad")
  output_file <- tempfile(fileext = ".rds")
  adata$write_h5ad(input_file)

  result <- scFlex::convert_anndata_to_sce(
    input = input_file,
    output = output_file
  )

  testthat::expect_true(file.exists(output_file))
  testthat::expect_false(is.null(result))

  sce <- readRDS(output_file)

  testthat::expect_true(
    "counts" %in% SummarizedExperiment::assayNames(sce)
  )
  testthat::expect_equal(
    unname(as.matrix(SummarizedExperiment::assay(sce, "counts"))),
    unname(t(counts)),
    tolerance = 0
  )
})


testthat::test_that("AnnData -> SCE rejects non-count-like X when counts layer is absent", {
  testthat::skip_on_cran()
  testthat::skip_if_not_installed("reticulate")

  reticulate::py_require("anndata>=0.10")
  ad <- reticulate::import("anndata", convert = FALSE)
  np <- reticulate::import("numpy", convert = FALSE)
  pd <- reticulate::import("pandas", convert = FALSE)

  cells <- c("c1", "c2")
  genes <- c("g1", "g2", "g3")

  normalized <- matrix(
    c(
      0.13, 1.27, 2.45,
      0.84, 1.18, 0.51
    ),
    nrow = 2,
    byrow = TRUE
  )

  adata <- ad$AnnData(
    X = np$array(normalized),
    obs = pd$DataFrame(index = cells),
    var = pd$DataFrame(index = genes)
  )

  input_file <- tempfile(fileext = ".h5ad")
  output_file <- tempfile(fileext = ".rds")
  adata$write_h5ad(input_file)

  result <- NULL
  testthat::expect_error(
    result <- scFlex::convert_anndata_to_sce(
      input = input_file,
      output = output_file
    ),
    regexp = "no valid raw counts matrix was found"
  )

  testthat::expect_null(result)
  testthat::expect_false(file.exists(output_file))
})


testthat::test_that("AnnData -> SCE prefers explicit counts layer over count-like X", {
  testthat::skip_on_cran()
  testthat::skip_if_not_installed("reticulate")
  testthat::skip_if_not_installed("SingleCellExperiment")
  testthat::skip_if_not_installed("SummarizedExperiment")

  reticulate::py_require("anndata>=0.10")
  ad <- reticulate::import("anndata", convert = FALSE)
  np <- reticulate::import("numpy", convert = FALSE)
  pd <- reticulate::import("pandas", convert = FALSE)

  cells <- c("c1", "c2")
  genes <- c("g1", "g2", "g3")

  x_counts <- matrix(
    c(
      9, 9, 9,
      8, 8, 8
    ),
    nrow = 2,
    byrow = TRUE
  )

  explicit_counts <- matrix(
    c(
      1, 0, 2,
      0, 3, 1
    ),
    nrow = 2,
    byrow = TRUE
  )

  adata <- ad$AnnData(
    X = np$array(x_counts)$astype("float64"),
    obs = pd$DataFrame(index = cells),
    var = pd$DataFrame(index = genes)
  )
  adata$layers$`__setitem__`(
    "counts",
    np$array(explicit_counts)$astype("float64")
  )

  input_file <- tempfile(fileext = ".h5ad")
  output_file <- tempfile(fileext = ".rds")
  adata$write_h5ad(input_file)

  scFlex::convert_anndata_to_sce(
    input = input_file,
    output = output_file
  )

  sce <- readRDS(output_file)

  testthat::expect_equal(
    unname(as.matrix(SummarizedExperiment::assay(sce, "counts"))),
    unname(t(explicit_counts)),
    tolerance = 0
  )
})
