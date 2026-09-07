# Regression tests for issues discovered during validation on real datasets

testthat::test_that("AnnData explicit counts layer accepts fractional non-negative counts", {
  testthat::skip_if_not_installed("reticulate")
  testthat::skip_if_not_installed("SeuratObject")
  testthat::skip_if_not_installed("SingleCellExperiment")
  testthat::skip_if_not_installed("SummarizedExperiment")

  reticulate::py_require("anndata>=0.10")

  ad <- reticulate::import("anndata", convert = FALSE)
  np <- reticulate::import("numpy", convert = FALSE)

  cells <- c("cell1", "cell2", "cell3")
  genes <- c("gene1", "gene2", "gene3", "gene4")

  counts <- matrix(
    c(
      0.0, 1.25, 0.0, 2.50,
      0.5, 0.00, 3.75, 0.00,
      1.1, 2.20, 0.0, 4.40
    ),
    nrow = length(cells),
    byrow = TRUE
  )

  normalized <- log1p(counts)

  adata <- ad$AnnData(X = np$array(normalized))
  adata$obs_names <- cells
  adata$var_names <- genes
  adata$layers$`__setitem__`("counts", np$array(counts))

  input_file <- tempfile(fileext = ".h5ad")
  seurat_file <- tempfile(fileext = ".rds")
  sce_file <- tempfile(fileext = ".rds")

  adata$write_h5ad(input_file)

  scTransit::convert_anndata_to_seurat(
    input = input_file,
    output = seurat_file
  )

  scTransit::convert_anndata_to_sce(
    input = input_file,
    output = sce_file
  )

  seu <- readRDS(seurat_file)
  sce <- readRDS(sce_file)

  seu_counts <- SeuratObject::LayerData(
    seu,
    assay = SeuratObject::DefaultAssay(seu),
    layer = "counts"
  )

  sce_counts <- SummarizedExperiment::assay(sce, "counts")

  testthat::expect_equal(
    unname(as.matrix(seu_counts)),
    unname(t(counts)),
    tolerance = 0
  )

  testthat::expect_equal(
    unname(as.matrix(sce_counts)),
    unname(t(counts)),
    tolerance = 0
  )
})


testthat::test_that("Loom explicit counts layer accepts fractional non-negative counts", {
  testthat::skip_if_not_installed("reticulate")
  testthat::skip_if_not_installed("SeuratObject")
  testthat::skip_if_not_installed("SingleCellExperiment")
  testthat::skip_if_not_installed("SummarizedExperiment")

  reticulate::py_require("anndata>=0.10")
  reticulate::py_require("loompy>=3.0")

  ad <- reticulate::import("anndata", convert = FALSE)
  np <- reticulate::import("numpy", convert = FALSE)

  cells <- c("cell1", "cell2", "cell3")
  genes <- c("gene1", "gene2", "gene3", "gene4")

  counts <- matrix(
    c(
      0.0, 1.25, 0.0, 2.50,
      0.5, 0.00, 3.75, 0.00,
      1.1, 2.20, 0.0, 4.40
    ),
    nrow = length(cells),
    byrow = TRUE
  )

  normalized <- log1p(counts)

  adata <- ad$AnnData(X = np$array(normalized))
  adata$obs_names <- cells
  adata$var_names <- genes
  adata$layers$`__setitem__`("counts", np$array(counts))
  adata$layers$`__setitem__`("logcounts", np$array(normalized))

  h5ad_file <- tempfile(fileext = ".h5ad")
  loom_file <- tempfile(fileext = ".loom")
  seurat_file <- tempfile(fileext = ".rds")
  sce_file <- tempfile(fileext = ".rds")

  adata$write_h5ad(h5ad_file)

  scTransit::convert_anndata_to_loom(
    input = h5ad_file,
    output = loom_file
  )

  suppressWarnings(
    scTransit::convert_loom_to_seurat(
      input = loom_file,
      output = seurat_file
    )
  )

  suppressWarnings(
    scTransit::convert_loom_to_sce(
      input = loom_file,
      output = sce_file
    )
  )

  seu <- readRDS(seurat_file)
  sce <- readRDS(sce_file)

  seu_counts <- SeuratObject::LayerData(
    seu,
    assay = SeuratObject::DefaultAssay(seu),
    layer = "counts"
  )

  sce_counts <- SummarizedExperiment::assay(sce, "counts")

  testthat::expect_equal(
    unname(as.matrix(seu_counts)),
    unname(t(counts)),
    tolerance = 0
  )

  testthat::expect_equal(
    unname(as.matrix(sce_counts)),
    unname(t(counts)),
    tolerance = 0
  )
})


testthat::test_that("inspect_sc reports Loom cells and features in correct orientation", {
  testthat::skip_if_not_installed("reticulate")
  testthat::skip_if_not_installed("hdf5r")

  reticulate::py_require("anndata>=0.10")
  reticulate::py_require("loompy>=3.0")

  ad <- reticulate::import("anndata", convert = FALSE)
  np <- reticulate::import("numpy", convert = FALSE)

  cells <- paste0("cell", 1:3)
  genes <- paste0("gene", 1:5)

  x <- matrix(
    seq_len(length(cells) * length(genes)),
    nrow = length(cells),
    ncol = length(genes)
  )

  adata <- ad$AnnData(X = np$array(x))
  adata$obs_names <- cells
  adata$var_names <- genes

  h5ad_file <- tempfile(fileext = ".h5ad")
  loom_file <- tempfile(fileext = ".loom")

  adata$write_h5ad(h5ad_file)

  scTransit::convert_anndata_to_loom(
    input = h5ad_file,
    output = loom_file
  )

  info <- suppressMessages(
    scTransit::inspect_sc(loom_file)
  )

  testthat::expect_identical(
    as.integer(info$num_cells),
    length(cells)
  )

  testthat::expect_identical(
    as.integer(info$num_features),
    length(genes)
  )
})


testthat::test_that("Loom round trip does not add original_feature_name when feature names are unique", {
  testthat::skip_if_not_installed("reticulate")
  testthat::skip_if_not_installed("SeuratObject")

  reticulate::py_require("anndata>=0.10")
  reticulate::py_require("loompy>=3.0")

  ad <- reticulate::import("anndata", convert = FALSE)
  np <- reticulate::import("numpy", convert = FALSE)

  cells <- paste0("cell", 1:3)
  genes <- paste0("gene", 1:4)

  counts <- matrix(
    c(
      1, 0, 2, 0,
      0, 3, 0, 1,
      2, 1, 0, 4
    ),
    nrow = length(cells),
    byrow = TRUE
  )

  adata <- ad$AnnData(X = np$array(log1p(counts)))
  adata$obs_names <- cells
  adata$var_names <- genes
  adata$layers$`__setitem__`("counts", np$array(counts))

  h5ad_file <- tempfile(fileext = ".h5ad")
  loom_file <- tempfile(fileext = ".loom")
  seurat_file <- tempfile(fileext = ".rds")

  adata$write_h5ad(h5ad_file)

  scTransit::convert_anndata_to_loom(
    input = h5ad_file,
    output = loom_file
  )

  suppressWarnings(
    scTransit::convert_loom_to_seurat(
      input = loom_file,
      output = seurat_file
    )
  )

  seu <- readRDS(seurat_file)

  feature_metadata <- seu[[SeuratObject::DefaultAssay(seu)]][[]]

  testthat::expect_false(
    "original_feature_name" %in% colnames(feature_metadata)
  )
})


testthat::test_that("Loom preserves original feature names when duplicates require uniquification", {
  testthat::skip_if_not_installed("reticulate")
  testthat::skip_if_not_installed("SeuratObject")

  reticulate::py_require("anndata>=0.10")
  reticulate::py_require("loompy>=3.0")

  ad <- reticulate::import("anndata", convert = FALSE)
  np <- reticulate::import("numpy", convert = FALSE)

  cells <- paste0("cell", 1:2)
  original_genes <- c("geneA", "geneA", "geneB")

  counts <- matrix(
    c(
      1, 2, 0,
      0, 3, 4
    ),
    nrow = length(cells),
    byrow = TRUE
  )

  adata <- ad$AnnData(X = np$array(log1p(counts)))
  adata$obs_names <- cells
  adata$var_names <- original_genes
  adata$layers$`__setitem__`("counts", np$array(counts))

  h5ad_file <- tempfile(fileext = ".h5ad")
  loom_file <- tempfile(fileext = ".loom")
  seurat_file <- tempfile(fileext = ".rds")

  adata$write_h5ad(h5ad_file)

  suppressWarnings(
    scTransit::convert_anndata_to_loom(
      input = h5ad_file,
      output = loom_file
    )
  )

  suppressWarnings(
    scTransit::convert_loom_to_seurat(
      input = loom_file,
      output = seurat_file
    )
  )

  seu <- readRDS(seurat_file)

  feature_metadata <- seu[[SeuratObject::DefaultAssay(seu)]][[]]

  testthat::expect_true(
    "original_feature_name" %in% colnames(feature_metadata)
  )

  testthat::expect_identical(
    as.character(feature_metadata$original_feature_name),
    original_genes
  )

  testthat::expect_equal(
    anyDuplicated(rownames(feature_metadata)),
    0
  )
})



test_that("AnnData conversions error when no valid counts source exists", {

  skip_if_not_installed("reticulate")
  skip_if_not_installed("Seurat")
  skip_if_not_installed("SingleCellExperiment")

  np <- reticulate::import("numpy", convert = FALSE)
  anndata <- reticulate::import("anndata", convert = FALSE)

  # Deliberately processed-looking matrix:
  # contains negative and non-integer values, so X must not be
  # interpreted as raw counts.
  x <- matrix(
    c(
      -1.2, 0.5, 2.3,
       0.1, -0.4, 1.7
    ),
    nrow = 2,
    byrow = TRUE
  )

  adata <- anndata$AnnData(
    X = np$array(x)
  )

  adata$obs_names <- c("cell1", "cell2")
  adata$var_names <- c("gene1", "gene2", "gene3")

  input <- tempfile(fileext = ".h5ad")
  seurat_output <- tempfile(fileext = ".RDS")
  sce_output <- tempfile(fileext = ".RDS")

  adata$write_h5ad(
    input,
    convert_strings_to_categoricals = FALSE
  )

  expect_error(
    convert_anndata_to_seurat(
      input = input,
      output = seurat_output
    ),
    "no valid raw counts matrix was found"
  )

  expect_false(
    file.exists(seurat_output)
  )

  expect_error(
    convert_anndata_to_sce(
      input = input,
      output = sce_output
    ),
    "no valid raw counts matrix was found"
  )

  expect_false(
    file.exists(sce_output)
  )
})
