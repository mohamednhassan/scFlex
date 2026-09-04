# Functional validation for Seurat <-> Loom conversions

testthat::test_that("Seurat -> Loom preserves core matrix, names, metadata, and reductions where supported", {
  testthat::skip_if_not_installed("SeuratObject")
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

  meta <- data.frame(
    donor = c("D1", "D1", "D2", "D2"),
    score = c(1.1, 2.2, 3.3, 4.4),
    row.names = colnames(counts)
  )

  seu <- SeuratObject::CreateSeuratObject(
    counts = counts,
    meta.data = meta,
    assay = "RNA"
  )

  data_mat <- log1p(counts)
  seu <- SeuratObject::SetAssayData(
    seu,
    assay = "RNA",
    layer = "data",
    new.data = data_mat
  )

  feature_meta <- data.frame(
    symbol = paste0("SYM", 1:5),
    gene_type = c("pc", "pc", "lnc", "pc", "lnc"),
    row.names = rownames(counts)
  )
  seu[["RNA"]][[]] <- feature_meta

  pca <- matrix(
    c(
      0.1, 0.2,
      0.3, 0.4,
      0.5, 0.6,
      0.7, 0.8
    ),
    nrow = 4,
    byrow = TRUE,
    dimnames = list(colnames(counts), c("PC_1", "PC_2"))
  )

  seu[["pca"]] <- SeuratObject::CreateDimReducObject(
    embeddings = pca,
    key = "PC_",
    assay = "RNA"
  )

  input_file <- tempfile(fileext = ".rds")
  output_file <- tempfile(fileext = ".loom")
  saveRDS(seu, input_file)

  result <- scTransit::convert_seurat_to_loom(
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

  # Loom primary matrix should contain Seurat normalized data when present.
  testthat::expect_equal(
    unname(as.matrix(loom_matrix)),
    unname(as.matrix(data_mat)),
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

  # Counts should be retained as a named layer when the converter supports it.
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

  # Reductions may be stored as column attributes.
  pca_keys <- ca_keys[grepl("pca|PC", ca_keys, ignore.case = TRUE)]
  testthat::expect_true(length(pca_keys) >= 1)
})


testthat::test_that("Loom -> Seurat preserves counts-like primary matrix, names, and metadata", {
  testthat::skip_if_not_installed("SeuratObject")
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

  result <- scTransit::convert_loom_to_seurat(
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

  seu <- readRDS(output_file)

  testthat::expect_s4_class(seu, "Seurat")
  testthat::expect_equal(dim(seu), c(4, 3))
  testthat::expect_identical(colnames(seu), cells)
  testthat::expect_identical(rownames(seu), genes)

  counts_rt <- SeuratObject::GetAssayData(
    seu,
    assay = SeuratObject::DefaultAssay(seu),
    layer = "counts"
  )

  testthat::expect_equal(
    unname(as.matrix(counts_rt)),
    unname(mat),
    tolerance = 0
  )

  testthat::expect_identical(
    as.character(seu$donor),
    c("A", "A", "B")
  )

  fm <- as.data.frame(seu[[SeuratObject::DefaultAssay(seu)]][[]])

  if ("symbol" %in% colnames(fm)) {
    testthat::expect_identical(
      as.character(fm$symbol),
      paste0("SYM", 1:4)
    )
  }
})


testthat::test_that("Seurat -> Loom -> Seurat round trip preserves supported core content", {
  testthat::skip_if_not_installed("SeuratObject")
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

  seu <- SeuratObject::CreateSeuratObject(
    counts = counts,
    meta.data = data.frame(
      sample = c("S1", "S1", "S2"),
      row.names = colnames(counts)
    )
  )

  input_file <- tempfile(fileext = ".rds")
  loom_file <- tempfile(fileext = ".loom")
  output_file <- tempfile(fileext = ".rds")

  saveRDS(seu, input_file)

  scTransit::convert_seurat_to_loom(
    input = input_file,
    output = loom_file
  )
  scTransit::convert_loom_to_seurat(
    input = loom_file,
    output = output_file
  )

  rt <- readRDS(output_file)

  testthat::expect_identical(colnames(rt), colnames(seu))
  testthat::expect_identical(rownames(rt), rownames(seu))

  counts_rt <- SeuratObject::GetAssayData(
    rt,
    assay = SeuratObject::DefaultAssay(rt),
    layer = "counts"
  )

  testthat::expect_equal(
    unname(as.matrix(counts_rt)),
    unname(as.matrix(counts)),
    tolerance = 0
  )

  testthat::expect_identical(
    as.character(rt$sample),
    c("S1", "S1", "S2")
  )
})




testthat::test_that("Seurat v5 split layers are supported during Seurat -> Loom", {
  testthat::skip_if_not_installed("SeuratObject")
  testthat::skip_if_not_installed("reticulate")

  reticulate::py_require("loompy>=3.0")
  loompy <- reticulate::import("loompy", convert = FALSE)

  counts1 <- Matrix::Matrix(
    matrix(
      c(
        1, 0,
        0, 2,
        3, 1
      ),
      nrow = 3,
      byrow = TRUE
    ),
    sparse = TRUE
  )
  counts2 <- Matrix::Matrix(
    matrix(
      c(
        2,
        1,
        0
      ),
      nrow = 3
    ),
    sparse = TRUE
  )

  rownames(counts1) <- rownames(counts2) <- paste0("g", 1:3)
  colnames(counts1) <- c("c1", "c2")
  colnames(counts2) <- "c3"

  assay <- SeuratObject::CreateAssay5Object(
    counts = list(
      sample1 = counts1,
      sample2 = counts2
    )
  )

  seu <- SeuratObject::CreateSeuratObject(
    counts = assay
  )
  SeuratObject::DefaultAssay(seu) <- "RNA"

  input_file <- tempfile(fileext = ".rds")
  output_file <- tempfile(fileext = ".loom")
  saveRDS(seu, input_file)

  scTransit::convert_seurat_to_loom(
    input = input_file,
    output = output_file
  )

  testthat::expect_true(file.exists(output_file))

  ds <- loompy$connect(output_file, mode = "r")
  on.exit(ds$close(), add = TRUE)

  shape <- as.integer(reticulate::py_to_r(ds$shape))
  testthat::expect_identical(shape, c(3L, 3L))

  cell_ids <- as.character(
    reticulate::py_to_r(ds$ca$`__getitem__`("CellID"))
  )
  testthat::expect_identical(cell_ids, c("c1", "c2", "c3"))
})
