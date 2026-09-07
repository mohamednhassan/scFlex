# Functional validation for Seurat <-> AnnData conversions

testthat::test_that("Seurat -> AnnData preserves core data components", {
  testthat::skip_if_not_installed("reticulate")
  testthat::skip_if_not_installed("SeuratObject")
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
  data_mat <- Matrix::Matrix(log1p(counts_dense), sparse = TRUE)

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

  seu <- SeuratObject::CreateSeuratObject(
    counts = counts,
    meta.data = metadata
  )

  assay_name <- SeuratObject::DefaultAssay(seu)

  seu <- SeuratObject::SetAssayData(
    seu,
    assay = assay_name,
    layer = "data",
    new.data = data_mat
  )

  seu[[assay_name]][[]] <- feature_metadata

  seu[["pca"]] <- SeuratObject::CreateDimReducObject(
    embeddings = pca,
    key = "PC_",
    assay = assay_name
  )

  input_file <- tempfile(fileext = ".rds")
  output_file <- tempfile(fileext = ".h5ad")
  saveRDS(seu, input_file)

  result <- scFlex::convert_seurat_to_anndata(
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

  x <- reticulate::py_to_r(adata$X)
  counts_layer <- reticulate::py_to_r(
    adata$layers$`__getitem__`("counts")
  )

  testthat::expect_equal(
    unname(as.matrix(x)),
    unname(t(as.matrix(data_mat))),
    tolerance = 1e-12
  )

  testthat::expect_equal(
    unname(as.matrix(counts_layer)),
    unname(t(as.matrix(counts))),
    tolerance = 0
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

  obsm_keys <- unlist(
    reticulate::iterate(adata$obsm$keys(), as.character),
    use.names = FALSE
  )

  testthat::expect_true(any(obsm_keys %in% c("X_pca", "pca")))

  pca_key <- obsm_keys[obsm_keys %in% c("X_pca", "pca")][1]
  pca_out <- reticulate::py_to_r(
    adata$obsm$`__getitem__`(pca_key)
  )

  testthat::expect_equal(
    unname(as.matrix(pca_out)),
    unname(pca),
    tolerance = 0
  )
})


testthat::test_that("AnnData -> Seurat preserves core data components", {
  testthat::skip_if_not_installed("reticulate")
  testthat::skip_if_not_installed("SeuratObject")

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

  normalized <- log1p(counts) + 0.125

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
    X = np$array(normalized),
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
  output_file <- tempfile(fileext = ".rds")
  adata$write_h5ad(input_file)

  result <- scFlex::convert_anndata_to_seurat(
    input = input_file,
    output = output_file
  )

  testthat::expect_true(file.exists(output_file))
  testthat::expect_identical(
    normalizePath(result, mustWork = TRUE),
    normalizePath(output_file, mustWork = TRUE)
  )

  seu <- readRDS(output_file)

  testthat::expect_s4_class(seu, "Seurat")
  testthat::expect_identical(rownames(seu), genes)
  testthat::expect_identical(colnames(seu), cells)

  assay_name <- SeuratObject::DefaultAssay(seu)

  seu_counts <- SeuratObject::LayerData(
    seu[[assay_name]],
    layer = "counts"
  )

  testthat::expect_equal(
    unname(as.matrix(seu_counts)),
    unname(t(counts)),
    tolerance = 0
  )

  layers <- SeuratObject::Layers(seu[[assay_name]])
  testthat::expect_true("data" %in% layers)

  seu_data <- SeuratObject::LayerData(
    seu[[assay_name]],
    layer = "data"
  )

  testthat::expect_equal(
    unname(as.matrix(seu_data)),
    unname(t(normalized)),
    tolerance = 1e-12
  )

  testthat::expect_identical(rownames(seu@meta.data), cells)
  testthat::expect_identical(
    as.character(seu$donor),
    c("A", "A", "B", "B")
  )
  testthat::expect_equal(
    unname(as.numeric(seu$score)),
    c(1.1, 2.2, 3.3, 4.4)
  )

  feature_meta <- as.data.frame(seu[[assay_name]][[]])
  testthat::expect_identical(rownames(feature_meta), genes)
  testthat::expect_identical(
    as.character(feature_meta$gene_type),
    c("pc", "pc", "lnc", "pc", "lnc")
  )
  testthat::expect_equal(
    as.numeric(feature_meta$gc),
    c(0.41, 0.52, 0.63, 0.44, 0.55)
  )

  reduction_names <- SeuratObject::Reductions(seu)
  testthat::expect_true(any(reduction_names %in% c("pca", "X_pca")))

  reduction_name <- reduction_names[
    reduction_names %in% c("pca", "X_pca")
  ][1]

  testthat::expect_equal(
    unname(SeuratObject::Embeddings(seu[[reduction_name]])),
    unname(pca),
    tolerance = 0
  )
})


testthat::test_that("Seurat -> AnnData -> Seurat round trip preserves supported content", {
  testthat::skip_if_not_installed("reticulate")
  testthat::skip_if_not_installed("SeuratObject")
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
  data_mat <- Matrix::Matrix(log1p(counts_dense), sparse = TRUE)

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

  original <- SeuratObject::CreateSeuratObject(
    counts = counts,
    meta.data = metadata
  )

  assay_name <- SeuratObject::DefaultAssay(original)

  original <- SeuratObject::SetAssayData(
    original,
    assay = assay_name,
    layer = "data",
    new.data = data_mat
  )

  original[[assay_name]][[]] <- feature_metadata

  original[["pca"]] <- SeuratObject::CreateDimReducObject(
    embeddings = pca,
    key = "PC_",
    assay = assay_name
  )

  seurat_file <- tempfile(fileext = ".rds")
  h5ad_file <- tempfile(fileext = ".h5ad")
  roundtrip_file <- tempfile(fileext = ".rds")

  saveRDS(original, seurat_file)

  scFlex::convert_seurat_to_anndata(
    input = seurat_file,
    output = h5ad_file
  )

  scFlex::convert_anndata_to_seurat(
    input = h5ad_file,
    output = roundtrip_file
  )

  roundtrip <- readRDS(roundtrip_file)
  roundtrip_assay <- SeuratObject::DefaultAssay(roundtrip)

  testthat::expect_identical(rownames(roundtrip), genes)
  testthat::expect_identical(colnames(roundtrip), cells)

  testthat::expect_equal(
    unname(as.matrix(
      SeuratObject::LayerData(
        roundtrip[[roundtrip_assay]],
        layer = "counts"
      )
    )),
    unname(as.matrix(
      SeuratObject::LayerData(
        original[[assay_name]],
        layer = "counts"
      )
    )),
    tolerance = 0
  )

  testthat::expect_equal(
    unname(as.matrix(
      SeuratObject::LayerData(
        roundtrip[[roundtrip_assay]],
        layer = "data"
      )
    )),
    unname(as.matrix(
      SeuratObject::LayerData(
        original[[assay_name]],
        layer = "data"
      )
    )),
    tolerance = 1e-12
  )

  testthat::expect_identical(
    as.character(roundtrip$sample),
    as.character(original$sample)
  )

  testthat::expect_identical(
    as.character(roundtrip$condition),
    as.character(original$condition)
  )

  original_feature_meta <- as.data.frame(original[[assay_name]][[]])
  roundtrip_feature_meta <- as.data.frame(roundtrip[[roundtrip_assay]][[]])

  testthat::expect_identical(
    as.character(roundtrip_feature_meta$symbol),
    as.character(original_feature_meta$symbol)
  )

  reduction_names <- SeuratObject::Reductions(roundtrip)
  testthat::expect_true(any(reduction_names %in% c("pca", "X_pca")))

  reduction_name <- reduction_names[
    reduction_names %in% c("pca", "X_pca")
  ][1]

  testthat::expect_equal(
    unname(SeuratObject::Embeddings(roundtrip[[reduction_name]])),
    unname(pca),
    tolerance = 0
  )
})


testthat::test_that("counts-only Seurat -> AnnData keeps counts without inventing normalized values", {
  testthat::skip_if_not_installed("reticulate")
  testthat::skip_if_not_installed("SeuratObject")
  testthat::skip_if_not_installed("Matrix")

  reticulate::py_require("anndata>=0.10")
  ad <- reticulate::import("anndata", convert = FALSE)

  counts <- Matrix::Matrix(
    matrix(
      c(
        1, 0,
        0, 3,
        2, 1
      ),
      nrow = 3,
      byrow = TRUE,
      dimnames = list(
        paste0("gene", 1:3),
        paste0("cell", 1:2)
      )
    ),
    sparse = TRUE
  )

  seu <- SeuratObject::CreateSeuratObject(counts = counts)

  input_file <- tempfile(fileext = ".rds")
  output_file <- tempfile(fileext = ".h5ad")
  saveRDS(seu, input_file)

  scFlex::convert_seurat_to_anndata(
    input = input_file,
    output = output_file
  )

  adata <- ad$read_h5ad(output_file)

  layer_keys <- unlist(
    reticulate::iterate(adata$layers$keys(), as.character),
    use.names = FALSE
  )

  testthat::expect_true("counts" %in% layer_keys)

  counts_layer <- reticulate::py_to_r(
    adata$layers$`__getitem__`("counts")
  )

  testthat::expect_equal(
    unname(as.matrix(counts_layer)),
    unname(t(as.matrix(counts))),
    tolerance = 0
  )

  x <- reticulate::py_to_r(adata$X)

  testthat::expect_equal(
    unname(as.matrix(x)),
    unname(t(as.matrix(counts))),
    tolerance = 0
  )
})


testthat::test_that("AnnData -> Seurat uses integer-like X as counts when counts layer is absent", {
  testthat::skip_if_not_installed("reticulate")
  testthat::skip_if_not_installed("SeuratObject")

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

  result <- scFlex::convert_anndata_to_seurat(
    input = input_file,
    output = output_file
  )

  testthat::expect_true(file.exists(output_file))
  testthat::expect_false(is.null(result))

  seu <- readRDS(output_file)
  assay_name <- SeuratObject::DefaultAssay(seu)

  testthat::expect_equal(
    unname(as.matrix(
      SeuratObject::LayerData(
        seu[[assay_name]],
        layer = "counts"
      )
    )),
    unname(t(counts)),
    tolerance = 0
  )
})


testthat::test_that("AnnData -> Seurat rejects non-count-like X when counts layer is absent", {
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
    result <- scFlex::convert_anndata_to_seurat(
      input = input_file,
      output = output_file
    ),
    regexp = "no valid raw counts matrix was found"
  )

  testthat::expect_null(result)
  testthat::expect_false(file.exists(output_file))
})


testthat::test_that("AnnData -> Seurat prefers explicit counts layer over count-like X", {
  testthat::skip_if_not_installed("reticulate")
  testthat::skip_if_not_installed("SeuratObject")

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

  scFlex::convert_anndata_to_seurat(
    input = input_file,
    output = output_file
  )

  seu <- readRDS(output_file)
  assay_name <- SeuratObject::DefaultAssay(seu)

  testthat::expect_equal(
    unname(as.matrix(
      SeuratObject::LayerData(
        seu[[assay_name]],
        layer = "counts"
      )
    )),
    unname(t(explicit_counts)),
    tolerance = 0
  )
})


testthat::test_that("split Seurat v5 layers are joined during Seurat -> AnnData", {
  testthat::skip_if_not_installed("reticulate")
  testthat::skip_if_not_installed("SeuratObject")
  testthat::skip_if_not_installed("Matrix")

  testthat::skip_if(
    utils::packageVersion("SeuratObject") < "5.0.0",
    "Split Assay5 layers require SeuratObject >= 5.0.0"
  )

  reticulate::py_require("anndata>=0.10")
  ad <- reticulate::import("anndata", convert = FALSE)

  genes <- paste0("gene", 1:4)
  cells <- paste0("cell", 1:6)

  counts_dense <- matrix(
    c(
      1, 0, 2, 0, 3, 1,
      0, 2, 0, 4, 0, 1,
      3, 1, 0, 0, 2, 0,
      0, 0, 1, 2, 1, 5
    ),
    nrow = length(genes),
    byrow = TRUE,
    dimnames = list(genes, cells)
  )

  counts <- Matrix::Matrix(counts_dense, sparse = TRUE)
  data_mat <- Matrix::Matrix(log1p(counts_dense), sparse = TRUE)

  seu <- SeuratObject::CreateSeuratObject(counts = counts)
  assay_name <- SeuratObject::DefaultAssay(seu)

  seu <- SeuratObject::SetAssayData(
    seu,
    assay = assay_name,
    layer = "data",
    new.data = data_mat
  )

  groups <- factor(
    c("sampleA", "sampleA", "sampleA", "sampleB", "sampleB", "sampleB")
  )

  seu[[assay_name]] <- split(
    seu[[assay_name]],
    f = groups
  )

  split_layers <- SeuratObject::Layers(seu[[assay_name]])

  testthat::expect_true(any(grepl("^counts\\.", split_layers)))
  testthat::expect_true(any(grepl("^data\\.", split_layers)))

  input_file <- tempfile(fileext = ".rds")
  output_file <- tempfile(fileext = ".h5ad")
  saveRDS(seu, input_file)

  scFlex::convert_seurat_to_anndata(
    input = input_file,
    output = output_file
  )

  adata <- ad$read_h5ad(output_file)

  counts_out <- reticulate::py_to_r(
    adata$layers$`__getitem__`("counts")
  )
  x_out <- reticulate::py_to_r(adata$X)

  testthat::expect_equal(
    unname(as.matrix(counts_out)),
    unname(t(counts_dense)),
    tolerance = 0
  )

  testthat::expect_equal(
    unname(as.matrix(x_out)),
    unname(t(as.matrix(data_mat))),
    tolerance = 1e-12
  )
})
