# Functional validation for Seurat <-> SingleCellExperiment conversions

testthat::test_that("Seurat -> SCE preserves core data components", {
  testthat::skip_if_not_installed("SeuratObject")
  testthat::skip_if_not_installed("SingleCellExperiment")
  testthat::skip_if_not_installed("SummarizedExperiment")
  testthat::skip_if_not_installed("S4Vectors")
  testthat::skip_if_not_installed("Matrix")

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

  embeddings <- matrix(
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
    embeddings = embeddings,
    key = "PC_",
    assay = assay_name
  )

  input_file <- tempfile(fileext = ".rds")
  output_file <- tempfile(fileext = ".rds")
  saveRDS(seu, input_file)

  result <- scFlex::convert_seurat_to_sce(
    input = input_file,
    output = output_file
  )

  testthat::expect_identical(
    normalizePath(result, mustWork = TRUE),
    normalizePath(output_file, mustWork = TRUE)
  )
  testthat::expect_true(file.exists(output_file))

  sce <- readRDS(output_file)

  testthat::expect_s4_class(sce, "SingleCellExperiment")
  testthat::expect_identical(rownames(sce), genes)
  testthat::expect_identical(colnames(sce), cells)
  testthat::expect_identical(dim(sce), c(length(genes), length(cells)))

  sce_counts <- SummarizedExperiment::assay(sce, "counts")
  sce_logcounts <- SummarizedExperiment::assay(sce, "logcounts")

  testthat::expect_equal(
    as.matrix(sce_counts),
    as.matrix(counts),
    tolerance = 0
  )
  testthat::expect_equal(
    as.matrix(sce_logcounts),
    as.matrix(data_mat),
    tolerance = 1e-12
  )

  sce_metadata <- as.data.frame(SummarizedExperiment::colData(sce))
  testthat::expect_identical(rownames(sce_metadata), cells)
  testthat::expect_identical(as.character(sce_metadata$donor), metadata$donor)
  testthat::expect_identical(
    as.character(sce_metadata$condition),
    as.character(metadata$condition)
  )
  testthat::expect_equal(sce_metadata$score, metadata$score)

  sce_feature_metadata <- as.data.frame(SummarizedExperiment::rowData(sce))
  testthat::expect_identical(rownames(sce_feature_metadata), genes)
  testthat::expect_identical(
    as.character(sce_feature_metadata$gene_type),
    feature_metadata$gene_type
  )
  testthat::expect_equal(sce_feature_metadata$gc, feature_metadata$gc)

  testthat::expect_true(
    "pca" %in% SingleCellExperiment::reducedDimNames(sce)
  )
  testthat::expect_equal(
    as.matrix(SingleCellExperiment::reducedDim(sce, "pca")),
    embeddings,
    tolerance = 0
  )
})


testthat::test_that("SCE -> Seurat preserves core data components", {
  testthat::skip_if_not_installed("SeuratObject")
  testthat::skip_if_not_installed("SingleCellExperiment")
  testthat::skip_if_not_installed("SummarizedExperiment")
  testthat::skip_if_not_installed("S4Vectors")
  testthat::skip_if_not_installed("Matrix")

  genes <- paste0("gene", 1:6)
  cells <- paste0("cell", 1:5)

  counts_dense <- matrix(
    c(
      2, 0, 1, 4, 0,
      0, 3, 0, 1, 2,
      5, 1, 0, 0, 2,
      0, 2, 3, 0, 1,
      1, 0, 4, 2, 0,
      3, 1, 0, 1, 5
    ),
    nrow = length(genes),
    byrow = TRUE,
    dimnames = list(genes, cells)
  )

  counts <- Matrix::Matrix(counts_dense, sparse = TRUE)
  logcounts <- Matrix::Matrix(log1p(counts_dense), sparse = TRUE)

  metadata <- data.frame(
    donor = c("A", "A", "B", "B", "C"),
    group = factor(c("X", "X", "Y", "Y", "X")),
    age = c(21, 22, 31, 32, 41),
    row.names = cells
  )

  feature_metadata <- data.frame(
    chromosome = c("1", "1", "2", "2", "3", "X"),
    marker = c(TRUE, FALSE, TRUE, FALSE, TRUE, FALSE),
    row.names = genes
  )

  embeddings <- matrix(
    c(
      -0.5, 0.1,
      -0.2, 0.2,
       0.0, 0.3,
       0.2, 0.4,
       0.5, 0.5
    ),
    nrow = length(cells),
    byrow = TRUE,
    dimnames = list(cells, c("UMAP_1", "UMAP_2"))
  )

  sce <- SingleCellExperiment::SingleCellExperiment(
    assays = list(
      counts = counts,
      logcounts = logcounts
    ),
    colData = S4Vectors::DataFrame(metadata),
    rowData = S4Vectors::DataFrame(feature_metadata)
  )

  SingleCellExperiment::reducedDim(sce, "umap") <- embeddings

  input_file <- tempfile(fileext = ".rds")
  output_file <- tempfile(fileext = ".rds")
  saveRDS(sce, input_file)

  result <- scFlex::convert_sce_to_seurat(
    input = input_file,
    output = output_file
  )

  testthat::expect_identical(
    normalizePath(result, mustWork = TRUE),
    normalizePath(output_file, mustWork = TRUE)
  )
  testthat::expect_true(file.exists(output_file))

  seu <- readRDS(output_file)

  testthat::expect_s4_class(seu, "Seurat")
  testthat::expect_identical(rownames(seu), genes)
  testthat::expect_identical(colnames(seu), cells)
  testthat::expect_equal(dim(seu), c(length(genes), length(cells)))

  assay_name <- SeuratObject::DefaultAssay(seu)

  seu_counts <- SeuratObject::LayerData(
    seu[[assay_name]],
    layer = "counts"
  )
  seu_data <- SeuratObject::LayerData(
    seu[[assay_name]],
    layer = "data"
  )

  testthat::expect_equal(
    as.matrix(seu_counts),
    as.matrix(counts),
    tolerance = 0
  )
  testthat::expect_equal(
    as.matrix(seu_data),
    as.matrix(logcounts),
    tolerance = 1e-12
  )

  testthat::expect_identical(rownames(seu@meta.data), cells)
  testthat::expect_identical(as.character(seu$donor), metadata$donor)
  testthat::expect_identical(
    as.character(seu$group),
    as.character(metadata$group)
  )
  testthat::expect_equal(unname(seu$age), metadata$age)

  seu_feature_metadata <- as.data.frame(seu[[assay_name]][[]])
  testthat::expect_identical(rownames(seu_feature_metadata), genes)
  testthat::expect_identical(
    as.character(seu_feature_metadata$chromosome),
    feature_metadata$chromosome
  )
  testthat::expect_identical(
    as.logical(seu_feature_metadata$marker),
    feature_metadata$marker
  )

  testthat::expect_true("umap" %in% SeuratObject::Reductions(seu))
  testthat::expect_equal(
    SeuratObject::Embeddings(seu[["umap"]]),
    embeddings,
    tolerance = 0
  )
})


testthat::test_that("Seurat -> SCE -> Seurat round trip preserves supported content", {
  testthat::skip_if_not_installed("SeuratObject")
  testthat::skip_if_not_installed("SingleCellExperiment")
  testthat::skip_if_not_installed("SummarizedExperiment")
  testthat::skip_if_not_installed("S4Vectors")
  testthat::skip_if_not_installed("Matrix")

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

  embeddings <- matrix(
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
    embeddings = embeddings,
    key = "PC_",
    assay = assay_name
  )

  seurat_file <- tempfile(fileext = ".rds")
  sce_file <- tempfile(fileext = ".rds")
  roundtrip_file <- tempfile(fileext = ".rds")

  saveRDS(original, seurat_file)

  scFlex::convert_seurat_to_sce(
    input = seurat_file,
    output = sce_file
  )
  scFlex::convert_sce_to_seurat(
    input = sce_file,
    output = roundtrip_file
  )

  roundtrip <- readRDS(roundtrip_file)
  roundtrip_assay <- SeuratObject::DefaultAssay(roundtrip)

  testthat::expect_identical(rownames(roundtrip), genes)
  testthat::expect_identical(colnames(roundtrip), cells)

  testthat::expect_equal(
    as.matrix(
      SeuratObject::LayerData(original[[assay_name]], layer = "counts")
    ),
    as.matrix(
      SeuratObject::LayerData(roundtrip[[roundtrip_assay]], layer = "counts")
    ),
    tolerance = 0
  )

  testthat::expect_equal(
    as.matrix(
      SeuratObject::LayerData(original[[assay_name]], layer = "data")
    ),
    as.matrix(
      SeuratObject::LayerData(roundtrip[[roundtrip_assay]], layer = "data")
    ),
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

  original_feature_metadata <- as.data.frame(original[[assay_name]][[]])
  roundtrip_feature_metadata <- as.data.frame(roundtrip[[roundtrip_assay]][[]])

  testthat::expect_identical(
    as.character(roundtrip_feature_metadata$symbol),
    as.character(original_feature_metadata$symbol)
  )

  testthat::expect_true("pca" %in% SeuratObject::Reductions(roundtrip))
  testthat::expect_equal(
    SeuratObject::Embeddings(roundtrip[["pca"]]),
    SeuratObject::Embeddings(original[["pca"]]),
    tolerance = 0
  )
})


testthat::test_that("counts-only SCE -> Seurat does not invent normalized data", {
  testthat::skip_if_not_installed("SeuratObject")
  testthat::skip_if_not_installed("SingleCellExperiment")
  testthat::skip_if_not_installed("S4Vectors")
  testthat::skip_if_not_installed("Matrix")

  counts <- Matrix::Matrix(
    matrix(
      c(1, 0, 2, 3, 1, 0),
      nrow = 3,
      dimnames = list(
        paste0("gene", 1:3),
        paste0("cell", 1:2)
      )
    ),
    sparse = TRUE
  )

  sce <- SingleCellExperiment::SingleCellExperiment(
    assays = list(counts = counts)
  )

  input_file <- tempfile(fileext = ".rds")
  output_file <- tempfile(fileext = ".rds")
  saveRDS(sce, input_file)

  scFlex::convert_sce_to_seurat(
    input = input_file,
    output = output_file
  )

  seu <- readRDS(output_file)
  assay_name <- SeuratObject::DefaultAssay(seu)
  layers <- SeuratObject::Layers(seu[[assay_name]])

  testthat::expect_true("counts" %in% layers)

  if ("data" %in% layers) {
    data_layer <- SeuratObject::LayerData(
      seu[[assay_name]],
      layer = "data"
    )
    testthat::expect_true(
      nrow(data_layer) == 0L ||
        ncol(data_layer) == 0L ||
        length(data_layer@x) == 0L
    )
  }
})


testthat::test_that("split Seurat v5 count layers are joined during Seurat -> SCE", {
  testthat::skip_if_not_installed("SeuratObject")
  testthat::skip_if_not_installed("SingleCellExperiment")
  testthat::skip_if_not_installed("SummarizedExperiment")
  testthat::skip_if_not_installed("Matrix")

  testthat::skip_if(
    utils::packageVersion("SeuratObject") < "5.0.0",
    "Split Assay5 layers require SeuratObject >= 5.0.0"
  )

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

  seu <- SeuratObject::CreateSeuratObject(counts = counts)
  assay_name <- SeuratObject::DefaultAssay(seu)

  groups <- factor(
    c("sampleA", "sampleA", "sampleA", "sampleB", "sampleB", "sampleB")
  )

  seu[[assay_name]] <- split(
    seu[[assay_name]],
    f = groups
  )

  split_layers <- SeuratObject::Layers(seu[[assay_name]])
  testthat::expect_true(any(grepl("^counts\\.", split_layers)))

  input_file <- tempfile(fileext = ".rds")
  output_file <- tempfile(fileext = ".rds")
  saveRDS(seu, input_file)

  scFlex::convert_seurat_to_sce(
    input = input_file,
    output = output_file
  )

  sce <- readRDS(output_file)
  converted_counts <- SummarizedExperiment::assay(sce, "counts")

  testthat::expect_identical(rownames(converted_counts), genes)
  testthat::expect_identical(colnames(converted_counts), cells)
  testthat::expect_equal(
    as.matrix(converted_counts),
    counts_dense,
    tolerance = 0
  )
})
