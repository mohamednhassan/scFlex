# scFlex 0.1.1

* Added automatic detection and joining of split Seurat v5 `Assay5` layers such as `counts.sample1`, `counts.sample2`, and matching `data.*` layers.
* Counts-only Seurat objects can now be written to AnnData without normalization; counts are placed in `X` and `layers["counts"]` with an explicit conversion note.
* `inspect_sc()` now reports detected count and normalized-data layer groups for `Assay5` objects.
* `check_input_info()` is retained as a compatibility wrapper but points users to `inspect_sc()`.
* Improved assay-cell, feature-metadata, and reduction alignment during Seurat extraction.

# scFlex 0.1.0

* Initial development release.
* Adds a common `convert_sc()` interface for Seurat, SingleCellExperiment, AnnData, and Loom.
* Preserves aligned counts, normalized expression, cell/feature metadata, and compatible dimensional reductions.
* Adds Seurat classic Assay <-> Assay5 conversion helpers.
* Uses `reticulate::py_require()` instead of a hard-coded Conda environment.
