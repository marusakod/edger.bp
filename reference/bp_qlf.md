# End-to-End BPCells edgeR QLF Pipeline

Runs the supported BPCells edgeR quasi-likelihood workflow: preparation,
streamed or dense-reference dispersion estimation when needed, and
blockwise quasi-likelihood F-testing.

## Usage

``` r
bp_qlf(
  y,
  design,
  group = NULL,
  coef = ncol(design),
  contrast = NULL,
  filter = TRUE,
  dispersion = NULL,
  norm.method = "TMM",
  dispersion.method = c("stream", "dense"),
  block.size = 1000,
  workers = 1,
  max.elements = 1e+07,
  prior.count = 0.125,
  filter.args = list(),
  dispersion.args = list(),
  qlf.args = list()
)
```

## Arguments

- y:

  A genes-by-samples count matrix.

- design:

  Design matrix.

- group:

  Optional group factor for filtering and dispersion estimation.

- coef, contrast:

  Test specification passed to
  [`bp_glm_qlf()`](https://marusakod.github.io/edger.bp/reference/bp_glm_qlf.md).

- filter:

  Logical scalar. If `TRUE`, apply
  [`bp_filter_by_expr()`](https://marusakod.github.io/edger.bp/reference/bp_filter_by_expr.md).

- dispersion:

  Optional scalar or gene-level dispersion vector. If `NULL`,
  dispersions are estimated with `dispersion.method`.

- norm.method:

  Normalization method passed to
  [`bp_prepare_dge()`](https://marusakod.github.io/edger.bp/reference/bp_prepare_dge.md).

- dispersion.method:

  Dispersion estimator to use when `dispersion` is `NULL`. `"stream"`
  uses
  [`bp_estimate_disp_stream()`](https://marusakod.github.io/edger.bp/reference/bp_estimate_disp_stream.md).
  `"dense"` uses
  [`bp_estimate_disp_dense()`](https://marusakod.github.io/edger.bp/reference/bp_estimate_disp_dense.md)
  as a small-data oracle.

- block.size:

  Number of rows to process per block.

- workers:

  Number of row blocks to process in parallel for supported streaming
  steps.

- max.elements:

  Maximum number of matrix elements allowed for dense dispersion
  estimation and dense AveLogCPM calculation.

- prior.count:

  Prior count passed to
  [`bp_glm_qlf()`](https://marusakod.github.io/edger.bp/reference/bp_glm_qlf.md).

- filter.args, dispersion.args, qlf.args:

  Additional arguments passed to
  [`bp_prepare_dge()`](https://marusakod.github.io/edger.bp/reference/bp_prepare_dge.md),
  the selected dispersion estimator, and
  [`bp_glm_qlf()`](https://marusakod.github.io/edger.bp/reference/bp_glm_qlf.md),
  respectively.

## Value

A list with `prepare`, `dispersion`, and `qlf` components.
