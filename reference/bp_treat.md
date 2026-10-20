# End-to-End BPCells edgeR Treat Pipeline

Runs the supported BPCells edgeR quasi-likelihood Treat workflow:
preparation, streamed or dense-reference dispersion estimation when
needed, and blockwise Treat testing. The modern QL path is used by
default; pass `treat.args = list(legacy = TRUE)` for legacy QL Treat
comparisons.

## Usage

``` r
bp_treat(
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
  lfc = log2(1.2),
  null = c("interval", "worst.case"),
  filter.args = list(),
  dispersion.args = list(),
  treat.args = list()
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
  [`bp_glm_treat()`](https://marusakod.github.io/edger.bp/reference/bp_glm_treat.md).

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
  [`bp_glm_treat()`](https://marusakod.github.io/edger.bp/reference/bp_glm_treat.md).

- lfc:

  Log2 fold-change threshold passed to
  [`bp_glm_treat()`](https://marusakod.github.io/edger.bp/reference/bp_glm_treat.md).

- null:

  Null hypothesis passed to
  [`bp_glm_treat()`](https://marusakod.github.io/edger.bp/reference/bp_glm_treat.md).

- filter.args, dispersion.args, treat.args:

  Additional arguments passed to
  [`bp_prepare_dge()`](https://marusakod.github.io/edger.bp/reference/bp_prepare_dge.md),
  the selected dispersion estimator, and
  [`bp_glm_treat()`](https://marusakod.github.io/edger.bp/reference/bp_glm_treat.md),
  respectively.

## Value

A list with `prepare`, `dispersion`, and `treat` components.
