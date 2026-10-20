# Prepare BPCells Counts for edgeR GLM Workflows

Runs the canonical first-stage edgeR workflow on a BPCells-backed count
matrix: optional expression filtering, recomputation of post-filter
library sizes, and TMM normalization.

## Usage

``` r
bp_prepare_dge(
  y,
  design = NULL,
  group = NULL,
  filter = TRUE,
  norm.method = c("TMM", "none"),
  block.size = 1000,
  workers = 1,
  verbose = FALSE,
  ...
)
```

## Arguments

- y:

  A genes-by-samples count matrix.

- design:

  Optional design matrix passed to
  [`bp_filter_by_expr()`](https://marusakod.github.io/edger.bp/reference/bp_filter_by_expr.md).

- group:

  Optional group factor passed to
  [`bp_filter_by_expr()`](https://marusakod.github.io/edger.bp/reference/bp_filter_by_expr.md).

- filter:

  Logical scalar. If `TRUE`, apply
  [`bp_filter_by_expr()`](https://marusakod.github.io/edger.bp/reference/bp_filter_by_expr.md).

- norm.method:

  Normalization method passed to
  [`bp_calc_norm_factors()`](https://marusakod.github.io/edger.bp/reference/bp_calc_norm_factors.md).

- block.size:

  Number of rows to process per block while filtering.

- workers:

  Number of row blocks or normalization columns to process in parallel
  for supported steps.

- verbose:

  Logical scalar. If `TRUE`, emit low-noise progress messages from
  supported preparation steps.

- ...:

  Additional arguments passed to
  [`bp_filter_by_expr()`](https://marusakod.github.io/edger.bp/reference/bp_filter_by_expr.md).

## Value

A list with the filtered matrix, keep mask, library sizes, normalization
factors, effective library sizes, log offsets, and sample metadata.
