# Dense Reference Dispersion Estimation

Estimates edgeR dispersions after materializing a BPCells-backed matrix
as a dense R matrix. This is intended as a correctness baseline and
small-data implementation for checking the streamed likelihood-grid
implementation.

## Usage

``` r
bp_estimate_disp_dense(
  y,
  design = NULL,
  group = NULL,
  lib.size = NULL,
  offset = NULL,
  prior.df = NULL,
  trend.method = "locfit",
  tagwise = TRUE,
  span = NULL,
  legacy.span = FALSE,
  min.row.sum = 5,
  grid.length = 21,
  grid.range = c(-10, 10),
  robust = FALSE,
  winsor.tail.p = c(0.05, 0.1),
  tol = 1e-06,
  max.elements = 1e+07,
  ...
)
```

## Arguments

- y:

  A genes-by-samples count matrix.

- design:

  Optional design matrix.

- group:

  Optional group factor.

- lib.size:

  Optional library sizes used when `offset` is `NULL`.

- offset:

  Optional sample-level log-offset vector. Row-constant matrices are
  accepted and canonicalized to vectors. Gene-specific offset matrices
  are not yet supported.

- prior.df, trend.method, tagwise, span, legacy.span, min.row.sum,
  grid.length, grid.range, robust, winsor.tail.p, tol:

  Arguments passed to
  [`edgeR::estimateDisp()`](https://rdrr.io/pkg/edgeR/man/estimateDisp.html).

- max.elements:

  Maximum number of matrix elements allowed before densifying.

- ...:

  Additional arguments passed to
  [`edgeR::estimateDisp()`](https://rdrr.io/pkg/edgeR/man/estimateDisp.html).

## Value

The list returned by
[`edgeR::estimateDisp.default()`](https://rdrr.io/pkg/edgeR/man/estimateDisp.html).
