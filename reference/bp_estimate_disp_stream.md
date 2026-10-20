# Streamed Dispersion Estimation

Estimates edgeR dispersions for BPCells-backed matrices without
materializing the full count matrix. The implementation mirrors the
audited
[`edgeR::estimateDisp.default()`](https://rdrr.io/pkg/edgeR/man/estimateDisp.html)
GLM design path for edgeR 4.10.1 by building the adjusted-profile
likelihood grid in dense row blocks.

## Usage

``` r
bp_estimate_disp_stream(
  y,
  design,
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
  weights = NULL,
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

  Design matrix. The streamed implementation currently supports the GLM
  design path only.

- group:

  Optional group factor. Validated for consistency with edgeR but
  otherwise unused when `design` is supplied.

- lib.size:

  Optional library sizes. If `NULL`, column sums of `y` are used,
  matching
  [`edgeR::estimateDisp.default()`](https://rdrr.io/pkg/edgeR/man/estimateDisp.html).

- offset:

  Optional sample-level log-offset vector. Row-constant matrices are
  accepted and canonicalized to vectors. Gene-specific offset matrices
  are not yet supported.

- prior.df, trend.method, tagwise, span, legacy.span, min.row.sum,
  grid.length, grid.range, robust, winsor.tail.p, tol:

  Arguments corresponding to
  [`edgeR::estimateDisp()`](https://rdrr.io/pkg/edgeR/man/estimateDisp.html).
  `robust = TRUE` is not yet supported by the streamed path.

- weights:

  Observation weights. Currently unsupported.

- block.size:

  Number of rows to process per block.

- workers:

  Number of row blocks to process in parallel.

- verbose:

  Logical scalar. If `TRUE`, emit low-noise progress messages.

- ...:

  Reserved for API compatibility with
  [`edgeR::estimateDisp()`](https://rdrr.io/pkg/edgeR/man/estimateDisp.html).

## Value

A list matching the core dispersion fields returned by
[`edgeR::estimateDisp.default()`](https://rdrr.io/pkg/edgeR/man/estimateDisp.html).
