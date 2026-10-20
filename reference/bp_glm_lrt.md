# Blockwise edgeR GLM Likelihood Ratio Test

Fits negative-binomial GLMs and runs likelihood ratio tests on
BPCells-backed counts by processing dense row blocks with public edgeR
wrappers.

## Usage

``` r
bp_glm_lrt(
  y,
  design,
  dispersion,
  offset = NULL,
  lib.size = NULL,
  coef = ncol(design),
  contrast = NULL,
  prior.count = 0.125,
  block.size = 1000,
  workers = 1,
  keep.fitted = FALSE,
  ave.log.cpm = NULL,
  ...
)
```

## Arguments

- y:

  A genes-by-samples count matrix.

- design:

  Design matrix.

- dispersion:

  Scalar or gene-level dispersion vector.

- offset:

  Optional sample-level log-offset vector. Row-constant matrices are
  accepted and canonicalized to vectors. Gene-specific offset matrices
  are not yet supported.

- lib.size:

  Optional library sizes used when `offset` is `NULL`.

- coef, contrast:

  Test specification passed to
  [`edgeR::glmLRT()`](https://rdrr.io/pkg/edgeR/man/glmLRT.html).

- prior.count:

  Prior count passed to
  [`edgeR::glmFit()`](https://rdrr.io/pkg/edgeR/man/glmfit.html).

- block.size:

  Number of rows to process per block.

- workers:

  Number of row blocks to process in parallel.

- keep.fitted:

  If `TRUE`, keep fitted values in the returned object.

- ave.log.cpm:

  Optional precomputed average log CPM vector. If supplied, it replaces
  the blockwise `logCPM` column from
  [`edgeR::glmLRT()`](https://rdrr.io/pkg/edgeR/man/glmLRT.html).

- ...:

  Additional arguments passed to
  [`edgeR::glmFit()`](https://rdrr.io/pkg/edgeR/man/glmfit.html).

## Value

A list containing the combined LRT table, GLM coefficients, deviance,
residual degrees of freedom, dispersions, and test metadata.
