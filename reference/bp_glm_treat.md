# Blockwise edgeR GLM Treat Test

Fits quasi-likelihood GLMs and runs Treat tests on BPCells-backed counts
by processing dense row blocks while preserving edgeR's global empirical
Bayes shrinkage step. Both edgeR's legacy QL path (`legacy = TRUE`) and
modern QL path (`legacy = FALSE`) are supported.

## Usage

``` r
bp_glm_treat(
  y,
  design,
  dispersion,
  offset = NULL,
  lib.size = NULL,
  coef = ncol(design),
  contrast = NULL,
  prior.count = 0.125,
  abundance.trend = TRUE,
  ave.log.cpm = NULL,
  block.size = 1000,
  workers = 1,
  robust = FALSE,
  winsor.tail.p = c(0.05, 0.1),
  lfc = log2(1.2),
  null = c("interval", "worst.case"),
  legacy = FALSE,
  weights = NULL,
  keep.unit.mat = FALSE,
  ...
)
```

## Arguments

- y:

  A genes-by-samples count matrix.

- design:

  Design matrix.

- dispersion:

  Scalar or gene-level negative-binomial dispersion vector.

- offset:

  Optional sample-level log-offset vector. Row-constant matrices are
  accepted and canonicalized to vectors. Gene-specific offset matrices
  are not yet supported.

- lib.size:

  Optional library sizes used by
  [`edgeR::glmFit()`](https://rdrr.io/pkg/edgeR/man/glmfit.html) when
  `offset` is `NULL`.

- coef, contrast:

  Test specification passed to
  [`edgeR::glmTreat()`](https://rdrr.io/pkg/edgeR/man/glmTreat.html).

- prior.count:

  Prior count passed to
  [`edgeR::glmFit()`](https://rdrr.io/pkg/edgeR/man/glmfit.html).

- abundance.trend:

  Logical scalar passed to the QL empirical Bayes trend.

- ave.log.cpm:

  Optional precomputed average log CPM vector.

- block.size:

  Number of rows to process per block.

- workers:

  Number of row blocks to process in parallel.

- robust, winsor.tail.p:

  Robust empirical Bayes options. `robust = TRUE` is not yet supported.

- lfc:

  Log2 fold-change threshold passed to
  [`edgeR::glmTreat()`](https://rdrr.io/pkg/edgeR/man/glmTreat.html).

- null:

  Null hypothesis passed to
  [`edgeR::glmTreat()`](https://rdrr.io/pkg/edgeR/man/glmTreat.html).

- legacy:

  Logical scalar. If `FALSE`, use edgeR's modern QL adjustment with
  streamed computation of the global average QL dispersion.

- weights:

  Observation weights. Currently unsupported.

- keep.unit.mat:

  Unit-level QL diagnostic matrices are not yet supported.

- ...:

  Additional arguments passed to
  [`edgeR::glmFit()`](https://rdrr.io/pkg/edgeR/man/glmfit.html).

## Value

A list containing the combined Treat table, GLM coefficients,
quasi-likelihood shrinkage fields, and test metadata.
