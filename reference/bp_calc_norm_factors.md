# Streaming TMM Normalization Factors

Computes TMM normalization factors for a matrix-like object while
avoiding a full dense realization of the genes-by-samples count matrix.

## Usage

``` r
bp_calc_norm_factors(
  y,
  lib.size = NULL,
  method = c("TMM", "none"),
  refColumn = NULL,
  logratioTrim = 0.3,
  sumTrim = 0.05,
  doWeighting = TRUE,
  Acutoff = -1e+10,
  workers = 1,
  verbose = FALSE
)
```

## Arguments

- y:

  A genes-by-samples count matrix.

- lib.size:

  Optional library sizes.

- method:

  Normalization method. Currently `"TMM"` or `"none"`.

- refColumn:

  Optional 1-based reference column.

- logratioTrim, sumTrim, doWeighting, Acutoff:

  TMM parameters with the same meaning as in
  [`edgeR::normLibSizes()`](https://rdrr.io/pkg/edgeR/man/calcNormFactors.html).

- workers:

  Number of columns to process in parallel for TMM.

- verbose:

  Logical scalar. If `TRUE`, emit low-noise progress messages.

## Value

A numeric vector of normalization factors, one per sample.
