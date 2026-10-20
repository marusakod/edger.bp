# Streaming edgeR-Compatible Expression Filter

Reproduces the default logic of
[`edgeR::filterByExpr()`](https://rdrr.io/pkg/edgeR/man/filterByExpr.html)
while processing counts in row blocks.

## Usage

``` r
bp_filter_by_expr(
  y,
  design = NULL,
  group = NULL,
  lib.size = NULL,
  min.count = 10,
  min.total.count = 15,
  large.n = 10,
  min.prop = 0.7,
  block.size = 1000,
  workers = 1
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

  Optional library sizes.

- min.count, min.total.count, large.n, min.prop:

  Filtering parameters.

- block.size:

  Number of rows to process per block.

- workers:

  Number of row blocks to process in parallel.

## Value

A logical vector indicating retained rows.
