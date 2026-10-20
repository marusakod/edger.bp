# Compute Library Sizes

Computes column sums for a count matrix without forcing the whole matrix
to a dense object.

## Usage

``` r
bp_library_sizes(y)
```

## Arguments

- y:

  A genes-by-samples count matrix.

## Value

A numeric vector of library sizes, one per sample.
