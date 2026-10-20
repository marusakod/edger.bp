# Ensure BPCells Storage Order

Returns a BPCells matrix with the requested physical storage order. If
the input already has that order, it is returned unchanged. Otherwise, a
companion BPCells directory is written with
[`BPCells::transpose_storage_order()`](https://bnprks.github.io/BPCells/reference/transpose_storage_order.html).

## Usage

``` r
bp_ensure_storage_order(
  y,
  order = c("row", "col"),
  outdir = NULL,
  tmpdir = tempdir(),
  load_bytes = 4194304L,
  sort_bytes = 1073741824L
)
```

## Arguments

- y:

  A BPCells-backed matrix.

- order:

  Desired physical storage order, either `"row"` or `"col"`.

- outdir:

  Directory where the companion BPCells matrix should be written when
  conversion is required. Required if `y` is not already in the desired
  order.

- tmpdir:

  Temporary directory passed to
  [`BPCells::transpose_storage_order()`](https://bnprks.github.io/BPCells/reference/transpose_storage_order.html).

- load_bytes, sort_bytes:

  Tuning arguments passed to
  [`BPCells::transpose_storage_order()`](https://bnprks.github.io/BPCells/reference/transpose_storage_order.html).

## Value

A BPCells matrix with the requested physical storage order.
