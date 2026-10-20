# Check Audited Package Versions

Reports whether the installed R, edgeR, limma, and BPCells versions
match one complete stack audited for exact streamed edgeR compatibility.

## Usage

``` r
bp_check_versions(error = FALSE)
```

## Arguments

- error:

  Logical scalar. If `TRUE`, throw an error when the installed versions
  do not match a complete audited stack. If `FALSE`, return the check
  table without error.

## Value

A data frame with package names, installed versions, expected versions,
individual version support, complete-stack match status, and the matched
stack name. A version can be individually supported while `ok` is
`FALSE` when the installed versions form an unaudited mixture.
