# Apply a block callback using the active worker pool

Uses the current PSOCK scope, or the usual numeric-worker backend
outside a scope. Only required callback bindings are sent to socket
workers.

## Usage

``` r
bp_block_lapply(x, workers, FUN, context = NULL)
```

## Arguments

- x:

  List or vector of tasks.

- workers:

  Number of workers, matching an active pool's size.

- FUN:

  Function of one task. Captured disk-backed matrices must equal
  `pool$counts`. Do not launch nested parallel jobs inside the callback.

- context:

  Optional label for error and progress messages.

## Value

A list of results in input order. Task errors include the task and
worker identity.
