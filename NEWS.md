# edger.bp 0.1.2

- Export reusable PSOCK worker pools for streamed dispersion, LRT, QLF and
  Treat fitting. Workers open their own BPCells matrices and reuse them across
  passes, with BLAS and OpenMP restricted to one thread.
- Transfer only callback bindings that workers need; reopen recognized
  BPCells counts in each worker rather than serialize external pointers.
- Check that workers load the same package installation as the parent and
  report callback errors with their pass context.

# edger.bp 0.1.1

- Reuse the fitted modern quasi-likelihood alternative model during Treat
  testing, avoiding one redundant blockwise GLM fit while preserving results.

# edger.bp 0.1.0

- Initial public release.
