# Conda recipes

The recipes build `r-edger-bp` 0.1.1 from the tagged GitHub source archive.

- `r44` targets the audited R 4.4.3, edgeR 4.4.0, limma 3.62.0, and
  BPCells 0.3.1 stack. All dependencies are available from conda-forge and
  Bioconda.
- `r45` targets the audited R 4.5.3, edgeR 4.8.2, limma 3.66.0, and
  BPCells 0.3.1 stack. It requires an R 4.5 build of `r-bpcells` in an
  additional channel.

Example build:

```bash
conda-build conda/r44 \
  --override-channels \
  --channel conda-forge \
  --channel bioconda
```

Both recipes run `bp_check_versions(error = TRUE)` as their package test.
