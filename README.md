# edger.bp Conda channel

This branch is an indexed Conda channel for tested `edger.bp` releases.

Channel URL:

```text
https://raw.githubusercontent.com/wmacnair/edger.bp/conda-channel
```

Create an R 4.4 environment:

```bash
micromamba create -n edger-bp-r44 \
  --override-channels \
  --channel https://raw.githubusercontent.com/wmacnair/edger.bp/conda-channel \
  --channel conda-forge \
  --channel bioconda \
  --strict-channel-priority \
  r-edger-bp=0.1.0=r44_0
```

Create an R 4.5 environment:

```bash
micromamba create -n edger-bp-r45 \
  --override-channels \
  --channel https://raw.githubusercontent.com/wmacnair/edger.bp/conda-channel \
  --channel conda-forge \
  --channel bioconda \
  --strict-channel-priority \
  r-edger-bp=0.1.0=r45_0
```

The channel contains the corresponding BPCells builds so both environments
solve correctly with strict channel priority.
