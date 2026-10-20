# edger.bp

Memory-efficient edgeR-compatible workflows for BPCells-backed count matrices.

Current high-level workflows are `bp_lrt()`, `bp_qlf()`, and `bp_treat()`.

## Agreement with edgeR

The [edgeR comparison vignette](https://marusakod.github.io/edger.bp/articles/edgeR-functionality.html)
runs native edgeR and edger.bp on the same simulated counts (115 genes and six
samples). The table below shows the maximum absolute differences for filtering,
normalization, dispersion estimation, and three statistical tests. These results
illustrate numerical agreement for this example under an audited dependency
stack; the vignette describes the supported inputs and workflows.

<!-- edger-comparison:start -->

|Quantity              |Maximum absolute difference |
|:---------------------|:---------------------------|
|filter keep vector    |0.00e+00                    |
|TMM norm factors      |0.00e+00                    |
|tagwise dispersion    |0.00e+00                    |
|LRT p-values          |0.00e+00                    |
|modern QLF p-values   |1.04e-13                    |
|modern Treat p-values |5.04e-14                    |

Generated with edger.bp 0.1.2, R 4.5.3, edgeR 4.8.2, limma 3.66.0, and BPCells 0.3.1.

<!-- edger-comparison:end -->

## Installation

### Conda

Create the audited R 4.4 environment:

```bash
micromamba create -n edger-bp \
  --override-channels \
  --channel https://raw.githubusercontent.com/wmacnair/edger.bp/conda-channel \
  --channel conda-forge \
  --channel bioconda \
  --strict-channel-priority \
  r-edger-bp=0.1.1=r44_0
```

An audited R 4.5 build is also available as `r-edger-bp=0.1.2=r45_0`.
See [the Conda channel documentation](https://github.com/marusakod/edger.bp/blob/main/conda/README.md) for details.

### R source package

Install the R package dependencies first:

```r
options(repos = c(
    bnprks = "https://bnprks.r-universe.dev",
    CRAN = "https://cloud.r-project.org"
))

if (!requireNamespace("BiocManager", quietly = TRUE)) {
    install.packages("BiocManager")
}

BiocManager::install("edgeR")
install.packages("BPCells")
```

Then install `edger.bp` from GitHub:

```r
if (!requireNamespace("remotes", quietly = TRUE)) {
    install.packages("remotes")
}

remotes::install_github("marusakod/edger.bp", build_vignettes = FALSE)
```

To install from a local checkout instead:

```r
remotes::install_local(".", build_vignettes = FALSE)
```

See [Using edger.bp](https://marusakod.github.io/edger.bp/articles/usage.html)
for a runnable package workflow.

See [edgeR functionality in edger.bp](https://marusakod.github.io/edger.bp/articles/edgeR-functionality.html)
for the complete test matrix and guidance on choosing between LRT, modern or legacy QLF, and
modern or legacy Treat.

## License

`edger.bp` is distributed under the GNU General Public License version 3.
