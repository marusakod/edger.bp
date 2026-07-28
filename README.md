# edger.bp

Memory-efficient edgeR-compatible workflows for BPCells-backed count matrices.

Current high-level workflows are `bp_lrt()`, `bp_qlf()`, and `bp_treat()`.

## Installation

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

remotes::install_github("wmacnair/edger.bp", build_vignettes = FALSE)
```

To install from a local checkout instead:

```r
remotes::install_local(".", build_vignettes = FALSE)
```

See `vignette("usage", package = "edger.bp")` for a runnable package workflow.

See `vignette("edgeR-functionality", package = "edger.bp")` for the complete
test matrix and guidance on choosing between LRT, modern or legacy QLF, and
modern or legacy Treat.

## License

`edger.bp` is distributed under the GNU General Public License version 3.
