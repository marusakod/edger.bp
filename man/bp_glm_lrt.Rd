\name{bp_glm_lrt}
\alias{bp_glm_lrt}
\title{Blockwise edgeR GLM Likelihood Ratio Test}
\usage{
bp_glm_lrt(
  y,
  design,
  dispersion,
  offset = NULL,
  lib.size = NULL,
  coef = ncol(design),
  contrast = NULL,
  prior.count = 0.125,
  block.size = 1000,
  workers = 1,
  keep.fitted = FALSE,
  ave.log.cpm = NULL,
  ...
)
}
\arguments{
\item{y}{A genes-by-samples count matrix.}
\item{design}{Design matrix.}
\item{dispersion}{Scalar or gene-level dispersion vector.}
\item{offset}{Optional sample-level log-offset vector. Row-constant matrices
are accepted and canonicalized to vectors. Gene-specific offset matrices are
not yet supported.}
\item{lib.size}{Optional library sizes used when \code{offset} is \code{NULL}.}
\item{coef, contrast}{Test specification passed to \code{edgeR::glmLRT()}.}
\item{prior.count}{Prior count passed to \code{edgeR::glmFit()}.}
\item{block.size}{Number of rows to process per block.}
\item{workers}{Number of row blocks to process in parallel.}
\item{keep.fitted}{If \code{TRUE}, keep fitted values in the returned object.}
\item{ave.log.cpm}{Optional precomputed average log CPM vector. If supplied,
it replaces the blockwise \code{logCPM} column from \code{edgeR::glmLRT()}.}
\item{...}{Additional arguments passed to \code{edgeR::glmFit()}.}
}
\value{
A list containing the combined LRT table, GLM coefficients, deviance,
residual degrees of freedom, dispersions, and test metadata.
}
\description{
Fits negative-binomial GLMs and runs likelihood ratio tests on BPCells-backed
counts by processing dense row blocks with public edgeR wrappers.
}
