\name{bp_glm_treat}
\alias{bp_glm_treat}
\title{Blockwise edgeR GLM Treat Test}
\usage{
bp_glm_treat(
  y,
  design,
  dispersion,
  offset = NULL,
  lib.size = NULL,
  coef = ncol(design),
  contrast = NULL,
  prior.count = 0.125,
  abundance.trend = TRUE,
  ave.log.cpm = NULL,
  block.size = 1000,
  workers = 1,
  robust = FALSE,
  winsor.tail.p = c(0.05, 0.1),
  lfc = log2(1.2),
  null = c("interval", "worst.case"),
  legacy = FALSE,
  weights = NULL,
  keep.unit.mat = FALSE,
  ...
)
}
\arguments{
\item{y}{A genes-by-samples count matrix.}
\item{design}{Design matrix.}
\item{dispersion}{Scalar or gene-level negative-binomial dispersion vector.}
\item{offset}{Optional sample-level log-offset vector. Row-constant matrices
are accepted and canonicalized to vectors. Gene-specific offset matrices are
not yet supported.}
\item{lib.size}{Optional library sizes used by \code{edgeR::glmFit()} when
\code{offset} is \code{NULL}.}
\item{coef, contrast}{Test specification passed to \code{edgeR::glmTreat()}.}
\item{prior.count}{Prior count passed to \code{edgeR::glmFit()}.}
\item{abundance.trend}{Logical scalar passed to the QL empirical Bayes trend.}
\item{ave.log.cpm}{Optional precomputed average log CPM vector.}
\item{block.size}{Number of rows to process per block.}
\item{workers}{Number of row blocks to process in parallel.}
\item{robust, winsor.tail.p}{Robust empirical Bayes options.
\code{robust = TRUE} is not yet supported.}
\item{lfc}{Log2 fold-change threshold passed to \code{edgeR::glmTreat()}.}
\item{null}{Null hypothesis passed to \code{edgeR::glmTreat()}.}
\item{legacy}{Logical scalar. If \code{FALSE}, use edgeR's modern QL
adjustment with streamed computation of the global average QL dispersion.}
\item{weights}{Observation weights. Currently unsupported.}
\item{keep.unit.mat}{Unit-level QL diagnostic matrices are not yet supported.}
\item{...}{Additional arguments passed to \code{edgeR::glmFit()}.}
}
\value{
A list containing the combined Treat table, GLM coefficients,
quasi-likelihood shrinkage fields, and test metadata.
}
\description{
Fits quasi-likelihood GLMs and runs Treat tests on BPCells-backed counts by
processing dense row blocks while preserving edgeR's global empirical Bayes
shrinkage step. Both edgeR's legacy QL path (\code{legacy = TRUE}) and modern QL
path (\code{legacy = FALSE}) are supported.
}
