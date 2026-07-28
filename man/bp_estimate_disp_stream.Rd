\name{bp_estimate_disp_stream}
\alias{bp_estimate_disp_stream}
\title{Streamed Dispersion Estimation}
\usage{
bp_estimate_disp_stream(
  y,
  design,
  group = NULL,
  lib.size = NULL,
  offset = NULL,
  prior.df = NULL,
  trend.method = "locfit",
  tagwise = TRUE,
  span = NULL,
  legacy.span = FALSE,
  min.row.sum = 5,
  grid.length = 21,
  grid.range = c(-10, 10),
  robust = FALSE,
  winsor.tail.p = c(0.05, 0.1),
  tol = 1e-06,
  weights = NULL,
  block.size = 1000,
  workers = 1,
  verbose = FALSE,
  ...
)
}
\arguments{
\item{y}{A genes-by-samples count matrix.}
\item{design}{Design matrix. The streamed implementation currently supports
the GLM design path only.}
\item{group}{Optional group factor. Validated for consistency with edgeR but
otherwise unused when \code{design} is supplied.}
\item{lib.size}{Optional library sizes. If \code{NULL}, column sums of
\code{y} are used, matching \code{edgeR::estimateDisp.default()}.}
\item{offset}{Optional sample-level log-offset vector. Row-constant matrices
are accepted and canonicalized to vectors. Gene-specific offset matrices are
not yet supported.}
\item{prior.df, trend.method, tagwise, span, legacy.span, min.row.sum, grid.length, grid.range, robust, winsor.tail.p, tol}{Arguments corresponding to
\code{edgeR::estimateDisp()}. \code{robust = TRUE} is not yet supported by the
streamed path.}
\item{weights}{Observation weights. Currently unsupported.}
\item{block.size}{Number of rows to process per block.}
\item{workers}{Number of row blocks to process in parallel.}
\item{verbose}{Logical scalar. If \code{TRUE}, emit low-noise progress messages.}
\item{...}{Reserved for API compatibility with \code{edgeR::estimateDisp()}.}
}
\value{
A list matching the core dispersion fields returned by
\code{edgeR::estimateDisp.default()}.
}
\description{
Estimates edgeR dispersions for BPCells-backed matrices without materializing
the full count matrix. The implementation mirrors the audited
\code{edgeR::estimateDisp.default()} GLM design path for edgeR 4.10.1 by
building the adjusted-profile likelihood grid in dense row blocks.
}
