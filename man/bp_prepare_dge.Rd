\name{bp_prepare_dge}
\alias{bp_prepare_dge}
\title{Prepare BPCells Counts for edgeR GLM Workflows}
\usage{
bp_prepare_dge(
  y,
  design = NULL,
  group = NULL,
  filter = TRUE,
  norm.method = c("TMM", "none"),
  block.size = 1000,
  workers = 1,
  verbose = FALSE,
  ...
)
}
\arguments{
\item{y}{A genes-by-samples count matrix.}
\item{design}{Optional design matrix passed to \code{bp_filter_by_expr()}.}
\item{group}{Optional group factor passed to \code{bp_filter_by_expr()}.}
\item{filter}{Logical scalar. If \code{TRUE}, apply \code{bp_filter_by_expr()}.}
\item{norm.method}{Normalization method passed to \code{bp_calc_norm_factors()}.}
\item{block.size}{Number of rows to process per block while filtering.}
\item{workers}{Number of row blocks or normalization columns to process in
parallel for supported steps.}
\item{verbose}{Logical scalar. If \code{TRUE}, emit low-noise progress messages
from supported preparation steps.}
\item{...}{Additional arguments passed to \code{bp_filter_by_expr()}.}
}
\value{
A list with the filtered matrix, keep mask, library sizes, normalization
factors, effective library sizes, log offsets, and sample metadata.
}
\description{
Runs the canonical first-stage edgeR workflow on a BPCells-backed count matrix:
optional expression filtering, recomputation of post-filter library sizes, and
TMM normalization.
}
