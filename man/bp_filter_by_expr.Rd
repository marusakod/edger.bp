\name{bp_filter_by_expr}
\alias{bp_filter_by_expr}
\title{Streaming edgeR-Compatible Expression Filter}
\usage{
bp_filter_by_expr(
  y,
  design = NULL,
  group = NULL,
  lib.size = NULL,
  min.count = 10,
  min.total.count = 15,
  large.n = 10,
  min.prop = 0.7,
  block.size = 1000,
  workers = 1
)
}
\arguments{
\item{y}{A genes-by-samples count matrix.}
\item{design}{Optional design matrix.}
\item{group}{Optional group factor.}
\item{lib.size}{Optional library sizes.}
\item{min.count, min.total.count, large.n, min.prop}{Filtering parameters.}
\item{block.size}{Number of rows to process per block.}
\item{workers}{Number of row blocks to process in parallel.}
}
\value{
A logical vector indicating retained rows.
}
\description{
Reproduces the default logic of \code{edgeR::filterByExpr()} while processing
counts in row blocks.
}
