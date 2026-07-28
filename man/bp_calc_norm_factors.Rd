\name{bp_calc_norm_factors}
\alias{bp_calc_norm_factors}
\title{Streaming TMM Normalization Factors}
\usage{
bp_calc_norm_factors(
  y,
  lib.size = NULL,
  method = c("TMM", "none"),
  refColumn = NULL,
  logratioTrim = 0.3,
  sumTrim = 0.05,
  doWeighting = TRUE,
  Acutoff = -1e+10,
  workers = 1,
  verbose = FALSE
)
}
\arguments{
\item{y}{A genes-by-samples count matrix.}
\item{lib.size}{Optional library sizes.}
\item{method}{Normalization method. Currently \code{"TMM"} or \code{"none"}.}
\item{refColumn}{Optional 1-based reference column.}
\item{logratioTrim, sumTrim, doWeighting, Acutoff}{TMM parameters with the same meaning as in \code{edgeR::normLibSizes()}.}
\item{workers}{Number of columns to process in parallel for TMM.}
\item{verbose}{Logical scalar. If \code{TRUE}, emit low-noise progress messages.}
}
\value{
A numeric vector of normalization factors, one per sample.
}
\description{
Computes TMM normalization factors for a matrix-like object while avoiding a
full dense realization of the genes-by-samples count matrix.
}
