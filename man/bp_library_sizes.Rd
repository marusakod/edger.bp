\name{bp_library_sizes}
\alias{bp_library_sizes}
\title{Compute Library Sizes}
\usage{
bp_library_sizes(y)
}
\arguments{
\item{y}{A genes-by-samples count matrix.}
}
\value{
A numeric vector of library sizes, one per sample.
}
\description{
Computes column sums for a count matrix without forcing the whole matrix to a
dense object.
}
