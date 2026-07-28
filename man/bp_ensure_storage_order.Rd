\name{bp_ensure_storage_order}
\alias{bp_ensure_storage_order}
\title{Ensure BPCells Storage Order}
\usage{
bp_ensure_storage_order(
  y,
  order = c("row", "col"),
  outdir = NULL,
  tmpdir = tempdir(),
  load_bytes = 4194304L,
  sort_bytes = 1073741824L
)
}
\arguments{
\item{y}{A BPCells-backed matrix.}
\item{order}{Desired physical storage order, either \code{"row"} or
\code{"col"}.}
\item{outdir}{Directory where the companion BPCells matrix should be written
when conversion is required. Required if \code{y} is not already in the
desired order.}
\item{tmpdir}{Temporary directory passed to
\code{BPCells::transpose_storage_order()}.}
\item{load_bytes, sort_bytes}{Tuning arguments passed to
\code{BPCells::transpose_storage_order()}.}
}
\value{
A BPCells matrix with the requested physical storage order.
}
\description{
Returns a BPCells matrix with the requested physical storage order. If the input
already has that order, it is returned unchanged. Otherwise, a companion BPCells
directory is written with \code{BPCells::transpose_storage_order()}.
}
