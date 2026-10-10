\name{bp_start_psock_pool}
\alias{bp_start_psock_pool}
\alias{bp_stop_psock_pool}
\alias{bp_with_psock_pool}
\title{Reusable PSOCK workers for disk-backed fitting}
\usage{
bp_start_psock_pool(matrix_dir, rows = TRUE, cols = TRUE,
                    workers = 1L, verbose = TRUE, timeout = 3600)
bp_stop_psock_pool(pool)
bp_with_psock_pool(pool, code)
}
\arguments{
\item{matrix_dir}{Native BPCells matrix directory, readable by every worker.}
\item{rows, cols}{Row and column indices in analysis order.}
\item{workers}{Number of independent R sessions.}
\item{verbose}{Print startup, pass and shutdown diagnostics.}
\item{timeout}{Socket timeout in seconds.}
\item{pool}{Pool returned by \code{bp_start_psock_pool}.}
\item{code}{Expression containing fits of \code{pool$counts}, using the pool's worker count.}
}
\description{Reuse one group of independent R sessions across streamed dispersion
estimation and multiple comparisons, avoiding repeated forks of the parent R session.}
\details{Each worker opens the matrix and selected rows/columns once. Count
matrices are not exported as dense data. Register
\code{on.exit(bp_stop_psock_pool(pool), add = TRUE)} immediately after creating
a pool inside a function. A pool belongs to its creating process; scopes cannot
be nested. Worker library threads are restricted to one. Numeric workers outside
a pool scope retain their existing serial/fork behaviour. A pool does not bypass
operating-system process limits. Callbacks referencing other matrices are rejected.
The timeout must allow the longest block/pass communication to finish.}
\value{The constructor returns a pool environment containing \code{counts}.
The scope returns the result of \code{code}; stopping returns invisible NULL.}
\examples{
\dontrun{
run_comparisons <- function(matrix_dir, design, offset) {
  pool <- bp_start_psock_pool(matrix_dir, workers = 4)
  on.exit(bp_stop_psock_pool(pool), add = TRUE)
  bp_with_psock_pool(pool, {
    dispersion <- bp_estimate_disp_stream(pool$counts, design,
      offset = offset, workers = 4)
    bp_glm_treat(pool$counts, design,
      dispersion = dispersion$tagwise.dispersion, offset = offset, workers = 4)
  })
}
}
}
