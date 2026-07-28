\name{bp_treat}
\alias{bp_treat}
\title{End-to-End BPCells edgeR Treat Pipeline}
\usage{
bp_treat(
  y,
  design,
  group = NULL,
  coef = ncol(design),
  contrast = NULL,
  filter = TRUE,
  dispersion = NULL,
  norm.method = "TMM",
  dispersion.method = c("stream", "dense"),
  block.size = 1000,
  workers = 1,
  max.elements = 1e+07,
  prior.count = 0.125,
  lfc = log2(1.2),
  null = c("interval", "worst.case"),
  filter.args = list(),
  dispersion.args = list(),
  treat.args = list()
)
}
\arguments{
\item{y}{A genes-by-samples count matrix.}
\item{design}{Design matrix.}
\item{group}{Optional group factor for filtering and dispersion estimation.}
\item{coef, contrast}{Test specification passed to \code{bp_glm_treat()}.}
\item{filter}{Logical scalar. If \code{TRUE}, apply \code{bp_filter_by_expr()}.}
\item{dispersion}{Optional scalar or gene-level dispersion vector. If
\code{NULL}, dispersions are estimated with \code{dispersion.method}.}
\item{norm.method}{Normalization method passed to \code{bp_prepare_dge()}.}
\item{dispersion.method}{Dispersion estimator to use when \code{dispersion} is
\code{NULL}. \code{"stream"} uses \code{bp_estimate_disp_stream()}.
\code{"dense"} uses \code{bp_estimate_disp_dense()} as a small-data oracle.}
\item{block.size}{Number of rows to process per block.}
\item{workers}{Number of row blocks to process in parallel for supported
streaming steps.}
\item{max.elements}{Maximum number of matrix elements allowed for dense
dispersion estimation and dense AveLogCPM calculation.}
\item{prior.count}{Prior count passed to \code{bp_glm_treat()}.}
\item{lfc}{Log2 fold-change threshold passed to \code{bp_glm_treat()}.}
\item{null}{Null hypothesis passed to \code{bp_glm_treat()}.}
\item{filter.args, dispersion.args, treat.args}{Additional arguments passed to
\code{bp_prepare_dge()}, the selected dispersion estimator, and
\code{bp_glm_treat()}, respectively.}
}
\value{
A list with \code{prepare}, \code{dispersion}, and \code{treat} components.
}
\description{
Runs the supported BPCells edgeR quasi-likelihood Treat workflow: preparation,
streamed or dense-reference dispersion estimation when needed, and blockwise
Treat testing. The modern QL path is used by default; pass
\code{treat.args = list(legacy = TRUE)} for legacy QL Treat comparisons.
}
