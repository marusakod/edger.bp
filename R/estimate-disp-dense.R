#' Dense Reference Dispersion Estimation
#'
#' Estimates edgeR dispersions after materializing a BPCells-backed matrix as a
#' dense R matrix. This is intended as a correctness baseline and small-data
#' implementation for checking the streamed likelihood-grid implementation.
#'
#' @param y A genes-by-samples count matrix.
#' @param design Optional design matrix.
#' @param group Optional group factor.
#' @param lib.size Optional library sizes. Used by `edgeR::estimateDisp()` when
#'   `offset` is `NULL`.
#' @param offset Optional sample-level log-offset vector of length `ncol(y)`.
#'   Row-constant matrices with dimensions `nrow(y)` by `ncol(y)` are accepted
#'   and canonicalized to vectors. Gene-specific offset matrices are not yet
#'   supported.
#' @param prior.df,trend.method,tagwise,span,legacy.span,min.row.sum,grid.length,grid.range,robust,winsor.tail.p,tol
#'   Arguments passed to `edgeR::estimateDisp()`.
#' @param max.elements Maximum number of matrix elements allowed before
#'   densifying. This guard prevents accidental full-scale dense conversion.
#' @param ... Additional arguments passed to `edgeR::estimateDisp()`.
#'
#' @return The list returned by `edgeR::estimateDisp.default()`.
#' @export
bp_estimate_disp_dense <- function(y,
                                   design = NULL,
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
                                   tol = 1e-6,
                                   max.elements = 1e7,
                                   ...) {
    dims <- .check_count_matrix(y)
    ntags <- dims[1L]
    nlibs <- dims[2L]
    if (ntags * nlibs > max.elements) {
        stop(
            "`y` has more than `max.elements` entries; dense dispersion estimation is only a reference/small-data path.",
            call. = FALSE
        )
    }

    lib.size <- .check_lib_size(lib.size, nlibs)
    offset <- .normalize_offset(offset, ntags, nlibs)
    if (!is.null(design)) {
        design <- as.matrix(design)
        if (nrow(design) != nlibs) {
            stop("`design` must have one row per sample.", call. = FALSE)
        }
    }

    counts <- .dense_slice(y, TRUE, TRUE)
    storage.mode(counts) <- "numeric"

    out <- edgeR::estimateDisp(
        counts,
        design = design,
        group = group,
        lib.size = lib.size,
        offset = offset,
        prior.df = prior.df,
        trend.method = trend.method,
        tagwise = tagwise,
        span = span,
        legacy.span = legacy.span,
        min.row.sum = min.row.sum,
        grid.length = grid.length,
        grid.range = grid.range,
        robust = robust,
        winsor.tail.p = winsor.tail.p,
        tol = tol,
        ...
    )
    class(out) <- c("BPCellsDenseDisp", class(out))
    out
}
