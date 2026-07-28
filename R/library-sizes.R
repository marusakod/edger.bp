#' Compute Library Sizes
#'
#' Computes column sums for a count matrix without forcing the whole matrix to a
#' dense object.
#'
#' @param y A genes-by-samples count matrix. BPCells matrices, sparse matrices,
#'   and ordinary matrix-like objects are supported.
#'
#' @return A numeric vector of library sizes, one per sample.
#' @export
bp_library_sizes <- function(y) {
    .check_count_matrix(y)
    out <- .matrix_col_sums(y)
    out <- as.numeric(out)
    names(out) <- colnames(y)
    out
}
