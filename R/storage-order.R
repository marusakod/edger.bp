#' Ensure BPCells Storage Order
#'
#' Returns a BPCells matrix with the requested physical storage order. If the
#' input already has that order, it is returned unchanged. Otherwise, a companion
#' BPCells directory is written with `BPCells::transpose_storage_order()`.
#'
#' @param y A BPCells-backed matrix.
#' @param order Desired physical storage order, either `"row"` or `"col"`.
#' @param outdir Directory where the companion BPCells matrix should be written
#'   when conversion is required. Required if `y` is not already in the desired
#'   order.
#' @param tmpdir Temporary directory passed to
#'   `BPCells::transpose_storage_order()`.
#' @param load_bytes,sort_bytes Tuning arguments passed to
#'   `BPCells::transpose_storage_order()`.
#'
#' @return A BPCells matrix with the requested physical storage order.
#' @export
bp_ensure_storage_order <- function(y,
                                    order = c("row", "col"),
                                    outdir = NULL,
                                    tmpdir = tempdir(),
                                    load_bytes = 4194304L,
                                    sort_bytes = 1073741824L) {
    order <- match.arg(order)
    current <- BPCells::storage_order(y)
    if (identical(current, order)) {
        return(y)
    }
    if (!identical(current, setdiff(c("row", "col"), order))) {
        stop("`y` must have BPCells storage order \"row\" or \"col\".", call. = FALSE)
    }
    if (is.null(outdir) || length(outdir) != 1L || is.na(outdir) || !nzchar(outdir)) {
        stop("`outdir` is required when writing a BPCells storage-order companion.", call. = FALSE)
    }
    if (file.exists(outdir)) {
        stop("`outdir` already exists; choose a new directory for the BPCells companion.", call. = FALSE)
    }

    BPCells::transpose_storage_order(
        y,
        outdir = outdir,
        tmpdir = tmpdir,
        load_bytes = load_bytes,
        sort_bytes = sort_bytes
    )
}
