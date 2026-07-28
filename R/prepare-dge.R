#' Prepare BPCells Counts for edgeR GLM Workflows
#'
#' Runs the canonical first-stage edgeR workflow on a BPCells-backed count
#' matrix: optional expression filtering, recomputation of post-filter library
#' sizes, and TMM normalization.
#'
#' @param y A genes-by-samples count matrix.
#' @param design Optional design matrix passed to [bp_filter_by_expr()].
#' @param group Optional group factor passed to [bp_filter_by_expr()].
#' @param filter Logical scalar. If `TRUE`, apply [bp_filter_by_expr()].
#' @param norm.method Normalization method passed to [bp_calc_norm_factors()].
#'   Currently `"TMM"` or `"none"`.
#' @param block.size Number of rows to process per block while filtering.
#' @param workers Number of row blocks or normalization columns to process in
#'   parallel for supported steps.
#' @param verbose Logical scalar. If `TRUE`, emit low-noise progress messages
#'   from supported preparation steps.
#' @param ... Additional arguments passed to [bp_filter_by_expr()].
#'
#' @return A list with the filtered matrix, keep mask, library sizes,
#'   normalization factors, effective library sizes, log offsets, and sample
#'   metadata.
#' @export
bp_prepare_dge <- function(y,
                           design = NULL,
                           group = NULL,
                           filter = TRUE,
                           norm.method = c("TMM", "none"),
                           block.size = 1000,
                           workers = 1,
                           verbose = FALSE,
                           ...) {
    dims <- .check_count_matrix(y)
    norm.method <- match.arg(norm.method)
    verbose <- .check_verbose(verbose)

    if (isTRUE(filter)) {
        .verbose_message(verbose, "DGE preparation: filtering genes")
        keep <- bp_filter_by_expr(
            y,
            design = design,
            group = group,
            block.size = block.size,
            workers = workers,
            ...
        )
    } else {
        keep <- rep(TRUE, dims[1L])
        names(keep) <- rownames(y)
    }

    .verbose_message(verbose, "DGE preparation: subsetting to retained genes")
    y.filtered <- y[which(keep), , drop = FALSE]
    .verbose_message(verbose, "DGE preparation: computing library sizes")
    lib.size <- bp_library_sizes(y.filtered)
    .verbose_message(verbose, "DGE preparation: computing ", norm.method, " normalization factors")
    norm.factors <- bp_calc_norm_factors(
        y.filtered,
        lib.size = lib.size,
        method = norm.method,
        workers = workers,
        verbose = verbose
    )
    effective.lib.size <- lib.size * norm.factors
    offset <- log(effective.lib.size)
    .verbose_message(verbose, "DGE preparation: finished")

    samples <- data.frame(
        lib.size = unname(lib.size),
        norm.factors = unname(norm.factors),
        effective.lib.size = unname(effective.lib.size),
        row.names = colnames(y.filtered),
        check.names = FALSE
    )

    list(
        counts = y.filtered,
        keep = keep,
        lib.size = lib.size,
        norm.factors = norm.factors,
        effective.lib.size = effective.lib.size,
        offset = offset,
        samples = samples
    )
}
