#' Streaming edgeR-Compatible Expression Filter
#'
#' Reproduces the default logic of [edgeR::filterByExpr()] while processing
#' counts in row blocks.
#'
#' @param y A genes-by-samples count matrix.
#' @param design Optional design matrix.
#' @param group Optional group factor.
#' @param lib.size Optional library sizes. If `NULL`, column sums of `y` are
#'   used. To mirror `filterByExpr.DGEList()` after normalization, pass effective
#'   library sizes explicitly.
#' @param min.count,min.total.count,large.n,min.prop Filtering parameters with
#'   the same meaning as in [edgeR::filterByExpr()].
#' @param block.size Number of rows to process per block.
#' @param workers Number of row blocks to process in parallel.
#'
#' @return A logical vector indicating retained rows.
#' @export
bp_filter_by_expr <- function(y,
                              design = NULL,
                              group = NULL,
                              lib.size = NULL,
                              min.count = 10,
                              min.total.count = 15,
                              large.n = 10,
                              min.prop = 0.7,
                              block.size = 1000,
                              workers = 1) {
    dims <- .check_count_matrix(y)
    ntags <- dims[1L]
    nlibs <- dims[2L]

    lib.size <- .check_lib_size(lib.size, nlibs)
    if (is.null(lib.size)) {
        lib.size <- bp_library_sizes(y)
    }

    if (is.null(group)) {
        if (is.null(design)) {
            message("No group or design set. Assuming all samples belong to one group.")
            min.sample.size <- nlibs
        } else {
            design <- as.matrix(design)
            min.sample.size <- 1 / max(stats::hat(design))
        }
    } else {
        group <- as.factor(group)
        group.tab <- tabulate(group)
        min.sample.size <- min(group.tab[group.tab > 0L])
    }

    if (min.sample.size > large.n) {
        min.sample.size <- large.n + (min.sample.size - large.n) * min.prop
    }

    cpm.cutoff <- min.count / stats::median(lib.size) * 1e6
    keep <- rep(FALSE, ntags)
    block.ids <- .row_block_indices(ntags, block.size)
    block.results <- .block_lapply(block.ids, workers, function(rows) {
        block <- .dense_slice(y, rows, TRUE)
        cpm <- t(t(block) / lib.size) * 1e6
        keep.cpm <- rowSums(cpm >= cpm.cutoff) >= (min.sample.size - 1e-14)
        keep.total <- rowSums(block) >= (min.total.count - 1e-14)
        list(rows = rows, keep = keep.cpm & keep.total)
    })
    for (block.result in block.results) {
        keep[block.result$rows] <- block.result$keep
    }

    names(keep) <- rownames(y)
    keep
}
