#' End-to-End BPCells edgeR QLF Pipeline
#'
#' Runs the supported BPCells edgeR quasi-likelihood workflow: preparation,
#' streamed or dense-reference dispersion estimation when needed, and blockwise
#' quasi-likelihood F-testing.
#'
#' @inheritParams bp_lrt
#' @param qlf.args Additional arguments passed to [bp_glm_qlf()].
#'
#' @return A list with `prepare`, `dispersion`, and `qlf` components.
#' @export
bp_qlf <- function(y,
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
                   max.elements = 1e7,
                   prior.count = 0.125,
                   filter.args = list(),
                   dispersion.args = list(),
                   qlf.args = list()) {
    design <- as.matrix(design)
    dispersion.method <- match.arg(dispersion.method)
    prepare.args <- c(
        list(
            y = y,
            design = design,
            group = group,
            filter = filter,
            norm.method = norm.method,
            block.size = block.size,
            workers = workers
        ),
        filter.args
    )
    prep <- do.call(bp_prepare_dge, prepare.args)

    disp.out <- NULL
    ave.log.cpm <- NULL
    if (is.null(dispersion)) {
        estimate.defaults <- list(
            y = prep$counts,
            design = design,
            group = group,
            lib.size = prep$effective.lib.size,
            offset = prep$offset
        )
        if (dispersion.method == "stream") {
            estimate.defaults$block.size <- block.size
            estimate.defaults$workers <- workers
            estimate.args <- utils::modifyList(estimate.defaults, dispersion.args)
            disp.out <- do.call(bp_estimate_disp_stream, estimate.args)
            ave.log.cpm <- .stream_ave_log_cpm(
                prep$counts,
                lib.size = estimate.args$lib.size,
                dispersion = disp.out$common.dispersion,
                block.size = estimate.args$block.size,
                workers = estimate.args$workers
            )
        } else {
            estimate.defaults$max.elements <- max.elements
            estimate.args <- utils::modifyList(estimate.defaults, dispersion.args)
            disp.out <- do.call(bp_estimate_disp_dense, estimate.args)
            ave.log.cpm <- .estimate_ave_log_cpm(
                prep$counts,
                lib.size = estimate.args$lib.size,
                dispersion = disp.out$common.dispersion,
                max.elements = max.elements
            )
        }
        dispersion <- .dispersion_from_estimate(disp.out)
    }

    qlf.call.args <- c(
        list(
            y = prep$counts,
            design = design,
            dispersion = dispersion,
            offset = prep$offset,
            coef = coef,
            contrast = contrast,
            prior.count = prior.count,
            block.size = block.size,
            workers = workers,
            ave.log.cpm = ave.log.cpm
        ),
        qlf.args
    )
    qlf <- do.call(bp_glm_qlf, qlf.call.args)

    out <- list(
        prepare = prep,
        dispersion = disp.out,
        qlf = qlf
    )
    class(out) <- c("BPCellsQLF", "list")
    out
}
