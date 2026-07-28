#' End-to-End BPCells edgeR LRT Pipeline
#'
#' Runs the current supported BPCells edgeR workflow: preparation, streamed or
#' dense-reference dispersion estimation when needed, and blockwise likelihood
#' ratio testing.
#'
#' @param y A genes-by-samples count matrix.
#' @param design Design matrix.
#' @param group Optional group factor for filtering and dispersion estimation.
#' @param coef,contrast Test specification passed to [bp_glm_lrt()].
#' @param filter Logical scalar. If `TRUE`, apply [bp_filter_by_expr()].
#' @param dispersion Optional scalar or gene-level dispersion vector. If `NULL`,
#'   dispersions are estimated with `dispersion.method`.
#' @param norm.method Normalization method passed to [bp_prepare_dge()].
#' @param dispersion.method Dispersion estimator to use when `dispersion` is
#'   `NULL`. `"stream"` uses [bp_estimate_disp_stream()]. `"dense"` uses
#'   [bp_estimate_disp_dense()] as a small-data oracle.
#' @param block.size Number of rows to process per block.
#' @param workers Number of row blocks to process in parallel for supported
#'   streaming steps.
#' @param max.elements Maximum number of matrix elements allowed for dense
#'   dispersion estimation and dense AveLogCPM calculation.
#' @param prior.count Prior count passed to [bp_glm_lrt()].
#' @param filter.args,dispersion.args,glm.args Additional arguments passed to
#'   [bp_prepare_dge()], the selected dispersion estimator, and
#'   [bp_glm_lrt()], respectively.
#'
#' @return A list with `prepare`, `dispersion`, and `lrt` components.
#' @export
bp_lrt <- function(y,
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
                   glm.args = list()) {
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

    lrt.args <- c(
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
        glm.args
    )
    lrt <- do.call(bp_glm_lrt, lrt.args)

    out <- list(
        prepare = prep,
        dispersion = disp.out,
        lrt = lrt
    )
    class(out) <- c("BPCellsLRT", "list")
    out
}

.dispersion_from_estimate <- function(dispersion) {
    if (!is.null(dispersion$tagwise.dispersion)) {
        return(dispersion$tagwise.dispersion)
    }
    if (!is.null(dispersion$trended.dispersion)) {
        return(dispersion$trended.dispersion)
    }
    dispersion$common.dispersion
}

.estimate_ave_log_cpm <- function(y, lib.size, dispersion, max.elements) {
    dims <- .check_count_matrix(y)
    if (prod(dims) > max.elements) {
        stop(
            "`y` has more than `max.elements` entries; dense AveLogCPM calculation is only available for the dense-reference path.",
            call. = FALSE
        )
    }
    counts <- .dense_slice(y, TRUE, TRUE)
    storage.mode(counts) <- "numeric"
    edgeR::aveLogCPM(
        counts,
        lib.size = lib.size,
        dispersion = dispersion
    )
}
