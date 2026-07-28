#' Blockwise edgeR GLM Likelihood Ratio Test
#'
#' Fits negative-binomial GLMs and runs likelihood ratio tests on BPCells-backed
#' counts by processing dense row blocks with public edgeR wrappers.
#'
#' @param y A genes-by-samples count matrix.
#' @param design Design matrix.
#' @param dispersion Scalar or gene-level dispersion vector. Dispersion
#'   estimation is intentionally not performed by this function.
#' @param offset Optional sample-level log-offset vector of length `ncol(y)`.
#'   Row-constant matrices with dimensions `nrow(y)` by `ncol(y)` are accepted
#'   and canonicalized to vectors. Gene-specific offset matrices are not yet
#'   supported.
#' @param lib.size Optional library sizes used by `edgeR::glmFit()` when
#'   `offset` is `NULL`.
#' @param coef,contrast Test specification passed to `edgeR::glmLRT()`.
#' @param prior.count Prior count passed to `edgeR::glmFit()`.
#' @param block.size Number of rows to process per block.
#' @param workers Number of row blocks to process in parallel.
#' @param keep.fitted Logical scalar. If `TRUE`, keep fitted values in the
#'   returned object. This may require `O(GK)` memory.
#' @param ave.log.cpm Optional precomputed average log CPM vector. If supplied,
#'   it replaces the blockwise `logCPM` column from `edgeR::glmLRT()`.
#' @param ... Additional arguments passed to `edgeR::glmFit()`.
#'
#' @return A list containing the combined LRT table, GLM coefficients,
#'   deviance, residual degrees of freedom, dispersions, and test metadata.
#' @export
bp_glm_lrt <- function(y,
                       design,
                       dispersion,
                       offset = NULL,
                       lib.size = NULL,
                       coef = ncol(design),
                       contrast = NULL,
                       prior.count = 0.125,
                       block.size = 1000,
                       workers = 1,
                       keep.fitted = FALSE,
                       ave.log.cpm = NULL,
                       ...) {
    dims <- .check_count_matrix(y)
    ntags <- dims[1L]
    nlibs <- dims[2L]
    design <- as.matrix(design)
    if (nrow(design) != nlibs) {
        stop("`design` must have one row per sample.", call. = FALSE)
    }

    if (is.null(offset) && is.null(lib.size)) {
        stop("Supply `offset` or `lib.size`; block-specific column sums are not valid global library sizes.", call. = FALSE)
    }
    lib.size <- .check_lib_size(lib.size, nlibs)
    offset <- .normalize_offset(offset, ntags, nlibs)
    if (!is.null(offset)) {
        .block_offset(offset, seq_len(min(ntags, 1L)), ntags, nlibs)
    }
    if (!is.null(ave.log.cpm)) {
        if (length(ave.log.cpm) != ntags || anyNA(ave.log.cpm)) {
            stop("`ave.log.cpm` must be a non-missing numeric vector with length nrow(y).", call. = FALSE)
        }
        ave.log.cpm <- as.numeric(ave.log.cpm)
    }

    block.ids <- .row_block_indices(ntags, block.size)
    block.results <- .block_lapply(seq_along(block.ids), workers, function(b) {
        rows <- block.ids[[b]]
        block <- .dense_slice(y, rows, TRUE)
        rownames(block) <- rownames(y)[rows]
        colnames(block) <- colnames(y)

        block.offset <- .block_offset(offset, rows, ntags, nlibs)
        block.dispersion <- .block_dispersion(dispersion, rows, ntags)

        fit <- edgeR::glmFit(
            block,
            design = design,
            dispersion = block.dispersion,
            offset = block.offset,
            lib.size = lib.size,
            prior.count = prior.count,
            ...
        )
        lrt <- edgeR::glmLRT(fit, coef = coef, contrast = contrast)

        list(
            rows = rows,
            coefficients = lrt$coefficients,
            unshrunk.coefficients = lrt$unshrunk.coefficients,
            fitted.values = if (isTRUE(keep.fitted)) lrt$fitted.values else NULL,
            deviance = lrt$deviance,
            df.residual = lrt$df.residual,
            dispersion = lrt$dispersion,
            df.test = rep_len(lrt$df.test, length(rows)),
            table = lrt$table,
            comparison = lrt$comparison
        )
    })

    coef.list <- lapply(block.results, `[[`, "coefficients")
    unshrunk.list <- lapply(block.results, `[[`, "unshrunk.coefficients")
    fitted.list <- if (isTRUE(keep.fitted)) lapply(block.results, `[[`, "fitted.values") else NULL
    table.list <- lapply(block.results, `[[`, "table")
    deviance <- rep_len(NA_real_, ntags)
    df.residual <- rep_len(NA_real_, ntags)
    dispersion.out <- rep_len(NA_real_, ntags)
    df.test.out <- rep_len(NA_real_, ntags)
    for (block.result in block.results) {
        rows <- block.result$rows
        deviance[rows] <- block.result$deviance
        df.residual[rows] <- block.result$df.residual
        dispersion.out[rows] <- block.result$dispersion
        df.test.out[rows] <- block.result$df.test
    }
    comparison <- block.results[[1L]]$comparison

    coefficients <- do.call(rbind, coef.list)
    unshrunk.coefficients <- do.call(rbind, unshrunk.list)
    table <- do.call(rbind, table.list)
    if (!is.null(ave.log.cpm) && "logCPM" %in% colnames(table)) {
        table$logCPM <- ave.log.cpm
    }
    rownames(coefficients) <- rownames(y)
    rownames(unshrunk.coefficients) <- rownames(y)
    rownames(table) <- rownames(y)
    names(deviance) <- rownames(y)
    names(df.residual) <- rownames(y)
    names(dispersion.out) <- rownames(y)
    names(df.test.out) <- rownames(y)

    out <- list(
        table = table,
        coefficients = coefficients,
        unshrunk.coefficients = unshrunk.coefficients,
        deviance = deviance,
        df.residual = df.residual,
        design = design,
        dispersion = dispersion.out,
        prior.count = prior.count,
        comparison = comparison,
        df.test = df.test.out
    )
    if (!is.null(ave.log.cpm)) {
        names(ave.log.cpm) <- rownames(y)
        out$AveLogCPM <- ave.log.cpm
    }
    if (isTRUE(keep.fitted)) {
        fitted.values <- do.call(rbind, fitted.list)
        rownames(fitted.values) <- rownames(y)
        colnames(fitted.values) <- colnames(y)
        out$fitted.values <- fitted.values
    }

    class(out) <- c("BPCellsDGELRT", "list")
    out
}
