#' Blockwise edgeR GLM Quasi-Likelihood F-Test
#'
#' Fits negative-binomial GLMs and runs quasi-likelihood F-tests on
#' BPCells-backed counts by processing dense row blocks while preserving
#' edgeR's global empirical Bayes shrinkage step. Both edgeR's legacy QL path
#' (`legacy = TRUE`) and modern QL path (`legacy = FALSE`) are supported.
#'
#' @param y A genes-by-samples count matrix.
#' @param design Design matrix.
#' @param dispersion Scalar or gene-level negative-binomial dispersion vector.
#' @param offset Optional sample-level log-offset vector of length `ncol(y)`.
#'   Row-constant matrices with dimensions `nrow(y)` by `ncol(y)` are accepted
#'   and canonicalized to vectors. Gene-specific offset matrices are not yet
#'   supported.
#' @param lib.size Optional library sizes used by `edgeR::glmFit()` when
#'   `offset` is `NULL`.
#' @param coef,contrast Test specification passed to `edgeR::glmQLFTest()`.
#' @param prior.count Prior count passed to `edgeR::glmFit()`.
#' @param abundance.trend Logical scalar passed to the QL empirical Bayes trend.
#' @param ave.log.cpm Optional precomputed average log CPM vector. If supplied,
#'   it is used for the abundance trend and output table.
#' @param block.size Number of rows to process per block.
#' @param workers Number of row blocks to process in parallel.
#' @param keep.fitted Logical scalar. If `TRUE`, keep fitted values in the
#'   returned object. This may require `O(GK)` memory.
#' @param robust,winsor.tail.p Robust empirical Bayes options. `robust = TRUE`
#'   is not yet supported by the blockwise QL path.
#' @param poisson.bound Passed to `edgeR::glmQLFTest()`.
#' @param legacy Logical scalar. If `FALSE`, use edgeR's modern QL adjustment
#'   with streamed computation of the global average QL dispersion.
#' @param weights Observation weights. Currently unsupported.
#' @param keep.unit.mat Logical scalar. Unit-level QL diagnostic matrices are
#'   not yet supported.
#' @param ... Additional arguments passed to `edgeR::glmFit()`.
#'
#' @return A list containing the combined QLF table, GLM coefficients,
#'   quasi-likelihood shrinkage fields, and test metadata.
#' @export
bp_glm_qlf <- function(y,
                       design,
                       dispersion,
                       offset = NULL,
                       lib.size = NULL,
                       coef = ncol(design),
                       contrast = NULL,
                       prior.count = 0.125,
                       abundance.trend = TRUE,
                       ave.log.cpm = NULL,
                       block.size = 1000,
                       workers = 1,
                       keep.fitted = FALSE,
                       robust = FALSE,
                       winsor.tail.p = c(0.05, 0.1),
                       poisson.bound = TRUE,
                       legacy = FALSE,
                       weights = NULL,
                       keep.unit.mat = FALSE,
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
    if (isTRUE(robust)) {
        stop("`robust = TRUE` is not yet supported by blockwise QL fitting.", call. = FALSE)
    }
    legacy <- .check_logical_scalar(legacy, "legacy")
    if (!is.null(weights)) {
        stop("Observation `weights` are not yet supported by blockwise QL fitting.", call. = FALSE)
    }
    if (isTRUE(keep.unit.mat)) {
        stop("`keep.unit.mat = TRUE` is not yet supported by blockwise QL fitting.", call. = FALSE)
    }
    abundance.trend <- .check_logical_scalar(abundance.trend, "abundance.trend")
    keep.fitted <- .check_logical_scalar(keep.fitted, "keep.fitted")
    poisson.bound <- .check_logical_scalar(poisson.bound, "poisson.bound")

    lib.size <- .check_lib_size(lib.size, nlibs)
    offset <- .normalize_offset(offset, ntags, nlibs)
    if (!is.null(offset)) {
        .block_offset(offset, seq_len(min(ntags, 1L)), ntags, nlibs)
    }
    if (is.null(ave.log.cpm)) {
        ave.lib.size <- if (!is.null(offset)) exp(as.numeric(offset)) else lib.size
        ave.log.cpm <- .stream_ave_log_cpm(
            y,
            lib.size = ave.lib.size,
            dispersion = dispersion,
            block.size = block.size,
            workers = workers
        )
    } else {
        if (length(ave.log.cpm) != ntags || anyNA(ave.log.cpm)) {
            stop("`ave.log.cpm` must be a non-missing numeric vector with length nrow(y).", call. = FALSE)
        }
        ave.log.cpm <- as.numeric(ave.log.cpm)
    }

    dispersion <- .check_ql_dispersion(dispersion, ntags)
    if (!legacy && max(dispersion) > 4) {
        message("Capping the NB dispersion(s) to be no larger than 4.")
        dispersion <- pmin(dispersion, 4)
    }

    ql.stats <- if (legacy) {
        .stream_ql_legacy_stats(
            y = y,
            design = design,
            dispersion = dispersion,
            offset = offset,
            lib.size = lib.size,
            ave.log.cpm = ave.log.cpm,
            prior.count = prior.count,
            block.size = block.size,
            workers = workers,
            ...
        )
    } else {
        .stream_ql_modern_stats(
            y = y,
            design = design,
            dispersion = dispersion,
            offset = offset,
            lib.size = lib.size,
            ave.log.cpm = ave.log.cpm,
            prior.count = prior.count,
            block.size = block.size,
            workers = workers,
            coef = coef,
            contrast = contrast,
            keep.fitted = keep.fitted,
            compute.test = TRUE,
            ...
        )
    }
    ql.covariate <- if (abundance.trend) ave.log.cpm else NULL
    ql.df <- if (legacy) ql.stats$df.residual.zeros else ql.stats$df.residual.adj
    s2.fit <- limma::squeezeVar(
        ql.stats$s2,
        df = ql.df,
        covariate = ql.covariate,
        robust = robust,
        winsor.tail.p = winsor.tail.p,
        legacy = legacy
    )

    test <- if (legacy) {
        .stream_qlf_legacy_test(
            y = y,
            design = design,
            dispersion = dispersion,
            offset = offset,
            lib.size = lib.size,
            ave.log.cpm = ave.log.cpm,
            prior.count = prior.count,
            ql.stats = ql.stats,
            s2.fit = s2.fit,
            coef = coef,
            contrast = contrast,
            block.size = block.size,
            workers = workers,
            keep.fitted = keep.fitted,
            poisson.bound = poisson.bound,
            ...
        )
    } else {
        .stream_qlf_modern_test(
            ql.stats = ql.stats,
            s2.fit = s2.fit,
            keep.fitted = keep.fitted
        )
    }

    ql.fields <- if (legacy) {
        list(
            df.residual.zeros = stats::setNames(ql.stats$df.residual.zeros, rownames(y))
        )
    } else {
        list(
            df.residual.adj = stats::setNames(ql.stats$df.residual.adj, rownames(y)),
            deviance.adj = stats::setNames(ql.stats$deviance.adj, rownames(y)),
            average.ql.dispersion = ql.stats$average.ql.dispersion
        )
    }
    out <- c(
        test,
        list(
            design = design,
            dispersion = if (length(dispersion) == 1L) rep(dispersion, ntags) else dispersion,
            prior.count = prior.count,
            AveLogCPM = stats::setNames(ave.log.cpm, rownames(y)),
            df.prior = s2.fit$df.prior,
            s2.prior = s2.fit$var.prior,
            s2.post = stats::setNames(s2.fit$var.post, rownames(y))
        ),
        ql.fields
    )
    class(out) <- c("BPCellsDGEQLF", "list")
    out
}

.stream_ql_legacy_stats <- function(y,
                                    design,
                                    dispersion,
                                    offset,
                                    lib.size,
                                    ave.log.cpm,
                                    prior.count,
                                    block.size,
                                    workers = 1,
                                    ...) {
    dims <- .check_count_matrix(y)
    ntags <- dims[1L]
    nlibs <- dims[2L]
    s2 <- rep_len(NA_real_, ntags)
    df.residual.zeros <- rep_len(NA_real_, ntags)
    residDF <- .edgeR_internal(".residDF")

    block.ids <- .row_block_indices(ntags, block.size)
    block.results <- .block_lapply(block.ids, workers, function(rows) {
        block <- .dense_slice(y, rows, TRUE)
        storage.mode(block) <- "numeric"
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
        zerofit <- (fit$fitted.values < 1e-4) & (fit$counts < 1e-4)
        block.df <- residDF(zerofit, fit$design)
        block.s2 <- fit$deviance / block.df
        block.s2[block.df == 0L] <- 0

        list(rows = rows, s2 = block.s2, df.residual.zeros = block.df)
    })
    for (block.result in block.results) {
        s2[block.result$rows] <- block.result$s2
        df.residual.zeros[block.result$rows] <- block.result$df.residual.zeros
    }

    list(
        s2 = s2,
        df.residual.zeros = df.residual.zeros
    )
}

.stream_ql_modern_stats <- function(y,
                                    design,
                                    dispersion,
                                    offset,
                                    lib.size,
                                    ave.log.cpm,
                                    prior.count,
                                    block.size,
                                    workers = 1,
                                    coef = ncol(design),
                                    contrast = NULL,
                                    keep.fitted = FALSE,
                                    keep.coefficients = FALSE,
                                    compute.test = FALSE,
                                    ...) {
    dims <- .check_count_matrix(y)
    ntags <- dims[1L]
    nlibs <- dims[2L]
    average.ql.dispersion <- .stream_ql_modern_average_dispersion(
        y = y,
        design = design,
        dispersion = dispersion,
        offset = offset,
        lib.size = lib.size,
        ave.log.cpm = ave.log.cpm,
        prior.count = prior.count,
        block.size = block.size,
        workers = workers,
        ...
    )

    block.ids <- .row_block_indices(ntags, block.size)
    s2 <- rep_len(NA_real_, ntags)
    df.residual.adj <- rep_len(NA_real_, ntags)
    deviance.adj <- rep_len(NA_real_, ntags)
    computeAdjVec <- .edgeR_internal(".cxx_compute_adj_vec")
    compressDispersions <- .edgeR_internal(".compressDispersions")

    block.results <- .block_lapply(block.ids, workers, function(rows) {
        block <- .dense_slice(y, rows, TRUE)
        storage.mode(block) <- "numeric"
        rownames(block) <- rownames(y)[rows]
        colnames(block) <- colnames(y)

        block.offset <- .block_offset(offset, rows, ntags, nlibs)
        block.dispersion <- .block_dispersion(dispersion, rows, ntags)
        fit <- edgeR::glmFit(
            block,
            design = design,
            dispersion = block.dispersion / average.ql.dispersion,
            offset = block.offset,
            lib.size = lib.size,
            prior.count = prior.count,
            ...
        )
        dispersion.mat <- compressDispersions(block, block.dispersion)
        weights.mat <- edgeR::makeCompressedMatrix(1, dim(block), byrow = TRUE)
        adj <- .Call(
            computeAdjVec,
            block,
            fit$fitted.values,
            design,
            dispersion.mat,
            average.ql.dispersion,
            weights.mat
        )

        test <- NULL
        if (isTRUE(compute.test)) {
            # The likelihood-ratio statistic and null fit do not depend on the
            # global QL squeeze. Compute them while the scaled-dispersion
            # alternative fit is already available, rather than refitting the
            # same alternative model in a later block pass.
            fit$AveLogCPM <- ave.log.cpm[rows]
            lrt <- edgeR::glmLRT(fit, coef = coef, contrast = contrast)
            # glmLRT returns the complete fitted object, including its dense
            # fitted-value matrix. Keep only the per-gene fields needed after
            # the global QL squeeze so completed blocks do not accumulate
            # O(GN) temporary storage in the parent process.
            test <- list(
                table = lrt$table,
                coefficients = lrt$coefficients,
                unshrunk.coefficients = lrt$unshrunk.coefficients,
                deviance = lrt$deviance,
                df.residual = lrt$df.residual,
                df.test = lrt$df.test,
                comparison = lrt$comparison
            )
        }

        list(
            rows = rows,
            s2 = adj$s2,
            df.residual.adj = adj$df,
            deviance.adj = adj$deviance,
            test = test,
            # Modern Treat uses this exact scaled-dispersion fit. Retaining
            # only its small coefficient matrices avoids repeating the full
            # genes-by-samples GLM optimization in the Treat pass.
            coefficients = if (isTRUE(keep.coefficients)) fit$coefficients else NULL,
            unshrunk.coefficients = if (isTRUE(keep.coefficients)) {
                fit$unshrunk.coefficients
            } else {
                NULL
            },
            fitted.values = if (isTRUE(compute.test) && isTRUE(keep.fitted)) {
                fit$fitted.values
            } else {
                NULL
            }
        )
    })

    for (block.result in block.results) {
        rows <- block.result$rows
        s2[rows] <- block.result$s2
        df.residual.adj[rows] <- block.result$df.residual.adj
        deviance.adj[rows] <- block.result$deviance.adj
    }

    out <- list(
        s2 = s2,
        df.residual.adj = df.residual.adj,
        deviance.adj = deviance.adj,
        average.ql.dispersion = average.ql.dispersion
    )
    if (isTRUE(keep.coefficients)) {
        out$coefficients <- do.call(
            rbind,
            lapply(block.results, `[[`, "coefficients")
        )
        out$unshrunk.coefficients <- do.call(
            rbind,
            lapply(block.results, `[[`, "unshrunk.coefficients")
        )
        rownames(out$coefficients) <- rownames(y)
        rownames(out$unshrunk.coefficients) <- rownames(y)
    }
    if (isTRUE(compute.test)) {
        test.results <- lapply(block.results, `[[`, "test")
        out$test <- list(
            table = do.call(rbind, lapply(test.results, `[[`, "table")),
            coefficients = do.call(rbind, lapply(test.results, `[[`, "coefficients")),
            unshrunk.coefficients = do.call(
                rbind,
                lapply(test.results, `[[`, "unshrunk.coefficients")
            ),
            deviance = unlist(lapply(test.results, `[[`, "deviance"), use.names = FALSE),
            df.residual = unlist(
                lapply(test.results, `[[`, "df.residual"),
                use.names = FALSE
            ),
            df.test = unlist(
                lapply(test.results, function(x) rep_len(x$df.test, nrow(x$table))),
                use.names = FALSE
            ),
            comparison = test.results[[1L]]$comparison
        )
        if (isTRUE(keep.fitted)) {
            out$test$fitted.values <- do.call(
                rbind,
                lapply(block.results, `[[`, "fitted.values")
            )
        }
        rownames(out$test$coefficients) <- rownames(y)
        rownames(out$test$unshrunk.coefficients) <- rownames(y)
        rownames(out$test$table) <- rownames(y)
        names(out$test$deviance) <- rownames(y)
        names(out$test$df.residual) <- rownames(y)
        names(out$test$df.test) <- rownames(y)
        if (isTRUE(keep.fitted)) {
            rownames(out$test$fitted.values) <- rownames(y)
            colnames(out$test$fitted.values) <- colnames(y)
        }
    }
    out
}

.stream_ql_modern_average_dispersion <- function(y,
                                                 design,
                                                 dispersion,
                                                 offset,
                                                 lib.size,
                                                 ave.log.cpm,
                                                 prior.count,
                                                 block.size,
                                                 workers = 1,
                                                 ...) {
    dims <- .check_count_matrix(y)
    ntags <- dims[1L]
    nlibs <- dims[2L]
    block.ids <- .row_block_indices(ntags, block.size)
    computeAdjVec <- .edgeR_internal(".cxx_compute_adj_vec")
    compressDispersions <- .edgeR_internal(".compressDispersions")

    # The two average-QL iterations use the same NB fit. Retain only its small
    # genes-by-coefficients matrix, then reconstruct fitted means during the
    # second count scan instead of repeating the GLM optimization. A zero prior
    # count also avoids computing reported/shrunken coefficients for this
    # intermediate fit; it does not alter its fitted means or deviance.
    first.results <- .block_lapply(seq_along(block.ids), workers, function(b) {
        rows <- block.ids[[b]]
        block <- .dense_slice(y, rows, TRUE)
        storage.mode(block) <- "numeric"
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
            prior.count = 0,
            ...
        )
        dispersion.mat <- compressDispersions(block, block.dispersion)
        weights.mat <- edgeR::makeCompressedMatrix(1, dim(block), byrow = TRUE)
        adj <- .Call(
            computeAdjVec,
            block,
            fit$fitted.values,
            design,
            dispersion.mat,
            1,
            weights.mat
        )

        list(
            rows = rows,
            s2 = adj$s2,
            df = adj$df,
            coefficients = fit$coefficients
        )
    })
    first <- .combine_ql_average_results(first.results, ntags)
    prior <- .ql_modern_compute_prior(ave.log.cpm, first$s2, first$df)

    second.results <- .block_lapply(seq_along(block.ids), workers, function(b) {
        rows <- block.ids[[b]]
        block <- .dense_slice(y, rows, TRUE)
        storage.mode(block) <- "numeric"
        rownames(block) <- rownames(y)[rows]
        colnames(block) <- colnames(y)

        block.offset <- .block_offset(offset, rows, ntags, nlibs)
        block.dispersion <- .block_dispersion(dispersion, rows, ntags)
        fitted.values <- .ql_fitted_from_coefficients(
            coefficients = first.results[[b]]$coefficients,
            design = design,
            offset = block.offset,
            lib.size = lib.size
        )
        dispersion.mat <- compressDispersions(block, block.dispersion)
        weights.mat <- edgeR::makeCompressedMatrix(1, dim(block), byrow = TRUE)
        adj <- .Call(
            computeAdjVec,
            block,
            fitted.values,
            design,
            dispersion.mat,
            prior,
            weights.mat
        )

        list(rows = rows, s2 = adj$s2, df = adj$df)
    })
    second <- .combine_ql_average_results(second.results, ntags)
    .ql_modern_compute_prior(ave.log.cpm, second$s2, second$df)
}

.combine_ql_average_results <- function(block.results, ntags) {
    s2 <- rep_len(NA_real_, ntags)
    df <- rep_len(NA_real_, ntags)
    for (block.result in block.results) {
        s2[block.result$rows] <- block.result$s2
        df[block.result$rows] <- block.result$df
    }
    list(s2 = s2, df = df)
}

.ql_fitted_from_coefficients <- function(coefficients, design, offset, lib.size) {
    finite.rows <- apply(coefficients, 1L, function(x) all(is.finite(x)))
    fitted.values <- matrix(0, nrow(coefficients), nrow(design))
    if (!any(finite.rows)) {
        return(fitted.values)
    }
    finite.coefficients <- coefficients[finite.rows, , drop = FALSE]
    design.group <- edgeR::designAsFactor(design)
    nprofiles <- nlevels(design.group)
    if (nprofiles < nrow(design)) {
        # One-vs-rest and other group designs repeat a small number of design
        # rows. Form each distinct linear predictor once, then expand it to
        # samples. General continuous-covariate designs fall through to the
        # ordinary matrix product below.
        first <- match(seq_len(nprofiles), as.integer(design.group))
        finite.fitted <- finite.coefficients %*% t(design[first, , drop = FALSE])
        finite.fitted <- finite.fitted[, as.integer(design.group), drop = FALSE]
    } else {
        finite.fitted <- finite.coefficients %*% t(design)
    }
    sample.offset <- if (is.null(offset)) log(lib.size) else as.numeric(offset)
    finite.fitted <- sweep(finite.fitted, 2L, sample.offset, `+`)
    fitted.values[finite.rows, ] <- exp(finite.fitted)
    fitted.values
}

.ql_modern_compute_prior <- function(ave.log.cpm, s2, df) {
    keep <- df > 1e-8
    if (!any(keep)) {
        stop("Cannot compute modern QL average dispersion; all adjusted residual degrees of freedom are zero.", call. = FALSE)
    }

    x <- ave.log.cpm[keep]
    y <- s2[keep]^(1 / 4)
    ord <- order(x)
    x <- x[ord]
    y <- y[ord]
    delta <- 0.01 * (x[length(x)] - x[1L])
    trend <- stats::lowess(x, y, f = 0.5, iter = 3L, delta = delta)$y
    prior <- as.numeric(stats::quantile(trend, probs = 0.9, type = 7, names = FALSE))
    if (prior < 1) {
        prior <- 1
    }
    prior^4
}

.stream_qlf_legacy_test <- function(y,
                                    design,
                                    dispersion,
                                    offset,
                                    lib.size,
                                    ave.log.cpm,
                                    prior.count,
                                    ql.stats,
                                    s2.fit,
                                    coef,
                                    contrast,
                                    block.size,
                                    workers = 1,
                                    keep.fitted,
                                    poisson.bound,
                                    ...) {
    dims <- .check_count_matrix(y)
    ntags <- dims[1L]
    nlibs <- dims[2L]
    block.ids <- .row_block_indices(ntags, block.size)
    coef.list <- vector("list", length(block.ids))
    unshrunk.list <- vector("list", length(block.ids))
    fitted.list <- if (isTRUE(keep.fitted)) vector("list", length(block.ids)) else NULL
    deviance <- rep_len(NA_real_, ntags)
    df.residual <- rep_len(NA_real_, ntags)
    df.total <- rep_len(NA_real_, ntags)
    df.test.out <- rep_len(NA_real_, ntags)
    table.list <- vector("list", length(block.ids))
    comparison <- NULL
    df.residual.total <- ntags * (nlibs - ncol(design))
    isBelowPoissonBound <- .edgeR_internal(".isBelowPoissonBound")

    block.results <- .block_lapply(seq_along(block.ids), workers, function(b) {
        rows <- block.ids[[b]]
        block <- .dense_slice(y, rows, TRUE)
        storage.mode(block) <- "numeric"
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
        fit$AveLogCPM <- ave.log.cpm[rows]
        fit$df.residual.zeros <- ql.stats$df.residual.zeros[rows]
        fit$df.prior <- .subset_prior_field(s2.fit$df.prior, rows, ntags)
        fit$s2.post <- s2.fit$var.post[rows]
        fit$s2.prior <- .subset_prior_field(s2.fit$var.prior, rows, ntags)
        fit$top.proportion <- NULL

        qlf <- edgeR::glmLRT(fit, coef = coef, contrast = contrast)
        F.stat <- qlf$table$LR / qlf$df.test / fit$s2.post
        block.df.total <- fit$df.prior + fit$df.residual.zeros
        block.df.total <- pmin(block.df.total, df.residual.total)
        F.pvalue <- stats::pf(
            F.stat,
            df1 = qlf$df.test,
            df2 = block.df.total,
            lower.tail = FALSE,
            log.p = FALSE
        )
        if (poisson.bound) {
            below.bound <- isBelowPoissonBound(fit)
            if (any(below.bound)) {
                pois.fit <- fit[below.bound, ]
                pois.fit <- edgeR::glmFit(
                    pois.fit$counts,
                    design = pois.fit$design,
                    offset = pois.fit$offset,
                    weights = pois.fit$weights,
                    start = pois.fit$unshrunk.coefficients,
                    dispersion = 0
                )
                pois.res <- edgeR::glmLRT(pois.fit, coef = coef, contrast = contrast)
                F.pvalue[below.bound] <- pmax(F.pvalue[below.bound], pois.res$table$PValue)
            }
        }
        qlf$table$LR <- qlf$table$PValue <- NULL
        qlf$table$F <- F.stat
        qlf$table$PValue <- F.pvalue
        qlf$df.total <- block.df.total

        list(
            rows = rows,
            coefficients = qlf$coefficients,
            unshrunk.coefficients = qlf$unshrunk.coefficients,
            fitted.values = if (isTRUE(keep.fitted)) fit$fitted.values else NULL,
            deviance = qlf$deviance,
            df.residual = qlf$df.residual,
            df.total = qlf$df.total,
            df.test = rep_len(qlf$df.test, length(rows)),
            table = qlf$table,
            comparison = qlf$comparison
        )
    })

    coef.list <- lapply(block.results, `[[`, "coefficients")
    unshrunk.list <- lapply(block.results, `[[`, "unshrunk.coefficients")
    fitted.list <- if (isTRUE(keep.fitted)) lapply(block.results, `[[`, "fitted.values") else NULL
    table.list <- lapply(block.results, `[[`, "table")
    for (block.result in block.results) {
        rows <- block.result$rows
        deviance[rows] <- block.result$deviance
        df.residual[rows] <- block.result$df.residual
        df.total[rows] <- block.result$df.total
        df.test.out[rows] <- block.result$df.test
    }
    comparison <- block.results[[1L]]$comparison

    coefficients <- do.call(rbind, coef.list)
    unshrunk.coefficients <- do.call(rbind, unshrunk.list)
    table <- do.call(rbind, table.list)
    rownames(coefficients) <- rownames(y)
    rownames(unshrunk.coefficients) <- rownames(y)
    rownames(table) <- rownames(y)
    names(deviance) <- rownames(y)
    names(df.residual) <- rownames(y)
    names(df.total) <- rownames(y)
    names(df.test.out) <- rownames(y)

    out <- list(
        table = table,
        coefficients = coefficients,
        unshrunk.coefficients = unshrunk.coefficients,
        deviance = deviance,
        df.residual = df.residual,
        comparison = comparison,
        df.test = df.test.out,
        df.total = df.total
    )
    if (isTRUE(keep.fitted)) {
        fitted.values <- do.call(rbind, fitted.list)
        rownames(fitted.values) <- rownames(y)
        colnames(fitted.values) <- colnames(y)
        out$fitted.values <- fitted.values
    }
    out
}

.stream_qlf_modern_test <- function(ql.stats, s2.fit, keep.fitted) {
    test <- ql.stats$test
    if (is.null(test)) {
        stop("Modern QL test statistics were not computed with the fit.", call. = FALSE)
    }

    F.stat <- test$table$LR / test$df.test / s2.fit$var.post
    df.residual.total <- sum(test$df.residual)
    df.total <- s2.fit$df.prior + ql.stats$df.residual.adj
    df.total <- pmin(df.total, df.residual.total)
    F.pvalue <- stats::pf(
        F.stat,
        df1 = test$df.test,
        df2 = df.total,
        lower.tail = FALSE,
        log.p = FALSE
    )
    test$table$LR <- test$table$PValue <- NULL
    test$table$F <- F.stat
    test$table$PValue <- F.pvalue
    test$df.total <- df.total
    names(test$df.total) <- names(test$df.residual)
    if (!isTRUE(keep.fitted)) {
        test$fitted.values <- NULL
    }
    test
}

.check_ql_dispersion <- function(dispersion, nrows) {
    if (is.null(dispersion)) {
        stop("`dispersion` must be supplied for blockwise QL fitting.", call. = FALSE)
    }
    if (anyNA(dispersion) || !is.numeric(dispersion)) {
        stop("`dispersion` must be a non-missing numeric scalar or vector.", call. = FALSE)
    }
    if (length(dispersion) == 1L) {
        return(as.numeric(dispersion))
    }
    if (length(dispersion) != nrows) {
        stop("`dispersion` must be scalar or have length nrow(y).", call. = FALSE)
    }
    as.numeric(dispersion)
}

.subset_prior_field <- function(x, rows, nrows) {
    if (length(x) == 1L) {
        return(x)
    }
    if (length(x) == nrows) {
        return(x[rows])
    }
    x
}
