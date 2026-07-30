#' Blockwise edgeR GLM Treat Test
#'
#' Fits quasi-likelihood GLMs and runs Treat tests on BPCells-backed counts by
#' processing dense row blocks while preserving edgeR's global empirical Bayes
#' shrinkage step. Both edgeR's legacy QL path (`legacy = TRUE`) and modern QL
#' path (`legacy = FALSE`) are supported.
#'
#' @inheritParams bp_glm_qlf
#' @param lfc Log2 fold-change threshold passed to [edgeR::glmTreat()].
#' @param null Null hypothesis passed to [edgeR::glmTreat()].
#'
#' @return A list containing the combined Treat table, GLM coefficients,
#'   quasi-likelihood shrinkage fields, and test metadata.
#' @export
bp_glm_treat <- function(y,
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
                         robust = FALSE,
                         winsor.tail.p = c(0.05, 0.1),
                         lfc = log2(1.2),
                         null = c("interval", "worst.case"),
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
        stop("`robust = TRUE` is not yet supported by blockwise Treat fitting.", call. = FALSE)
    }
    legacy <- .check_logical_scalar(legacy, "legacy")
    if (!is.null(weights)) {
        stop("Observation `weights` are not yet supported by blockwise Treat fitting.", call. = FALSE)
    }
    if (isTRUE(keep.unit.mat)) {
        stop("`keep.unit.mat = TRUE` is not yet supported by blockwise Treat fitting.", call. = FALSE)
    }
    abundance.trend <- .check_logical_scalar(abundance.trend, "abundance.trend")
    null <- match.arg(null)
    if (length(lfc) != 1L || is.na(lfc) || lfc < 0) {
        stop("`lfc` must be a non-negative numeric scalar.", call. = FALSE)
    }

    if (lfc == 0) {
        return(bp_glm_qlf(
            y = y,
            design = design,
            dispersion = dispersion,
            offset = offset,
            lib.size = lib.size,
            coef = coef,
            contrast = contrast,
            prior.count = prior.count,
            abundance.trend = abundance.trend,
            ave.log.cpm = ave.log.cpm,
            block.size = block.size,
            workers = workers,
            robust = robust,
            winsor.tail.p = winsor.tail.p,
            legacy = legacy,
            weights = weights,
            keep.unit.mat = keep.unit.mat,
            ...
        ))
    }

    lib.size <- .check_lib_size(lib.size, nlibs)
    offset <- .normalize_offset(offset, ntags, nlibs)
    if (is.null(offset)) {
        offset <- log(lib.size)
    }
    offset <- .check_sample_offset(offset, nlibs)
    if (is.null(ave.log.cpm)) {
        ave.log.cpm <- .stream_ave_log_cpm(
            y,
            lib.size = exp(offset),
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
            lib.size = NULL,
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
            lib.size = NULL,
            ave.log.cpm = ave.log.cpm,
            prior.count = prior.count,
            block.size = block.size,
            workers = workers,
            keep.coefficients = TRUE,
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
        .stream_treat_legacy_test(
            y = y,
            design = design,
            dispersion = dispersion,
            offset = offset,
            ave.log.cpm = ave.log.cpm,
            prior.count = prior.count,
            ql.stats = ql.stats,
            s2.fit = s2.fit,
            coef = coef,
            contrast = contrast,
            lfc = lfc,
            null = null,
            block.size = block.size,
            workers = workers,
            ...
        )
    } else {
        .stream_treat_modern_test(
            y = y,
            design = design,
            dispersion = dispersion,
            offset = offset,
            ave.log.cpm = ave.log.cpm,
            prior.count = prior.count,
            ql.stats = ql.stats,
            s2.fit = s2.fit,
            coef = coef,
            contrast = contrast,
            lfc = lfc,
            null = null,
            block.size = block.size,
            workers = workers,
            ...
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
            s2.post = stats::setNames(s2.fit$var.post, rownames(y)),
            lfc = lfc,
            null = null
        ),
        ql.fields
    )
    class(out) <- c("BPCellsDGETreat", "list")
    out
}

.stream_treat_legacy_test <- function(y,
                                      design,
                                      dispersion,
                                      offset,
                                      ave.log.cpm,
                                      prior.count,
                                      ql.stats,
                                      s2.fit,
                                      coef,
                                      contrast,
                                      lfc,
                                      null,
                                      block.size,
                                      workers = 1,
                                      ...) {
    dims <- .check_count_matrix(y)
    ntags <- dims[1L]
    nlibs <- dims[2L]
    block.ids <- .row_block_indices(ntags, block.size)
    treat.info <- .treat_design_info(design, coef, contrast)
    design.treat <- treat.info$design
    design0 <- design.treat[, -treat.info$coef, drop = FALSE]
    offset.adj <- lfc * log(2) * design.treat[, treat.info$coef]
    df.total.cap <- ntags * (nlibs - ncol(design))
    zscoreT <- utils::getFromNamespace("zscoreT", "limma")
    integratepnorm <- .edgeR_internal(".integratepnorm")
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
            prior.count = prior.count,
            ...
        )
        fit$AveLogCPM <- ave.log.cpm[rows]
        fit$df.residual.zeros <- ql.stats$df.residual.zeros[rows]
        fit$df.prior <- .subset_prior_field(s2.fit$df.prior, rows, ntags)
        fit$s2.post <- s2.fit$var.post[rows]
        fit$s2.prior <- .subset_prior_field(s2.fit$var.prior, rows, ntags)
        fit$top.proportion <- NULL

        logfc <- .treat_logfc(fit, treat.info)
        fit0 <- edgeR::glmFit(
            block,
            design = design0,
            offset = block.offset + offset.adj,
            dispersion = block.dispersion,
            prior.count = 0,
            ...
        )
        fit1 <- edgeR::glmFit(
            block,
            design = design.treat,
            offset = block.offset + offset.adj,
            dispersion = block.dispersion,
            prior.count = 0,
            ...
        )
        z.left <- sqrt(pmax(0, fit0$deviance - fit1$deviance))

        fit0 <- edgeR::glmFit(
            block,
            design = design0,
            offset = block.offset - offset.adj,
            dispersion = block.dispersion,
            prior.count = 0,
            ...
        )
        fit1 <- edgeR::glmFit(
            block,
            design = design.treat,
            offset = block.offset - offset.adj,
            dispersion = block.dispersion,
            prior.count = 0,
            ...
        )
        z.right <- sqrt(pmax(0, fit0$deviance - fit1$deviance))

        swap <- z.left > z.right
        if (any(swap)) {
            tmp <- z.left[swap]
            z.left[swap] <- z.right[swap]
            z.right[swap] <- tmp
        }

        block.df.total <- fit$df.prior + fit$df.residual.zeros
        block.df.total <- pmin(block.df.total, df.total.cap)
        z.left <- zscoreT(z.left / sqrt(fit$s2.post), df = block.df.total)
        z.right <- zscoreT(z.right / sqrt(fit$s2.post), df = block.df.total)

        within <- abs(logfc$unshrunk.logFC) <= lfc
        sgn <- 2 * within - 1
        z.left <- z.left * sgn

        p.value <- .treat_p_value(z.left, z.right, null, integratepnorm)

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
            pois.res <- edgeR::glmTreat(
                pois.fit,
                coef = treat.info$coef,
                contrast = treat.info$contrast,
                lfc = lfc,
                null = null
            )
            p.value[below.bound] <- pmax(p.value[below.bound], pois.res$table$PValue)
        }

        tab <- data.frame(
            logCPM = fit$AveLogCPM,
            PValue = p.value,
            row.names = rownames(block)
        )
        if (!is.null(logfc$logFC)) {
            tab <- cbind(
                logFC = logfc$logFC,
                unshrunk.logFC = logfc$unshrunk.logFC,
                tab
            )
        } else {
            tab <- cbind(logFC = logfc$unshrunk.logFC, tab)
        }

        list(
            rows = rows,
            coefficients = fit$coefficients,
            unshrunk.coefficients = fit$unshrunk.coefficients,
            deviance = fit$deviance,
            df.residual = fit$df.residual,
            df.total = block.df.total,
            table = tab
        )
    })

    coefficients <- do.call(rbind, lapply(block.results, `[[`, "coefficients"))
    unshrunk.coefficients <- do.call(rbind, lapply(block.results, `[[`, "unshrunk.coefficients"))
    table <- do.call(rbind, lapply(block.results, `[[`, "table"))
    deviance <- rep_len(NA_real_, ntags)
    df.residual <- rep_len(NA_real_, ntags)
    df.total <- rep_len(NA_real_, ntags)
    for (block.result in block.results) {
        rows <- block.result$rows
        deviance[rows] <- block.result$deviance
        df.residual[rows] <- block.result$df.residual
        df.total[rows] <- block.result$df.total
    }
    rownames(coefficients) <- rownames(y)
    rownames(unshrunk.coefficients) <- rownames(y)
    rownames(table) <- rownames(y)
    names(deviance) <- rownames(y)
    names(df.residual) <- rownames(y)
    names(df.total) <- rownames(y)

    list(
        table = table,
        coefficients = coefficients,
        unshrunk.coefficients = unshrunk.coefficients,
        deviance = deviance,
        df.residual = df.residual,
        df.total = df.total,
        comparison = treat.info$comparison
    )
}

.stream_treat_modern_test <- function(y,
                                      design,
                                      dispersion,
                                      offset,
                                      ave.log.cpm,
                                      prior.count,
                                      ql.stats,
                                      s2.fit,
                                      coef,
                                      contrast,
                                      lfc,
                                      null,
                                      block.size,
                                      workers = 1,
                                      ...) {
    dims <- .check_count_matrix(y)
    ntags <- dims[1L]
    nlibs <- dims[2L]
    block.ids <- .row_block_indices(ntags, block.size)
    treat.info <- .treat_design_info(design, coef, contrast)
    design.treat <- treat.info$design
    design0 <- design.treat[, -treat.info$coef, drop = FALSE]
    offset.adj <- lfc * log(2) * design.treat[, treat.info$coef]
    df.total.cap <- ntags * (nlibs - ncol(design))
    average.ql.dispersion <- ql.stats$average.ql.dispersion
    zscoreT <- utils::getFromNamespace("zscoreT", "limma")
    integratepnorm <- .edgeR_internal(".integratepnorm")

    block.results <- .block_lapply(seq_along(block.ids), workers, function(b) {
        rows <- block.ids[[b]]
        block <- .dense_slice(y, rows, TRUE)
        storage.mode(block) <- "numeric"
        rownames(block) <- rownames(y)[rows]
        colnames(block) <- colnames(y)

        block.offset <- .block_offset(offset, rows, ntags, nlibs)
        block.dispersion <- .block_dispersion(dispersion, rows, ntags)
        fit.dispersion <- block.dispersion / average.ql.dispersion
        # The modern QL statistics pass already fitted this alternative model
        # at `fit.dispersion`; reuse its coefficients instead of fitting it a
        # second time before the four threshold-constrained Treat fits.
        fit <- list(
            coefficients = ql.stats$coefficients[rows, , drop = FALSE],
            unshrunk.coefficients = ql.stats$unshrunk.coefficients[
                rows,
                ,
                drop = FALSE
            ],
            prior.count = prior.count
        )
        fit$dispersion <- if (length(block.dispersion) == 1L) {
            rep(block.dispersion, length(rows))
        } else {
            block.dispersion
        }
        fit$AveLogCPM <- ave.log.cpm[rows]
        fit$df.residual.adj <- ql.stats$df.residual.adj[rows]
        fit$deviance.adj <- ql.stats$deviance.adj[rows]
        fit$average.ql.dispersion <- average.ql.dispersion
        fit$df.prior <- .subset_prior_field(s2.fit$df.prior, rows, ntags)
        fit$s2.post <- s2.fit$var.post[rows]
        fit$s2.prior <- .subset_prior_field(s2.fit$var.prior, rows, ntags)
        fit$top.proportion <- NULL

        logfc <- .treat_logfc(fit, treat.info)
        fit0 <- edgeR::glmFit(
            block,
            design = design0,
            offset = block.offset + offset.adj,
            dispersion = fit.dispersion,
            prior.count = 0,
            ...
        )
        fit1 <- edgeR::glmFit(
            block,
            design = design.treat,
            offset = block.offset + offset.adj,
            dispersion = fit.dispersion,
            prior.count = 0,
            ...
        )
        z.left <- sqrt(pmax(0, fit0$deviance - fit1$deviance))

        fit0 <- edgeR::glmFit(
            block,
            design = design0,
            offset = block.offset - offset.adj,
            dispersion = fit.dispersion,
            prior.count = 0,
            ...
        )
        fit1 <- edgeR::glmFit(
            block,
            design = design.treat,
            offset = block.offset - offset.adj,
            dispersion = fit.dispersion,
            prior.count = 0,
            ...
        )
        z.right <- sqrt(pmax(0, fit0$deviance - fit1$deviance))

        swap <- z.left > z.right
        if (any(swap)) {
            tmp <- z.left[swap]
            z.left[swap] <- z.right[swap]
            z.right[swap] <- tmp
        }

        block.df.total <- fit$df.prior + fit$df.residual.adj
        block.df.total <- pmin(block.df.total, df.total.cap)
        z.left <- zscoreT(z.left / sqrt(fit$s2.post), df = block.df.total)
        z.right <- zscoreT(z.right / sqrt(fit$s2.post), df = block.df.total)

        within <- abs(logfc$unshrunk.logFC) <= lfc
        sgn <- 2 * within - 1
        z.left <- z.left * sgn

        p.value <- .treat_p_value(z.left, z.right, null, integratepnorm)

        tab <- data.frame(
            logCPM = fit$AveLogCPM,
            PValue = p.value,
            row.names = rownames(block)
        )
        if (!is.null(logfc$logFC)) {
            tab <- cbind(
                logFC = logfc$logFC,
                unshrunk.logFC = logfc$unshrunk.logFC,
                tab
            )
        } else {
            tab <- cbind(logFC = logfc$unshrunk.logFC, tab)
        }

        list(
            rows = rows,
            coefficients = fit$coefficients,
            unshrunk.coefficients = fit$unshrunk.coefficients,
            deviance = fit$deviance,
            df.residual = fit$df.residual,
            df.total = block.df.total,
            table = tab
        )
    })

    coefficients <- do.call(rbind, lapply(block.results, `[[`, "coefficients"))
    unshrunk.coefficients <- do.call(rbind, lapply(block.results, `[[`, "unshrunk.coefficients"))
    table <- do.call(rbind, lapply(block.results, `[[`, "table"))
    deviance <- rep_len(NA_real_, ntags)
    df.residual <- rep_len(NA_real_, ntags)
    df.total <- rep_len(NA_real_, ntags)
    for (block.result in block.results) {
        rows <- block.result$rows
        deviance[rows] <- block.result$deviance
        df.residual[rows] <- block.result$df.residual
        df.total[rows] <- block.result$df.total
    }
    rownames(coefficients) <- rownames(y)
    rownames(unshrunk.coefficients) <- rownames(y)
    rownames(table) <- rownames(y)
    names(deviance) <- rownames(y)
    names(df.residual) <- rownames(y)
    names(df.total) <- rownames(y)

    list(
        table = table,
        coefficients = coefficients,
        unshrunk.coefficients = unshrunk.coefficients,
        deviance = deviance,
        df.residual = df.residual,
        df.total = df.total,
        comparison = treat.info$comparison
    )
}

.treat_design_info <- function(design, coef, contrast) {
    coef.names <- colnames(design)
    if (is.null(coef.names)) {
        coef.names <- paste0("coef", seq_len(ncol(design)))
    }

    if (is.null(contrast)) {
        if (length(coef) > 1L) {
            coef <- unique(coef)
        }
        if (is.character(coef)) {
            check.coef <- coef %in% colnames(design)
            if (any(!check.coef)) {
                stop("One or more named coef arguments do not match a column of the design matrix.", call. = FALSE)
            }
            coef.name <- coef
            coef <- match(coef, colnames(design))
        } else {
            coef.name <- coef.names[coef]
        }
        if (length(coef) != 1L || is.na(coef) || coef < 1L || coef > ncol(design)) {
            stop("`coef` must identify one valid design column for Treat tests.", call. = FALSE)
        }
        return(list(
            coef = as.integer(coef),
            contrast = NULL,
            comparison = coef.name,
            design = design
        ))
    }

    contrast <- as.matrix(contrast)
    if (ncol(contrast) > 1L) {
        warning("Found more than one contrast vector in `contrast`. Using first column only.", call. = FALSE)
        contrast <- contrast[, 1L, drop = FALSE]
    }
    reform <- utils::getFromNamespace("contrastAsCoef", "limma")(
        design,
        contrast = contrast,
        first = TRUE
    )
    contrast <- drop(contrast)
    selected <- contrast != 0
    coef.name <- paste(paste(contrast[selected], coef.names[selected], sep = "*"), collapse = " ")
    list(
        coef = 1L,
        contrast = contrast,
        comparison = coef.name,
        design = reform$design
    )
}

.treat_logfc <- function(fit, treat.info) {
    shrunk <- fit$prior.count != 0
    if (is.null(treat.info$contrast)) {
        unshrunk.logFC <- fit$coefficients[, treat.info$coef, drop = FALSE] / log(2)
        if (shrunk) {
            logFC <- as.vector(unshrunk.logFC)
            unshrunk.logFC <- fit$unshrunk.coefficients[, treat.info$coef, drop = FALSE] / log(2)
        } else {
            logFC <- NULL
        }
    } else {
        unshrunk.logFC <- drop((fit$coefficients %*% treat.info$contrast) / log(2))
        if (shrunk) {
            logFC <- as.vector(unshrunk.logFC)
            unshrunk.logFC <- drop((fit$unshrunk.coefficients %*% treat.info$contrast) / log(2))
        } else {
            logFC <- NULL
        }
    }
    list(
        logFC = logFC,
        unshrunk.logFC = as.vector(unshrunk.logFC)
    )
}

.treat_p_value <- function(z.left, z.right, null, integratepnorm) {
    if (identical(null, "interval")) {
        cutoff <- 1.470402
        use.tails <- z.right + z.left > cutoff
        p.value <- rep_len(1, length(z.left))
        finite <- !is.na(use.tails)
        if (any(finite & use.tails)) {
            idx <- finite & use.tails
            p.value[idx] <- integratepnorm(-z.right[idx], -z.right[idx] + cutoff) +
                integratepnorm(z.left[idx] - cutoff, z.left[idx])
        }
        if (any(finite & !use.tails)) {
            idx <- finite & !use.tails
            p.value[idx] <- 2 * integratepnorm(-z.right[idx], z.left[idx])
        }
        return(pmin(pmax(p.value, 0), 1))
    }
    pmin(pmax(stats::pnorm(-z.right) + stats::pnorm(z.left), 0), 1)
}
