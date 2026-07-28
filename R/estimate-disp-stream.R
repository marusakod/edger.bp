#' Streamed Dispersion Estimation
#'
#' Estimates edgeR dispersions for BPCells-backed matrices without materializing
#' the full count matrix. The implementation mirrors the audited
#' `edgeR::estimateDisp.default()` GLM design path for edgeR 4.10.1 by building
#' the adjusted-profile likelihood grid in dense row blocks.
#'
#' @param y A genes-by-samples count matrix.
#' @param design Design matrix. The streamed implementation currently supports
#'   the GLM design path only.
#' @param group Optional group factor. Validated for consistency with edgeR but
#'   otherwise unused when `design` is supplied.
#' @param lib.size Optional library sizes. If `NULL`, column sums of `y` are
#'   used, matching `edgeR::estimateDisp.default()`.
#' @param offset Optional sample-level log-offset vector of length `ncol(y)`.
#'   Row-constant matrices with dimensions `nrow(y)` by `ncol(y)` are accepted
#'   and canonicalized to vectors. Gene-specific offset matrices are not yet
#'   supported.
#' @param prior.df,trend.method,tagwise,span,legacy.span,min.row.sum,grid.length,grid.range,robust,winsor.tail.p,tol
#'   Arguments corresponding to `edgeR::estimateDisp()`. `robust = TRUE` and
#'   observation weights are not yet supported by the streamed path.
#' @param weights Observation weights. Currently unsupported.
#' @param block.size Number of rows to process per block.
#' @param workers Number of row blocks to process in parallel.
#' @param verbose Logical scalar. If `TRUE`, emit low-noise progress messages.
#' @param ... Reserved for API compatibility with `edgeR::estimateDisp()`.
#'
#' @return A list matching the core dispersion fields returned by
#'   `edgeR::estimateDisp.default()`.
#' @export
bp_estimate_disp_stream <- function(y,
                                    design,
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
                                    weights = NULL,
                                    block.size = 1000,
                                    workers = 1,
                                    verbose = FALSE,
                                    ...) {
    .check_stream_disp_versions()
    dims <- .check_count_matrix(y)
    ntags <- dims[1L]
    nlibs <- dims[2L]
    workers <- .check_workers(workers)
    verbose <- .check_verbose(verbose)
    trend.method <- match.arg(
        trend.method,
        c("none", "loess", "locfit", "movingave", "locfit.mixed")
    )

    design <- as.matrix(design)
    if (nrow(design) != nlibs) {
        stop("`design` must have one row per sample.", call. = FALSE)
    }
    if (ncol(design) >= nlibs) {
        warning("No residual df: setting dispersion to NA", call. = FALSE)
        return(list(
            common.dispersion = NA_real_,
            trended.dispersion = NA_real_,
            tagwise.dispersion = NA_real_
        ))
    }
    if (!is.null(group) && length(group) != nlibs) {
        stop("Incorrect length of group.", call. = FALSE)
    }
    if (isTRUE(robust)) {
        stop("`robust = TRUE` is not yet supported by streamed dispersion estimation.", call. = FALSE)
    }
    if (!is.null(weights)) {
        stop("Observation `weights` are not yet supported by streamed dispersion estimation.", call. = FALSE)
    }
    if (length(grid.range) != 2L || anyNA(grid.range)) {
        stop("`grid.range` must be a numeric vector of length 2.", call. = FALSE)
    }
    grid.length <- as.integer(grid.length)
    if (length(grid.length) != 1L || is.na(grid.length) || grid.length < 2L) {
        stop("`grid.length` must be an integer of at least 2.", call. = FALSE)
    }

    if (ntags == 0L) {
        .verbose_message(verbose, "streamed dispersion: no genes supplied")
        return(list(span = span, prior.df = prior.df, prior.n = NULL))
    }

    if (is.null(lib.size)) {
        .verbose_message(verbose, "streamed dispersion: computing library sizes")
        lib.size <- bp_library_sizes(y)
    }
    lib.size <- .check_lib_size(lib.size, nlibs)
    offset <- .normalize_offset(offset, ntags, nlibs)
    if (is.null(offset)) {
        offset <- log(lib.size)
    }
    offset <- .check_sample_offset(offset, nlibs)

    sel <- as.vector(.matrix_row_sums(y) >= min.row.sum)
    spline.pts <- seq(from = grid.range[1L], to = grid.range[2L], length.out = grid.length)
    spline.disp <- 0.1 * 2^spline.pts
    .verbose_message(
        verbose,
        "streamed dispersion: building APL grid for ",
        sum(sel),
        " selected genes, ",
        grid.length,
        " grid points, ",
        workers,
        " worker(s)"
    )
    l0 <- .build_apl_grid_stream(
        y = y,
        design = design,
        offset = offset,
        sel = sel,
        spline.disp = spline.disp,
        block.size = block.size,
        workers = workers
    )
    .verbose_message(verbose, "streamed dispersion: estimating common dispersion")

    maximizeInterpolant <- .edgeR_internal("maximizeInterpolant")
    overall <- maximizeInterpolant(spline.pts, matrix(colSums(l0), nrow = 1L))
    common.dispersion <- 0.1 * 2^overall
    .verbose_message(
        verbose,
        "streamed dispersion: common dispersion = ",
        signif(common.dispersion, 4)
    )

    if (trend.method != "none") {
        .verbose_message(verbose, "streamed dispersion: computing AveLogCPM")
        AveLogCPM <- .stream_ave_log_cpm(
            y = y,
            lib.size = lib.size,
            dispersion = common.dispersion,
            block.size = block.size,
            workers = workers
        )
        .verbose_message(verbose, "streamed dispersion: fitting ", trend.method, " trend")
        WLEB <- .edgeR_internal("WLEB")
        wleb.args <- list(
            theta = spline.pts,
            loglik = l0,
            covariate = AveLogCPM[sel],
            trend.method = trend.method,
            span = span,
            overall = FALSE,
            individual = FALSE,
            m0.out = TRUE
        )
        if ("legacy.span" %in% names(formals(WLEB))) {
            wleb.args$legacy.span <- legacy.span
        }
        out.1 <- do.call(WLEB, wleb.args)
        span <- out.1$span
        m0 <- out.1$shared.loglik
        disp.trend <- 0.1 * 2^out.1$trend
        trended.dispersion <- rep(disp.trend[which.min(AveLogCPM[sel])], ntags)
        trended.dispersion[sel] <- disp.trend
    } else {
        .verbose_message(verbose, "streamed dispersion: skipping abundance trend")
        AveLogCPM <- NULL
        m0 <- matrix(colMeans(l0), ntags, length(spline.pts), byrow = TRUE)
        disp.trend <- common.dispersion
        trended.dispersion <- NULL
    }

    if (!tagwise) {
        .verbose_message(verbose, "streamed dispersion: finished without tagwise dispersions")
        out <- list(
            common.dispersion = common.dispersion,
            trended.dispersion = trended.dispersion
        )
        class(out) <- c("BPCellsStreamDisp", "list")
        return(out)
    }

    if (is.null(prior.df)) {
        .verbose_message(verbose, "streamed dispersion: estimating prior degrees of freedom")
        s2.info <- .stream_dispersion_s2(
            y = y,
            design = design,
            offset = offset,
            sel = sel,
            dispersion = disp.trend,
            block.size = block.size,
            workers = workers
        )
        s2.fit <- limma::squeezeVar(
            s2.info$s2,
            df = s2.info$df.residual,
            covariate = AveLogCPM[sel],
            robust = robust,
            winsor.tail.p = winsor.tail.p
        )
        prior.df <- s2.fit$df.prior
        .verbose_message(
            verbose,
            "streamed dispersion: prior.df = ",
            paste(signif(prior.df, 4), collapse = ",")
        )
    } else {
        .verbose_message(
            verbose,
            "streamed dispersion: using supplied prior.df = ",
            paste(signif(prior.df, 4), collapse = ",")
        )
    }

    ncoefs <- ncol(design)
    prior.n <- prior.df / (nlibs - ncoefs)
    if (trend.method != "none") {
        tagwise.dispersion <- trended.dispersion
    } else {
        tagwise.dispersion <- rep(common.dispersion, ntags)
    }

    too.large <- prior.n > 1e6
    if (!all(too.large)) {
        .verbose_message(verbose, "streamed dispersion: estimating tagwise dispersions")
        temp.n <- prior.n
        if (any(too.large)) {
            temp.n[too.large] <- 1e6
        }
        WLEB <- .edgeR_internal("WLEB")
        wleb.args <- list(
            theta = spline.pts,
            loglik = l0,
            prior.n = temp.n,
            covariate = AveLogCPM[sel],
            trend.method = trend.method,
            span = span,
            overall = FALSE,
            trend = FALSE,
            m0 = m0
        )
        if ("legacy.span" %in% names(formals(WLEB))) {
            wleb.args$legacy.span <- FALSE
        }
        out.2 <- do.call(WLEB, wleb.args)
        tagwise.dispersion[sel] <- 0.1 * 2^out.2$individual
    } else {
        .verbose_message(verbose, "streamed dispersion: skipping tagwise shrinkage because prior.n is very large")
    }

    out <- list(
        common.dispersion = common.dispersion,
        trended.dispersion = trended.dispersion,
        tagwise.dispersion = tagwise.dispersion,
        span = span,
        prior.df = prior.df,
        prior.n = prior.n
    )
    class(out) <- c("BPCellsStreamDisp", "list")
    .verbose_message(verbose, "streamed dispersion: finished")
    out
}

.build_apl_grid_stream <- function(y,
                                   design,
                                   offset,
                                   sel,
                                   spline.disp,
                                   block.size,
                                   workers = 1) {
    dims <- .check_count_matrix(y)
    ntags <- dims[1L]
    nlibs <- dims[2L]
    selected.rows <- which(sel)
    grid.length <- length(spline.disp)
    l0 <- matrix(0, length(selected.rows), grid.length)
    if (!length(selected.rows)) {
        return(l0)
    }

    comboGroups <- .edgeR_internal(".comboGroups")
    adjustedProfileLik <- .edgeR_internal("adjustedProfileLik")
    selected.blocks <- .row_block_indices(length(selected.rows), block.size)

    block.results <- .block_lapply(selected.blocks, workers, function(selected.pos) {
        rows <- selected.rows[selected.pos]
        block <- .dense_slice(y, rows, TRUE)
        storage.mode(block) <- "numeric"
        rownames(block) <- rownames(y)[rows]
        colnames(block) <- colnames(y)

        block.offset <- .block_offset(offset, rows, ntags, nlibs)
        glmfit <- edgeR::glmFit(
            block,
            design = design,
            offset = block.offset,
            dispersion = 0.05,
            prior.count = 0
        )
        zerofit <- (glmfit$counts < 1e-4 & glmfit$fitted.values < 1e-4)
        by.group <- comboGroups(zerofit)
        l0.block <- matrix(0, length(rows), grid.length)

        for (subg in by.group) {
            cur.nzero <- !zerofit[subg[1L], ]
            if (!any(cur.nzero)) {
                next
            }
            if (all(cur.nzero)) {
                redesign <- design
            } else {
                redesign <- design[cur.nzero, , drop = FALSE]
                QR <- qr(redesign)
                redesign <- redesign[, QR$pivot[seq_len(QR$rank)], drop = FALSE]
                if (nrow(redesign) == ncol(redesign)) {
                    next
                }
            }

            cury <- block[subg, cur.nzero, drop = FALSE]
            curo <- block.offset[cur.nzero]
            last.beta <- NULL
            for (i in seq_along(spline.disp)) {
                apl <- adjustedProfileLik(
                    spline.disp[i],
                    y = cury,
                    design = redesign,
                    offset = curo,
                    weights = NULL,
                    start = last.beta,
                    get.coef = TRUE
                )
                l0.block[subg, i] <- apl$apl
                last.beta <- apl$beta
            }
        }
        list(selected.pos = selected.pos, l0 = l0.block)
    })
    for (block.result in block.results) {
        l0[block.result$selected.pos, ] <- block.result$l0
    }

    l0
}

.stream_ave_log_cpm <- function(y,
                                lib.size,
                                dispersion,
                                block.size,
                                workers = 1) {
    dims <- .check_count_matrix(y)
    ntags <- dims[1L]
    nlibs <- dims[2L]
    lib.size <- .check_lib_size(lib.size, nlibs)
    block.ids <- .row_block_indices(ntags, block.size)
    ave <- rep_len(NA_real_, ntags)

    block.results <- .block_lapply(block.ids, workers, function(rows) {
        block <- .dense_slice(y, rows, TRUE)
        storage.mode(block) <- "numeric"
        block.dispersion <- .block_dispersion(dispersion, rows, ntags)
        list(
            rows = rows,
            ave = edgeR::aveLogCPM(
                block,
                lib.size = lib.size,
                dispersion = block.dispersion
            )
        )
    })
    for (block.result in block.results) {
        ave[block.result$rows] <- block.result$ave
    }

    names(ave) <- rownames(y)
    ave
}

.stream_dispersion_s2 <- function(y,
                                  design,
                                  offset,
                                  sel,
                                  dispersion,
                                  block.size,
                                  workers = 1) {
    dims <- .check_count_matrix(y)
    ntags <- dims[1L]
    nlibs <- dims[2L]
    selected.rows <- which(sel)
    s2 <- rep_len(NA_real_, length(selected.rows))
    df.residual <- rep_len(NA_real_, length(selected.rows))
    if (!length(selected.rows)) {
        return(list(s2 = s2, df.residual = df.residual))
    }

    residDF <- .edgeR_internal(".residDF")
    selected.blocks <- .row_block_indices(length(selected.rows), block.size)

    block.results <- .block_lapply(selected.blocks, workers, function(selected.pos) {
        rows <- selected.rows[selected.pos]
        block <- .dense_slice(y, rows, TRUE)
        storage.mode(block) <- "numeric"
        rownames(block) <- rownames(y)[rows]
        colnames(block) <- colnames(y)

        block.offset <- .block_offset(offset, rows, ntags, nlibs)
        block.dispersion <- if (length(dispersion) == 1L) {
            as.numeric(dispersion)
        } else {
            as.numeric(dispersion[selected.pos])
        }
        glmfit <- edgeR::glmFit(
            block,
            offset = block.offset,
            design = design,
            dispersion = block.dispersion,
            prior.count = 0
        )
        zerofit <- (glmfit$counts < 1e-4 & glmfit$fitted.values < 1e-4)
        block.df <- residDF(zerofit, design)
        block.s2 <- glmfit$deviance / block.df
        block.s2[block.df == 0] <- 0
        block.s2 <- pmax(block.s2, 0)

        list(
            selected.pos = selected.pos,
            s2 = block.s2,
            df.residual = block.df
        )
    })
    for (block.result in block.results) {
        s2[block.result$selected.pos] <- block.result$s2
        df.residual[block.result$selected.pos] <- block.result$df.residual
    }

    list(s2 = s2, df.residual = df.residual)
}

.check_sample_offset <- function(offset, nlibs) {
    if (!is.numeric(offset) || is.null(dim(offset))) {
        if (length(offset) != nlibs) {
            stop("Vector `offset` must have length equal to the number of samples.", call. = FALSE)
        }
        offset <- as.numeric(offset)
        if (anyNA(offset) || any(!is.finite(offset))) {
            stop("`offset` must contain finite values.", call. = FALSE)
        }
        return(offset)
    }
    stop("Gene-specific offset matrices are not yet supported; supply sample-level offsets as a vector.", call. = FALSE)
}

.edgeR_internal <- function(name) {
    utils::getFromNamespace(name, "edgeR")
}

.check_stream_disp_versions <- function() {
    versions <- bp_check_versions(error = FALSE)
    if (!all(versions$ok)) {
        stop(
            paste0(
                "Streamed dispersion estimation requires a complete audited ",
                "dependency stack; installed versions are ",
                paste(versions$package, versions$installed, sep = " ", collapse = ", "),
                ". Supported stacks are: ",
                .format_audited_stacks(),
                "."
            ),
            call. = FALSE
        )
    }
    invisible(TRUE)
}
