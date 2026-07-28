#' Streaming TMM Normalization Factors
#'
#' Computes TMM normalization factors for a matrix-like object while avoiding a
#' full dense realization of the genes-by-samples count matrix.
#'
#' @param y A genes-by-samples count matrix.
#' @param lib.size Optional library sizes. If `NULL`, column sums of `y` are
#'   used.
#' @param method Normalization method. Only `"TMM"` and `"none"` are currently
#'   implemented.
#' @param refColumn Optional 1-based reference column for TMM.
#' @param logratioTrim,sumTrim,doWeighting,Acutoff TMM parameters with the same
#'   meaning as in [edgeR::normLibSizes()].
#' @param workers Number of columns to process in parallel for TMM.
#' @param verbose Logical scalar. If `TRUE`, emit low-noise progress messages.
#'
#' @return A numeric vector of normalization factors, one per sample.
#' @export
bp_calc_norm_factors <- function(y,
                                 lib.size = NULL,
                                 method = c("TMM", "none"),
                                 refColumn = NULL,
                                 logratioTrim = 0.3,
                                 sumTrim = 0.05,
                                 doWeighting = TRUE,
                                 Acutoff = -1e10,
                                 workers = 1,
                                 verbose = FALSE) {
    dims <- .check_count_matrix(y)
    ntags <- dims[1L]
    nlibs <- dims[2L]
    workers <- .check_workers(workers)
    verbose <- .check_verbose(verbose)

    method <- match.arg(method)
    lib.size <- .check_lib_size(lib.size, nlibs)
    if (is.null(lib.size)) {
        .verbose_message(verbose, "TMM: computing library sizes for ", nlibs, " samples")
        lib.size <- bp_library_sizes(y)
    }

    if (method == "none" || ntags == 0L || nlibs == 1L) {
        .verbose_message(verbose, method, ": returning unit normalization factors")
        factors <- rep_len(1, nlibs)
        names(factors) <- colnames(y)
        return(factors)
    }

    .verbose_message(verbose, "TMM: finding genes with nonzero counts")
    positive.counts <- .matrix_row_sums(y > 0)
    norm.rows <- which(positive.counts != 0L)

    if (length(norm.rows) == 0L) {
        .verbose_message(verbose, "TMM: no nonzero genes found; returning unit normalization factors")
        factors <- rep_len(1, nlibs)
        names(factors) <- colnames(y)
        return(factors)
    }

    if (is.null(refColumn)) {
        .verbose_message(
            verbose,
            "TMM: selecting reference column across ",
            nlibs,
            " samples and ",
            length(norm.rows),
            " nonzero genes"
        )
        f75 <- .tmm_upper_quartiles(y, norm.rows, lib.size, workers = workers)
        if (stats::median(f75) < 1e-20) {
            refColumn <- .tmm_sqrt_reference(y, norm.rows, workers = workers)
        } else {
            refColumn <- which.min(abs(f75 - mean(f75)))
        }
    }

    if (length(refColumn) != 1L || is.na(refColumn) ||
        refColumn < 1L || refColumn > nlibs) {
        stop("`refColumn` must be a valid 1-based column index.", call. = FALSE)
    }

    ref.name <- colnames(y)
    ref.label <- if (!is.null(ref.name)) {
        paste0(" (", ref.name[[refColumn]], ")")
    } else {
        ""
    }
    .verbose_message(verbose, "TMM: using reference column ", refColumn, ref.label)
    .verbose_message(
        verbose,
        "TMM: computing normalization factors for ",
        nlibs,
        " samples with ",
        workers,
        " worker(s)"
    )
    ref <- as.numeric(.dense_slice(y, norm.rows, refColumn)[, 1L])
    factors <- unlist(.block_lapply(seq_len(nlibs), workers, function(j) {
        obs <- as.numeric(.dense_slice(y, norm.rows, j)[, 1L])
        .calc_tmm_factor(
            obs = obs,
            ref = ref,
            libsize.obs = lib.size[j],
            libsize.ref = lib.size[refColumn],
            logratioTrim = logratioTrim,
            sumTrim = sumTrim,
            doWeighting = doWeighting,
            Acutoff = Acutoff
        )
    }), use.names = FALSE)

    .verbose_message(verbose, "TMM: rescaling normalization factors")
    factors <- factors / exp(mean(log(factors)))
    names(factors) <- colnames(y)
    .verbose_message(verbose, "TMM: finished")
    factors
}

.tmm_upper_quartiles <- function(y, rows, lib.size, workers = 1) {
    nlibs <- ncol(y)
    out <- unlist(.block_lapply(seq_len(nlibs), workers, function(j) {
        column <- as.numeric(.dense_slice(y, rows, j)[, 1L])
        stats::quantile(column, probs = 0.75) / lib.size[j]
    }), use.names = FALSE)
    out
}

.tmm_sqrt_reference <- function(y, rows, workers = 1) {
    nlibs <- ncol(y)
    sums <- unlist(.block_lapply(seq_len(nlibs), workers, function(j) {
        column <- as.numeric(.dense_slice(y, rows, j)[, 1L])
        sum(sqrt(column))
    }), use.names = FALSE)
    which.max(sums)
}

.calc_tmm_factor <- function(obs,
                             ref,
                             libsize.obs,
                             libsize.ref,
                             logratioTrim,
                             sumTrim,
                             doWeighting,
                             Acutoff) {
    obs <- as.numeric(obs)
    ref <- as.numeric(ref)
    nO <- libsize.obs
    nR <- libsize.ref

    logR <- log2((obs / nO) / (ref / nR))
    absE <- (log2(obs / nO) + log2(ref / nR)) / 2
    v <- (nO - obs) / nO / obs + (nR - ref) / nR / ref

    fin <- is.finite(logR) & is.finite(absE) & (absE > Acutoff)
    logR <- logR[fin]
    absE <- absE[fin]
    v <- v[fin]

    if (max(abs(logR)) < 1e-6) {
        return(1)
    }

    n <- length(logR)
    loL <- floor(n * logratioTrim) + 1
    hiL <- n + 1 - loL
    loS <- floor(n * sumTrim) + 1
    hiS <- n + 1 - loS
    keep <- (rank(logR) >= loL & rank(logR) <= hiL) &
        (rank(absE) >= loS & rank(absE) <= hiS)

    if (doWeighting) {
        f <- sum(logR[keep] / v[keep], na.rm = TRUE) /
            sum(1 / v[keep], na.rm = TRUE)
    } else {
        f <- mean(logR[keep], na.rm = TRUE)
    }

    if (is.na(f)) {
        f <- 0
    }
    2^f
}
