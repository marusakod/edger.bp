.check_count_matrix <- function(y) {
    dims <- dim(y)
    if (length(dims) != 2L || anyNA(dims)) {
        stop("`y` must be a two-dimensional count matrix.", call. = FALSE)
    }
    invisible(dims)
}

.row_block_indices <- function(n, block.size) {
    block.size <- as.integer(block.size)
    if (length(block.size) != 1L || is.na(block.size) || block.size < 1L) {
        stop("`block.size` must be a positive integer.", call. = FALSE)
    }
    starts <- seq.int(1L, n, by = block.size)
    Map(seq.int, starts, pmin(starts + block.size - 1L, n))
}

.check_workers <- function(workers) {
    workers <- as.integer(workers)
    if (length(workers) != 1L || is.na(workers) || workers < 1L) {
        stop("`workers` must be a positive integer.", call. = FALSE)
    }
    workers
}

.check_verbose <- function(verbose) {
    .check_logical_scalar(verbose, "verbose")
}

.verbose_message <- function(verbose, ...) {
    if (isTRUE(verbose)) {
        message(format(Sys.time(), "%Y-%m-%d %H:%M:%S"), " | ", ...)
    }
    invisible(NULL)
}

.block_lapply <- function(x,
                          workers,
                          FUN,
                          ...,
                          .context = NULL,
                          .labels = NULL) {
    workers <- .check_workers(workers)
    n.tasks <- length(x)
    if (is.null(.labels)) {
        .labels <- sprintf("task %d of %d", seq_len(n.tasks), n.tasks)
    } else if (!is.character(.labels) || length(.labels) != n.tasks || anyNA(.labels)) {
        stop("`.labels` must be a character vector with one value per task.", call. = FALSE)
    }
    if (!is.null(.context) &&
        (!is.character(.context) || length(.context) != 1L || is.na(.context))) {
        stop("`.context` must be NULL or a character scalar.", call. = FALSE)
    }

    run.one <- function(i) {
        tryCatch(
            list(
                ok = TRUE,
                value = FUN(x[[i]], ...),
                task = i,
                pid = Sys.getpid()
            ),
            error = function(e) {
                condition.call <- conditionCall(e)
                list(
                    ok = FALSE,
                    value = NULL,
                    task = i,
                    pid = Sys.getpid(),
                    condition.class = paste(class(e), collapse = "/"),
                    condition.message = conditionMessage(e),
                    condition.call = if (is.null(condition.call)) {
                        NULL
                    } else {
                        paste(deparse(condition.call, width.cutoff = 500L), collapse = " ")
                    }
                )
            }
        )
    }

    task.ids <- seq_len(n.tasks)
    if (.bp_psock_state$worker && workers > 1L) {
        stop("Nested parallel workers are not supported inside a PSOCK task.")
    }
    if (!is.null(.bp_psock_state$pool)) {
        results <- .psock_block_lapply(x, workers, FUN, .context, list(...))
    } else if (workers == 1L || n.tasks <= 1L) {
        results <- lapply(task.ids, run.one)
    } else {
        results <- parallel::mclapply(
            task.ids,
            run.one,
            mc.cores = workers,
            mc.preschedule = TRUE
        )
    }

    ok <- vapply(
        results,
        function(result) is.list(result) && isTRUE(result$ok),
        logical(1)
    )
    if (!all(ok)) {
        failed <- which(!ok)
        details <- vapply(failed, function(i) {
            result <- results[[i]]
            if (!is.list(result) || is.null(result$ok)) {
                returned <- paste(as.character(result), collapse = " ")
                return(paste0(
                    "  ", .labels[[i]], "\n",
                    "    worker result: ", returned
                ))
            }

            lines <- c(
                paste0("  ", .labels[[i]]),
                paste0("    worker PID: ", result$pid),
                paste0("    condition: ", result$condition.class),
                paste0("    message: ", result$condition.message)
            )
            if (!is.null(result$condition.call)) {
                lines <- c(lines, paste0("    call: ", result$condition.call))
            }
            paste(lines, collapse = "\n")
        }, character(1))
        context <- if (is.null(.context) || !nzchar(.context)) {
            ""
        } else {
            paste0(" during ", .context)
        }
        stop(
            "Block processing failed", context, ":\n",
            paste(details, collapse = "\n"),
            call. = FALSE
        )
    }

    lapply(results, `[[`, "value")
}

.dense_slice <- function(y, i = TRUE, j = TRUE) {
    suppressWarnings(as.matrix(y[i, j, drop = FALSE]))
}

.check_lib_size <- function(lib.size, nlibs) {
    if (is.null(lib.size)) {
        return(NULL)
    }
    if (anyNA(lib.size)) {
        stop("NA library sizes are not permitted.", call. = FALSE)
    }
    if (length(lib.size) != nlibs) {
        if (length(lib.size) > 1L) {
            warning(
                "length(lib.size) does not match the number of samples; recycling",
                call. = FALSE
            )
        }
        lib.size <- rep_len(lib.size, nlibs)
    }
    as.numeric(lib.size)
}

.matrix_col_sums <- function(y) {
    BPCells::colSums(y)
}

.matrix_row_sums <- function(y) {
    BPCells::rowSums(y)
}

.block_offset <- function(offset, rows, nrows, nlibs) {
    if (is.null(offset)) {
        return(NULL)
    }
    if (is.vector(offset)) {
        if (length(offset) != nlibs) {
            stop("Vector `offset` must have length equal to the number of samples.", call. = FALSE)
        }
        return(as.numeric(offset))
    }
    offset <- as.matrix(offset)
    if (!identical(dim(offset), c(nrows, nlibs))) {
        stop("Matrix `offset` must have dimensions nrow(y) by ncol(y).", call. = FALSE)
    }
    offset[rows, , drop = FALSE]
}

.normalize_offset <- function(offset, nrows, nlibs) {
    if (is.null(offset) || is.vector(offset)) {
        return(offset)
    }
    offset <- as.matrix(offset)
    if (!identical(dim(offset), c(nrows, nlibs))) {
        stop("Matrix `offset` must have dimensions nrow(y) by ncol(y).", call. = FALSE)
    }
    first.row <- offset[1L, ]
    if (all(offset == matrix(first.row, nrow = nrows, ncol = nlibs, byrow = TRUE))) {
        return(as.numeric(first.row))
    }
    stop("Gene-specific offset matrices are not yet supported; supply sample-level offsets as a vector.", call. = FALSE)
}

.block_dispersion <- function(dispersion, rows, nrows) {
    if (is.null(dispersion)) {
        stop("`dispersion` must be supplied for blockwise GLM fitting.", call. = FALSE)
    }
    if (length(dispersion) == 1L) {
        return(as.numeric(dispersion))
    }
    if (length(dispersion) != nrows) {
        stop("`dispersion` must be scalar or have length nrow(y).", call. = FALSE)
    }
    as.numeric(dispersion[rows])
}

.check_logical_scalar <- function(x, name) {
    if (length(x) != 1L || is.na(x)) {
        stop(sprintf("`%s` must be a logical scalar.", name), call. = FALSE)
    }
    isTRUE(x)
}
