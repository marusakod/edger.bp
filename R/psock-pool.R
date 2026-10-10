.bp_psock_state <- new.env(parent = emptyenv())
.bp_psock_state$pool <- NULL
.bp_psock_state$worker <- FALSE

#' Create a reusable disk-backed PSOCK worker pool
#'
#' Each worker opens the same read-only BPCells directory and row/column subset.
#' Use [bp_with_psock_pool()] around successive edgeR fits to reuse the workers.
#' Register `on.exit(bp_stop_psock_pool(pool), add = TRUE)` immediately after
#' creation. The pool is owned by the creating R process; do not export it.
#' @param matrix_dir Native BPCells matrix directory.
#' @param rows,cols Indices selecting the analysis matrix, in analysis order.
#' @param workers Number of independent R worker sessions.
#' @param verbose Print pool identities and per-pass task/thread diagnostics.
#' @param timeout Socket timeout in seconds (also bounds an unresponsive call).
#' @return An environment representing the pool; `pool$counts` is the matching
#'   matrix in the parent session.
#' @export
bp_start_psock_pool <- function(matrix_dir, rows = TRUE, cols = TRUE,
                                workers = 1L, verbose = TRUE, timeout = 3600) {
    workers <- .check_workers(workers)
    if (.bp_psock_state$worker || !is.null(.bp_psock_state$pool))
        stop("Nested PSOCK pools are not supported.")
    if (length(timeout) != 1L || !is.finite(timeout) || timeout <= 0)
        stop("timeout must be a positive finite number of seconds")
    matrix_dir <- normalizePath(matrix_dir, mustWork = TRUE)
    pool <- new.env(parent = emptyenv())
    pool$owner <- Sys.getpid()
    pool$closed <- FALSE
    pool$workers <- workers
    pool$verbose <- isTRUE(verbose)
    pool$pass <- 0L
    pool$counts <- BPCells::open_matrix_dir(matrix_dir)[rows, cols, drop = FALSE]
    pool$cluster <- list()
    success <- FALSE
    on.exit(if (!success) bp_stop_psock_pool(pool), add = TRUE)
    # Set these before Rscript starts, not after a BLAS library initializes.
    keys <- c("OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "OMP_NUM_THREADS",
              "VECLIB_MAXIMUM_THREADS", "NUMEXPR_NUM_THREADS")
    previous <- Sys.getenv(keys, unset = NA_character_)
    on.exit({
        Sys.unsetenv(keys[is.na(previous)])
        do.call(Sys.setenv, as.list(previous[!is.na(previous)]))
    }, add = TRUE)
    do.call(Sys.setenv, as.list(setNames(rep("1", length(keys)), keys)))
    # Sequential startup lets us close already-created workers on partial failure.
    for (i in seq_len(workers)) {
        node <- parallel::makePSOCKcluster(1L, rscript = file.path(R.home("bin"), "Rscript"),
            setup_strategy = "sequential", timeout = timeout)
        pool$cluster <- structure(c(unclass(pool$cluster), unclass(node)),
                                  class = c("SOCKcluster", "cluster"))
    }
    init <- function(libs, package_path, path, rows, cols) {
        .libPaths(libs)
        loadNamespace("edger.bp")
        if (normalizePath(find.package("edger.bp")) != package_path)
            stop("PSOCK worker loaded a different edger.bp installation")
        state <- get(".bp_psock_state", asNamespace("edger.bp"))
        state$worker <- TRUE
        state$counts <- BPCells::open_matrix_dir(path)[rows, cols, drop = FALSE]
        if (requireNamespace("data.table", quietly = TRUE)) data.table::setDTthreads(1L)
        Sys.getpid()
    }
    environment(init) <- baseenv()
    pool$pids <- unlist(parallel::clusterCall(pool$cluster, init, .libPaths(),
        normalizePath(find.package("edger.bp")), matrix_dir, rows, cols))
    success <- TRUE
    .verbose_message(pool$verbose, "PSOCK pool started: ", workers,
        " workers; PIDs ", paste(pool$pids, collapse = ","))
    .psock_parent_diagnostics(pool)
    pool
}

#' @rdname bp_start_psock_pool
#' @param pool Pool returned by `bp_start_psock_pool()`.
#' @export
bp_stop_psock_pool <- function(pool) {
    if (pool$owner != Sys.getpid()) stop("Only the pool owner may stop workers.")
    if (pool$closed) return(invisible(NULL))
    pool$closed <- TRUE
    for (node in pool$cluster) {
        tryCatch(parallel::stopCluster(structure(list(node), class = c("SOCKcluster", "cluster"))),
            error = function(e) warning("PSOCK shutdown: ", conditionMessage(e), call. = FALSE))
    }
    .verbose_message(pool$verbose, "PSOCK pool stopped after ", pool$pass, " passes")
    invisible(NULL)
}

#' @rdname bp_start_psock_pool
#' @param code Expression containing fits of `pool$counts`. Existing numeric
#'   `workers` arguments should equal the pool size. Nested pools are rejected.
#' @export
bp_with_psock_pool <- function(pool, code) {
    if (pool$closed || pool$owner != Sys.getpid()) stop("PSOCK pool is closed or has another owner.")
    if (!is.null(.bp_psock_state$pool)) stop("Nested PSOCK scopes are not supported.")
    .bp_psock_state$pool <- pool
    on.exit(.bp_psock_state$pool <- NULL, add = TRUE)
    force(code)
}

# Capture only free bindings, never the full calling frame (which may contain
# previous large fits). Namespace functions keep their namespace references.
# BPCells matrices and registered native symbols are rebuilt inside workers.
.psock_pack_function <- function(fun, pool) {
    # Source references retain parser environments and whole source files.
    # Workers need executable code, not those potentially large environments.
    fun <- utils::removeSource(fun)
    env <- environment(fun)
    # codetools does not list a lexically inherited `...` as a free variable.
    uses_dots <- "..." %in% all.names(body(fun))
    globals <- unique(unlist(suppressWarnings(codetools::findGlobals(fun, merge = FALSE))))
    bindings <- list()
    dots <- if (uses_dots && !"..." %in% names(formals(fun))) eval(quote(list(...)), env) else list()
    for (name in setdiff(globals, "...")) {
        if (!exists(name, env, inherits = TRUE)) next
        value <- get(name, env, inherits = TRUE)
        if (inherits(value, "IterableMatrix")) {
            if (!identical(value, pool$counts))
                stop("PSOCK callback matrix differs from pool$counts: ", name)
            value <- structure(list(kind = "counts"), class = "bp_psock_binding")
        } else if (inherits(value, "NativeSymbolInfo")) {
            package <- value$dll[["name"]]
            ns <- asNamespace(package)
            symbols <- ls(ns, all.names = TRUE)
            symbols <- symbols[vapply(symbols, function(n) identical(get(n, ns), value), logical(1))]
            if (!length(symbols)) stop("Cannot reconstruct native symbol in worker: ", value$name)
            value <- structure(list(kind = "native", binding = symbols[[1L]],
                package = package), class = "bp_psock_binding")
        } else if (is.function(value) && !is.primitive(value) &&
                   !isNamespace(environment(value)) && !identical(environment(value), baseenv())) {
            value <- structure(list(kind = "function", packet = .psock_pack_function(value, pool)),
                               class = "bp_psock_binding")
        }
        bindings[name] <- list(value)
    }
    if (uses_dots) {
        fml <- formals(fun)
        if (!"..." %in% names(fml)) formals(fun) <- c(fml, alist(... = ))
    }
    environment(fun) <- baseenv()
    list(fun = fun, bindings = bindings, dots = dots)
}

.psock_unpack_function <- function(packet) {
    env <- new.env(parent = baseenv())
    for (name in names(packet$bindings)) {
        value <- packet$bindings[[name]]
        if (inherits(value, "bp_psock_binding")) {
            value <- switch(value$kind,
                counts = .bp_psock_state$counts,
                native = get(value$binding, asNamespace(value$package), inherits = FALSE),
                "function" = .psock_unpack_function(value$packet)$fun)
        }
        assign(name, value, envir = env)
    }
    environment(packet$fun) <- env
    packet
}

.psock_task <- function(item) {
    tryCatch(list(ok = TRUE, value = do.call(.bp_psock_state$packet$fun,
        c(list(item$value), .bp_psock_state$packet$dots)), task = item$index, pid = Sys.getpid()),
        error = function(e) list(ok = FALSE, value = NULL, task = item$index,
            pid = Sys.getpid(), condition.class = paste(class(e), collapse = "/"),
            condition.message = conditionMessage(e),
            condition.call = paste(deparse(conditionCall(e)), collapse = " ")))
}

.psock_snapshot <- function() {
    status <- if (file.exists("/proc/self/status")) readLines("/proc/self/status") else character()
    list(pid = Sys.getpid(), threads = sub("^Threads:\\s*", "", grep("^Threads:", status, value = TRUE)),
         rss = grep("^VmRSS:", status, value = TRUE))
}

.psock_parent_diagnostics <- function(pool) {
    if (!pool$verbose) return(invisible(NULL))
    read <- function(path) if (file.exists(path)) {
        tryCatch(readLines(path, warn = FALSE), error = function(e) character())
    } else character()
    limit <- grep("^Max processes", read("/proc/self/limits"), value = TRUE)
    .verbose_message(TRUE, "PSOCK parent ", Sys.getpid(), "; ", paste(limit, collapse = "; "))
    cg <- sub("^0::", "", grep("^0::", read("/proc/self/cgroup"), value = TRUE))
    if (length(cg) == 1L) {
        path <- paste0("/sys/fs/cgroup", cg)
        while (startsWith(path, "/sys/fs/cgroup") && path != "/sys/fs") {
            current <- read(file.path(path, "pids.current"))
            if (length(current)) .verbose_message(TRUE, path, ": pids.current=", current,
                " pids.max=", read(file.path(path, "pids.max")),
                " events=", paste(read(file.path(path, "pids.events")), collapse = ";"))
            path <- dirname(path)
        }
    }
    invisible(NULL)
}

.psock_install <- function(packet) {
    .bp_psock_state$packet <- .psock_unpack_function(packet)
    NULL
}

.psock_apply_index <- function(i, task_values) {
    .psock_task(list(index = i, value = task_values[[i]]))
}

.psock_finish <- function() {
    .bp_psock_state$packet <- NULL
    .psock_snapshot()
}

.psock_block_lapply <- function(x, workers, FUN, context, extra_args = list()) {
    pool <- .bp_psock_state$pool
    if (pool$closed || workers != pool$workers) stop("PSOCK pool size/worker count mismatch or closed pool.")
    packet <- .psock_pack_function(FUN, pool)
    packet$dots <- c(packet$dots, extra_args)
    pool$pass <- pool$pass + 1L
    label <- paste0("PSOCK pass ", pool$pass, if (!is.null(context)) paste0(" (", context, ")"))
    .verbose_message(pool$verbose, label, ": ", length(x), " tasks")
    .psock_parent_diagnostics(pool)
    tryCatch({
        parallel::clusterCall(pool$cluster, .psock_install, packet)
        result <- parallel::parLapply(pool$cluster, seq_along(x), .psock_apply_index, task_values = x)
        snapshots <- parallel::clusterCall(pool$cluster, .psock_finish)
        if (!identical(unlist(lapply(snapshots, `[[`, "pid")), pool$pids))
            stop("PSOCK worker identities changed")
        .verbose_message(pool$verbose, label, " complete; worker threads: ",
            paste(vapply(snapshots, function(s) paste(s$threads, collapse = ""), ""), collapse = ","))
        result
    }, error = function(e) {
        .psock_parent_diagnostics(pool)
        bp_stop_psock_pool(pool)
        stop(label, " failed: ", conditionMessage(e), call. = FALSE)
    })
}

#' Apply a callback using the active pool (or the usual numeric-worker backend)
#' @param x Tasks.
#' @param workers Number of workers.
#' @param FUN Function of one task. In a pool scope, disk-backed matrix bindings
#'   must equal `pool$counts`; callbacks must not launch nested parallel work.
#' @param context Optional error/progress label.
#' @return List of results in task order.
#' @export
bp_block_lapply <- function(x, workers, FUN, context = NULL) {
    .block_lapply(x, workers, FUN, .context = context)
}
