#' Check Audited Package Versions
#'
#' Reports whether the installed R, edgeR, limma, and BPCells versions match
#' one complete stack audited for exact streamed edgeR compatibility.
#'
#' @param error Logical scalar. If `TRUE`, throw an error when the installed
#'   versions do not match a complete audited stack. If `FALSE`, return the
#'   check table without error.
#'
#' @return A data frame with package names, installed versions, expected
#'   versions, individual version support, complete-stack match status, and
#'   the matched stack name. A version can be individually supported while
#'   `ok` is `FALSE` when the installed versions form an unaudited mixture.
#' @export
bp_check_versions <- function(error = FALSE) {
    if (length(error) != 1L || is.na(error)) {
        stop("`error` must be a logical scalar.", call. = FALSE)
    }
    error <- isTRUE(error)

    out <- .audit_installed_versions(.installed_dependency_versions())

    if (error && !all(out$ok)) {
        stop(
            paste0(
                "Installed versions do not match any complete audited stack: ",
                paste(out$package, out$installed, sep = " ", collapse = ", "),
                ". Supported stacks are: ",
                .format_audited_stacks(),
                "."
            ),
            call. = FALSE
        )
    }

    out
}

.installed_dependency_versions <- function() {
    c(
        R = paste(R.version$major, R.version$minor, sep = "."),
        edgeR = as.character(utils::packageVersion("edgeR")),
        limma = as.character(utils::packageVersion("limma")),
        BPCells = as.character(utils::packageVersion("BPCells"))
    )
}

.audited_stacks <- function(include.candidate = FALSE) {
    stacks <- data.frame(
        stack = c("development", "conda-r45", "conda-scprocess"),
        R = c("4.6.1", "4.5.3", "4.4.3"),
        edgeR = c("4.10.1", "4.8.2", "4.4.0"),
        limma = c("3.68.4", "3.66.0", "3.62.0"),
        BPCells = c("0.3.1", "0.3.1", "0.3.1"),
        stringsAsFactors = FALSE,
        check.names = FALSE
    )
    if (isTRUE(include.candidate)) {
        candidate <- .installed_dependency_versions()
        packages <- setdiff(names(stacks), "stack")
        already.audited <- any(vapply(
            seq_len(nrow(stacks)),
            function(i) {
                identical(
                    unname(candidate[packages]),
                    unname(unlist(stacks[i, packages], use.names = FALSE))
                )
            },
            logical(1)
        ))
        if (!already.audited) {
            stacks <- rbind(
                stacks,
                data.frame(
                    stack = "candidate",
                    as.list(candidate),
                    stringsAsFactors = FALSE,
                    check.names = FALSE
                )
            )
        }
    }
    stacks
}

.audit_installed_versions <- function(
    installed,
    stacks = .audited_stacks(
        include.candidate = identical(
            Sys.getenv("EDGER_BP_COMPATIBILITY_TEST"),
            "candidate"
        )
    )
) {
    packages <- setdiff(names(stacks), "stack")
    if (
        is.null(names(installed)) ||
        !setequal(names(installed), packages) ||
        anyDuplicated(names(installed))
    ) {
        stop(
            "`installed` must contain one named version for each audited dependency.",
            call. = FALSE
        )
    }
    installed <- stats::setNames(as.character(installed[packages]), packages)

    matches <- vapply(
        seq_len(nrow(stacks)),
        function(i) {
            identical(
                unname(installed),
                unname(unlist(stacks[i, packages], use.names = FALSE))
            )
        },
        logical(1)
    )
    matched <- which(matches)
    matched.stack <- if (length(matched) == 1L) stacks$stack[matched] else NA_character_

    if (length(matched) == 1L) {
        expected <- stats::setNames(
            unlist(stacks[matched, packages], use.names = FALSE),
            packages
        )
    } else {
        expected <- vapply(
            packages,
            function(package) paste(unique(stacks[[package]]), collapse = " | "),
            character(1)
        )
    }
    individually.supported <- vapply(
        packages,
        function(package) installed[[package]] %in% stacks[[package]],
        logical(1)
    )

    data.frame(
        package = packages,
        installed = unname(installed),
        expected = unname(expected),
        version_supported = unname(individually.supported),
        ok = rep(length(matched) == 1L, length(packages)),
        stack = rep(matched.stack, length(packages)),
        stringsAsFactors = FALSE
    )
}

.format_audited_stacks <- function(stacks = .audited_stacks()) {
    packages <- setdiff(names(stacks), "stack")
    paste(
        vapply(
            seq_len(nrow(stacks)),
            function(i) {
                paste0(
                    stacks$stack[i],
                    " (",
                    paste(packages, unlist(stacks[i, packages]), sep = " ", collapse = ", "),
                    ")"
                )
            },
            character(1)
        ),
        collapse = "; "
    )
}
