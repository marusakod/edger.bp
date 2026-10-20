# Run from the package root after installing this checkout and its dependencies.
# Evaluate the vignette itself so the README table uses the same example.
stopifnot(all(edger.bp::bp_check_versions(error = TRUE)$ok))

vignette_env <- new.env(parent = globalenv())
rendered_md <- tempfile(fileext = ".md")
knitr::opts_chunk$set(error = FALSE)
tryCatch(
    knitr::knit(
        "vignettes/edgeR-functionality.Rmd",
        output = rendered_md,
        envir = vignette_env,
        quiet = TRUE
    ),
    finally = unlink(rendered_md)
)
comparison <- vignette_env$comparison
stopifnot(
    is.data.frame(comparison),
    nrow(comparison) == 6L,
    all(is.finite(comparison$max_abs_diff)),
    all(comparison$max_abs_diff < 1e-6)
)

comparison$max_abs_diff <- formatC(
    comparison$max_abs_diff, format = "e", digits = 2
)
table_md <- knitr::kable(
    comparison,
    format = "pipe",
    col.names = c("Quantity", "Maximum absolute difference")
)
versions <- sprintf(
    "Generated with edger.bp %s, R %s, edgeR %s, limma %s, and BPCells %s.",
    packageVersion("edger.bp"),
    getRversion(),
    packageVersion("edgeR"),
    packageVersion("limma"),
    packageVersion("BPCells")
)

readme <- readLines("README.md")
start <- which(readme == "<!-- edger-comparison:start -->")
end <- which(readme == "<!-- edger-comparison:end -->")
stopifnot(length(start) == 1L, length(end) == 1L, start < end)
writeLines(
    c(readme[seq_len(start)], "", table_md, "", versions, "",
      readme[seq.int(end, length(readme))]),
    "README.md"
)
