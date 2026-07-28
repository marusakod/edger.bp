test_that("bp_check_versions reports the active audited stack", {
    observed <- bp_check_versions()

    expect_s3_class(observed, "data.frame")
    expect_identical(observed$package, c("R", "edgeR", "limma", "BPCells"))
    expect_true(all(observed$version_supported))
    expect_true(all(observed$ok))
    expect_length(unique(observed$stack), 1L)
    expect_false(is.na(unique(observed$stack)))
})

test_that("all audited stacks are accepted as complete tuples", {
    stacks <- edger.bp:::.audited_stacks(include.candidate = FALSE)

    for (i in seq_len(nrow(stacks))) {
        installed <- unlist(stacks[i, -1L], use.names = TRUE)
        observed <- edger.bp:::.audit_installed_versions(installed, stacks)

        expect_true(all(observed$version_supported))
        expect_true(all(observed$ok))
        expect_identical(unique(observed$stack), stacks$stack[i])
    }
})

test_that("an unaudited mixture does not pass stack validation", {
    installed <- c(
        R = "4.4.3",
        edgeR = "4.10.1",
        limma = "3.62.0",
        BPCells = "0.3.1"
    )
    observed <- edger.bp:::.audit_installed_versions(
        installed,
        edger.bp:::.audited_stacks(include.candidate = FALSE)
    )

    expect_true(all(observed$version_supported))
    expect_false(any(observed$ok))
    expect_true(all(is.na(observed$stack)))
})

test_that("bp_check_versions validates error argument", {
    expect_error(bp_check_versions(error = NA), "`error` must be a logical scalar")
})
