test_that("bp_library_sizes matches dense column sums", {
    fixture <- make_counts()
    bp <- make_bp_matrix(fixture$counts)

    expect_equal(
        bp_library_sizes(bp),
        as.numeric(colSums(fixture$counts)),
        tolerance = 0,
        ignore_attr = TRUE
    )
})
