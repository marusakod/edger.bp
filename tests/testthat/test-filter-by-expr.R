test_that("bp_filter_by_expr matches edgeR with a design", {
    fixture <- make_counts(seed = 2)
    bp <- make_bp_matrix(fixture$counts)

    observed <- bp_filter_by_expr(
        bp,
        design = fixture$design,
        block.size = 13
    )
    expected <- edgeR::filterByExpr(
        fixture$counts,
        design = fixture$design
    )

    expect_identical(observed, expected)
})

test_that("bp_filter_by_expr matches edgeR with groups and supplied library sizes", {
    fixture <- make_counts(seed = 3)
    bp <- make_bp_matrix(fixture$counts)
    lib.size <- colSums(fixture$counts) * seq(0.8, 1.2, length.out = ncol(fixture$counts))

    observed <- bp_filter_by_expr(
        bp,
        group = fixture$group,
        lib.size = lib.size,
        block.size = 11
    )
    expected <- edgeR::filterByExpr(
        fixture$counts,
        group = fixture$group,
        lib.size = lib.size
    )

    expect_identical(observed, expected)
})
