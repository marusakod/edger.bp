test_that("bp_calc_norm_factors matches edgeR TMM", {
    fixture <- make_counts(seed = 4)
    bp <- make_bp_matrix(fixture$counts)

    observed <- bp_calc_norm_factors(bp, method = "TMM")
    expected <- edgeR::normLibSizes(fixture$counts, method = "TMM")

    expect_equal(observed, expected, tolerance = 1e-12)
})

test_that("bp_calc_norm_factors gives the same TMM factors in parallel", {
    fixture <- make_counts(seed = 44, nsamples = 8)
    bp <- make_bp_matrix(fixture$counts)

    serial <- bp_calc_norm_factors(bp, method = "TMM", workers = 1)
    parallel <- bp_calc_norm_factors(bp, method = "TMM", workers = 2)

    expect_equal(parallel, serial, tolerance = 1e-12)
})

test_that("bp_calc_norm_factors reports low-noise progress when verbose", {
    fixture <- make_counts(seed = 45, nsamples = 6)
    bp <- make_bp_matrix(fixture$counts)

    messages <- capture.output(
        observed <- bp_calc_norm_factors(bp, method = "TMM", verbose = TRUE),
        type = "message"
    )

    expect_length(observed, ncol(fixture$counts))
    expect_true(any(grepl("TMM: selecting reference column", messages, fixed = TRUE)))
    expect_true(any(grepl("TMM: computing normalization factors", messages, fixed = TRUE)))
    expect_true(any(grepl("TMM: finished", messages, fixed = TRUE)))
})

test_that("bp_calc_norm_factors matches edgeR TMM after filterByExpr workflow", {
    fixture <- make_counts(seed = 5)
    bp <- make_bp_matrix(fixture$counts)

    dge <- edgeR::DGEList(fixture$counts)
    keep <- edgeR::filterByExpr(dge, design = fixture$design)
    dge <- dge[keep, , keep.lib.sizes = FALSE]
    dge <- edgeR::normLibSizes(dge, method = "TMM")

    bp_keep <- bp_filter_by_expr(bp, design = fixture$design, block.size = 9)
    bp_filtered <- bp[bp_keep, , drop = FALSE]
    lib.size <- bp_library_sizes(bp_filtered)
    observed <- bp_calc_norm_factors(bp_filtered, lib.size = lib.size)

    expect_identical(bp_keep, keep)
    expect_equal(unname(lib.size), dge$samples$lib.size, tolerance = 0)
    expect_equal(unname(observed), dge$samples$norm.factors, tolerance = 1e-12)
})

test_that("bp_calc_norm_factors supports method none", {
    fixture <- make_counts(seed = 6)
    bp <- make_bp_matrix(fixture$counts)

    observed <- bp_calc_norm_factors(bp, method = "none")

    expect_equal(observed, rep(1, ncol(fixture$counts)), ignore_attr = TRUE)
})
