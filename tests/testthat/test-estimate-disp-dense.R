test_that("bp_estimate_disp_dense matches edgeR estimateDisp with prepared offsets", {
    fixture <- make_counts(ngenes = 120, nsamples = 6, seed = 14)
    bp <- make_bp_matrix(fixture$counts)
    prep <- bp_prepare_dge(bp, design = fixture$design, block.size = 17)
    dense <- fixture$counts[prep$keep, , drop = FALSE]

    observed <- bp_estimate_disp_dense(
        prep$counts,
        design = fixture$design,
        offset = prep$offset
    )
    expected <- edgeR::estimateDisp(
        dense,
        design = fixture$design,
        offset = prep$offset
    )

    expect_equal(observed$common.dispersion, expected$common.dispersion, tolerance = 1e-12)
    expect_equal(observed$trended.dispersion, expected$trended.dispersion, tolerance = 1e-12)
    expect_equal(observed$tagwise.dispersion, expected$tagwise.dispersion, tolerance = 1e-12)
    expect_equal(observed$span, expected$span, tolerance = 1e-12)
    expect_equal(observed$prior.df, expected$prior.df, tolerance = 1e-12)
    expect_equal(observed$prior.n, expected$prior.n, tolerance = 1e-12)
})

test_that("bp_estimate_disp_dense matches edgeR estimateDisp with group and library sizes", {
    fixture <- make_counts(ngenes = 90, nsamples = 6, seed = 15)
    bp <- make_bp_matrix(fixture$counts)
    keep <- bp_filter_by_expr(bp, group = fixture$group, block.size = 11)
    bp_filtered <- bp[keep, , drop = FALSE]
    dense <- fixture$counts[keep, , drop = FALSE]
    lib.size <- bp_library_sizes(bp_filtered)

    observed <- bp_estimate_disp_dense(
        bp_filtered,
        group = fixture$group,
        lib.size = lib.size,
        trend.method = "none"
    )
    expected <- edgeR::estimateDisp(
        dense,
        group = fixture$group,
        lib.size = lib.size,
        trend.method = "none"
    )

    expect_equal(observed$common.dispersion, expected$common.dispersion, tolerance = 1e-12)
    expect_equal(observed$tagwise.dispersion, expected$tagwise.dispersion, tolerance = 1e-12)
    expect_equal(observed$prior.df, expected$prior.df, tolerance = 1e-12)
    expect_equal(observed$prior.n, expected$prior.n, tolerance = 1e-12)
})

test_that("bp_estimate_disp_dense supports tagwise = FALSE", {
    fixture <- make_counts(ngenes = 100, nsamples = 6, seed = 16)
    bp <- make_bp_matrix(fixture$counts)
    prep <- bp_prepare_dge(bp, design = fixture$design, block.size = 12)
    dense <- fixture$counts[prep$keep, , drop = FALSE]

    observed <- bp_estimate_disp_dense(
        prep$counts,
        design = fixture$design,
        offset = prep$offset,
        tagwise = FALSE
    )
    expected <- edgeR::estimateDisp(
        dense,
        design = fixture$design,
        offset = prep$offset,
        tagwise = FALSE
    )

    expect_equal(observed$common.dispersion, expected$common.dispersion, tolerance = 1e-12)
    expect_equal(observed$trended.dispersion, expected$trended.dispersion, tolerance = 1e-12)
    expect_null(observed$tagwise.dispersion)
})

test_that("bp_estimate_disp_dense protects against accidental large dense conversion", {
    fixture <- make_counts(ngenes = 20, nsamples = 4, seed = 17)
    bp <- make_bp_matrix(fixture$counts)

    expect_error(
        bp_estimate_disp_dense(
            bp,
            design = fixture$design[seq_len(4), , drop = FALSE],
            max.elements = 10
        ),
        "dense dispersion estimation is only a reference/small-data path"
    )
})

test_that("bp_estimate_disp_dense rejects gene-specific offset matrices for now", {
    fixture <- make_counts(ngenes = 60, nsamples = 6, seed = 18)
    bp <- make_bp_matrix(fixture$counts)
    prep <- bp_prepare_dge(bp, design = fixture$design, block.size = 10)
    offset <- matrix(
        prep$offset,
        nrow = nrow(prep$counts),
        ncol = ncol(prep$counts),
        byrow = TRUE
    )
    offset[1, 1] <- offset[1, 1] + 0.01

    expect_error(
        bp_estimate_disp_dense(
            prep$counts,
            design = fixture$design,
            offset = offset
        ),
        "Gene-specific offset matrices are not yet supported"
    )
})
