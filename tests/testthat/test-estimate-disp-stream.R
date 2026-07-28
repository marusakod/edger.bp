test_that("bp_estimate_disp_stream matches edgeR estimateDisp with prepared offsets", {
    fixture <- make_counts(ngenes = 120, nsamples = 6, seed = 24)
    bp <- make_bp_matrix(fixture$counts)
    prep <- bp_prepare_dge(bp, design = fixture$design, block.size = 17)
    dense <- fixture$counts[prep$keep, , drop = FALSE]

    observed <- bp_estimate_disp_stream(
        prep$counts,
        design = fixture$design,
        offset = prep$offset,
        block.size = 13
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

test_that("bp_estimate_disp_stream matches edgeR with explicit library sizes", {
    fixture <- make_counts(ngenes = 100, nsamples = 6, seed = 25)
    bp <- make_bp_matrix(fixture$counts)
    prep <- bp_prepare_dge(bp, design = fixture$design, block.size = 12)
    dense <- fixture$counts[prep$keep, , drop = FALSE]

    observed <- bp_estimate_disp_stream(
        prep$counts,
        design = fixture$design,
        lib.size = prep$effective.lib.size,
        offset = prep$offset,
        block.size = 9
    )
    expected <- edgeR::estimateDisp(
        dense,
        design = fixture$design,
        lib.size = prep$effective.lib.size,
        offset = prep$offset
    )

    expect_equal(observed$common.dispersion, expected$common.dispersion, tolerance = 1e-12)
    expect_equal(observed$trended.dispersion, expected$trended.dispersion, tolerance = 1e-12)
    expect_equal(observed$tagwise.dispersion, expected$tagwise.dispersion, tolerance = 1e-12)
    expect_equal(observed$prior.df, expected$prior.df, tolerance = 1e-12)
    expect_equal(observed$prior.n, expected$prior.n, tolerance = 1e-12)
})

test_that("bp_estimate_disp_stream supports tagwise = FALSE", {
    fixture <- make_counts(ngenes = 90, nsamples = 6, seed = 26)
    bp <- make_bp_matrix(fixture$counts)
    prep <- bp_prepare_dge(bp, design = fixture$design, block.size = 10)
    dense <- fixture$counts[prep$keep, , drop = FALSE]

    observed <- bp_estimate_disp_stream(
        prep$counts,
        design = fixture$design,
        offset = prep$offset,
        tagwise = FALSE,
        block.size = 11
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

test_that("bp_estimate_disp_stream supports trend.method = none when all rows are selected", {
    fixture <- make_counts(ngenes = 70, nsamples = 6, seed = 27)
    bp <- make_bp_matrix(fixture$counts)
    prep <- bp_prepare_dge(bp, design = fixture$design, filter = FALSE, block.size = 8)

    observed <- bp_estimate_disp_stream(
        prep$counts,
        design = fixture$design,
        lib.size = prep$effective.lib.size,
        offset = prep$offset,
        trend.method = "none",
        min.row.sum = 0,
        block.size = 7
    )
    expected <- edgeR::estimateDisp(
        fixture$counts,
        design = fixture$design,
        lib.size = prep$effective.lib.size,
        offset = prep$offset,
        trend.method = "none",
        min.row.sum = 0
    )

    expect_equal(observed$common.dispersion, expected$common.dispersion, tolerance = 1e-12)
    expect_null(observed$trended.dispersion)
    expect_equal(observed$tagwise.dispersion, expected$tagwise.dispersion, tolerance = 1e-12)
    expect_equal(observed$prior.df, expected$prior.df, tolerance = 1e-12)
    expect_equal(observed$prior.n, expected$prior.n, tolerance = 1e-12)
})

test_that("bp_estimate_disp_stream reports phase progress when verbose", {
    fixture <- make_counts(ngenes = 70, nsamples = 6, seed = 29)
    bp <- make_bp_matrix(fixture$counts)
    prep <- bp_prepare_dge(bp, design = fixture$design, block.size = 8)

    messages <- capture.output(
        observed <- bp_estimate_disp_stream(
            prep$counts,
            design = fixture$design,
            offset = prep$offset,
            block.size = 7,
            verbose = TRUE
        ),
        type = "message"
    )

    expect_s3_class(observed, "BPCellsStreamDisp")
    expect_true(any(grepl("streamed dispersion: building APL grid", messages, fixed = TRUE)))
    expect_true(any(grepl("streamed dispersion: common dispersion", messages, fixed = TRUE)))
    expect_true(any(grepl("streamed dispersion: estimating tagwise dispersions", messages, fixed = TRUE)))
    expect_true(any(grepl("streamed dispersion: finished", messages, fixed = TRUE)))
})

test_that("bp_estimate_disp_stream rejects unsupported options", {
    fixture <- make_counts(ngenes = 40, nsamples = 6, seed = 28)
    bp <- make_bp_matrix(fixture$counts)

    expect_error(
        bp_estimate_disp_stream(
            bp,
            design = fixture$design,
            robust = TRUE
        ),
        "`robust = TRUE` is not yet supported"
    )
    expect_error(
        bp_estimate_disp_stream(
            bp,
            design = fixture$design,
            weights = matrix(1, nrow = nrow(fixture$counts), ncol = ncol(fixture$counts))
        ),
        "Observation `weights` are not yet supported"
    )
})
