test_that("block workers are prescheduled instead of forked per task", {
    pids <- .block_lapply(seq_len(40L), workers = 2L, function(i) Sys.getpid())

    expect_lte(length(unique(unlist(pids))), 2L)
})

test_that("streamed filtering is invariant to worker count", {
    fixture <- make_counts(ngenes = 90, nsamples = 6, seed = 50)
    bp <- make_bp_matrix(fixture$counts)

    serial <- bp_filter_by_expr(
        bp,
        design = fixture$design,
        group = fixture$group,
        block.size = 9,
        workers = 1
    )
    parallel <- bp_filter_by_expr(
        bp,
        design = fixture$design,
        group = fixture$group,
        block.size = 9,
        workers = 2
    )

    expect_identical(parallel, serial)
})

test_that("streamed dispersion is invariant to worker count", {
    fixture <- make_counts(ngenes = 100, nsamples = 6, seed = 51)
    bp <- make_bp_matrix(fixture$counts)
    prep <- bp_prepare_dge(bp, design = fixture$design, group = fixture$group, block.size = 10)

    serial <- bp_estimate_disp_stream(
        prep$counts,
        design = fixture$design,
        offset = prep$offset,
        block.size = 8,
        workers = 1
    )
    parallel <- bp_estimate_disp_stream(
        prep$counts,
        design = fixture$design,
        offset = prep$offset,
        block.size = 8,
        workers = 2
    )

    expect_equal(parallel$common.dispersion, serial$common.dispersion, tolerance = 1e-12)
    expect_equal(parallel$trended.dispersion, serial$trended.dispersion, tolerance = 1e-12)
    expect_equal(parallel$tagwise.dispersion, serial$tagwise.dispersion, tolerance = 1e-12)
})

test_that("bp_lrt is invariant to worker count", {
    fixture <- make_counts(ngenes = 100, nsamples = 6, seed = 52)
    bp <- make_bp_matrix(fixture$counts)

    serial <- bp_lrt(
        bp,
        design = fixture$design,
        group = fixture$group,
        coef = 2,
        block.size = 10,
        workers = 1
    )
    parallel <- bp_lrt(
        bp,
        design = fixture$design,
        group = fixture$group,
        coef = 2,
        block.size = 10,
        workers = 2
    )

    expect_equal(parallel$dispersion$tagwise.dispersion, serial$dispersion$tagwise.dispersion, tolerance = 1e-12)
    expect_equal(parallel$lrt$table, serial$lrt$table, tolerance = 1e-12)
    expect_equal(parallel$lrt$coefficients, serial$lrt$coefficients, tolerance = 1e-12)
})

test_that("bp_qlf is invariant to worker count", {
    fixture <- make_counts(ngenes = 90, nsamples = 6, seed = 53)
    bp <- make_bp_matrix(fixture$counts)

    serial <- bp_qlf(
        bp,
        design = fixture$design,
        group = fixture$group,
        coef = 2,
        block.size = 9,
        workers = 1
    )
    parallel <- bp_qlf(
        bp,
        design = fixture$design,
        group = fixture$group,
        coef = 2,
        block.size = 9,
        workers = 2
    )

    expect_equal(parallel$dispersion$tagwise.dispersion, serial$dispersion$tagwise.dispersion, tolerance = 1e-12)
    expect_equal(parallel$qlf$table, serial$qlf$table, tolerance = 1e-12)
    expect_equal(parallel$qlf$s2.post, serial$qlf$s2.post, tolerance = 1e-12, ignore_attr = TRUE)
})
