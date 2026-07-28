test_that("bp_qlf matches dense edgeR legacy QLF workflow", {
    fixture <- make_counts(ngenes = 120, nsamples = 6, seed = 44)
    bp <- make_bp_matrix(fixture$counts)

    observed <- bp_qlf(
        bp,
        design = fixture$design,
        group = fixture$group,
        coef = 2,
        block.size = 11,
        qlf.args = list(legacy = TRUE)
    )

    dge <- edgeR::DGEList(fixture$counts)
    keep <- edgeR::filterByExpr(dge, design = fixture$design, group = fixture$group)
    dge <- dge[keep, , keep.lib.sizes = FALSE]
    dge <- edgeR::normLibSizes(dge)
    dge <- edgeR::estimateDisp(dge, design = fixture$design)
    fit <- edgeR::glmQLFit(
        dge,
        design = fixture$design,
        dispersion = dge$tagwise.dispersion,
        legacy = TRUE
    )
    expected <- edgeR::glmQLFTest(fit, coef = 2)

    expect_identical(observed$prepare$keep, keep)
    expect_equal(observed$dispersion$tagwise.dispersion, dge$tagwise.dispersion, tolerance = 1e-12)
    expect_equal(observed$qlf$table, expected$table, tolerance = 1e-12)
    expect_equal(observed$qlf$coefficients, expected$coefficients, tolerance = 1e-12)
    expect_equal(observed$qlf$s2.post, fit$s2.post, tolerance = 1e-12, ignore_attr = TRUE)
})

test_that("bp_qlf supports supplied dispersion without estimating it", {
    fixture <- make_counts(ngenes = 90, nsamples = 6, seed = 45)
    bp <- make_bp_matrix(fixture$counts)
    dispersion <- 0.12

    observed <- bp_qlf(
        bp,
        design = fixture$design,
        group = fixture$group,
        dispersion = dispersion,
        coef = 2,
        block.size = 9
    )

    prep <- bp_prepare_dge(bp, design = fixture$design, group = fixture$group, block.size = 9)
    expected <- bp_glm_qlf(
        prep$counts,
        design = fixture$design,
        dispersion = dispersion,
        offset = prep$offset,
        coef = 2,
        block.size = 9,
        legacy = FALSE
    )

    expect_null(observed$dispersion)
    expect_equal(observed$qlf$table, expected$table, tolerance = 1e-12)
    expect_equal(observed$qlf$coefficients, expected$coefficients, tolerance = 1e-12)
})

test_that("bp_qlf defaults to dense edgeR modern QLF workflow", {
    fixture <- make_counts(ngenes = 115, nsamples = 6, seed = 48)
    bp <- make_bp_matrix(fixture$counts)

    observed <- bp_qlf(
        bp,
        design = fixture$design,
        group = fixture$group,
        coef = 2,
        block.size = 10
    )

    dge <- edgeR::DGEList(fixture$counts)
    keep <- edgeR::filterByExpr(dge, design = fixture$design, group = fixture$group)
    dge <- dge[keep, , keep.lib.sizes = FALSE]
    dge <- edgeR::normLibSizes(dge)
    dge <- edgeR::estimateDisp(dge, design = fixture$design)
    fit <- edgeR::glmQLFit(
        dge,
        design = fixture$design,
        dispersion = dge$tagwise.dispersion,
        legacy = FALSE
    )
    expected <- edgeR::glmQLFTest(fit, coef = 2)

    expect_identical(observed$prepare$keep, keep)
    expect_equal(observed$dispersion$tagwise.dispersion, dge$tagwise.dispersion, tolerance = 1e-12)
    expect_equal(observed$qlf$table, expected$table, tolerance = 1e-10)
    expect_equal(observed$qlf$coefficients, expected$coefficients, tolerance = 1e-10)
    expect_equal(observed$qlf$average.ql.dispersion, fit$average.ql.dispersion, tolerance = 1e-10)
    expect_equal(observed$qlf$s2.post, fit$s2.post, tolerance = 1e-10, ignore_attr = TRUE)
})
