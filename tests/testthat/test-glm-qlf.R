test_that("bp_glm_qlf matches edgeR legacy glmQLFit and glmQLFTest", {
    fixture <- make_counts(ngenes = 100, nsamples = 6, seed = 41)
    bp <- make_bp_matrix(fixture$counts)
    prep <- bp_prepare_dge(bp, design = fixture$design, group = fixture$group, block.size = 10)
    dense <- fixture$counts[prep$keep, , drop = FALSE]
    disp <- edgeR::estimateDisp(dense, design = fixture$design, offset = prep$offset)
    ave <- edgeR::aveLogCPM(
        dense,
        lib.size = prep$effective.lib.size,
        dispersion = disp$common.dispersion
    )

    observed <- bp_glm_qlf(
        prep$counts,
        design = fixture$design,
        dispersion = disp$tagwise.dispersion,
        offset = prep$offset,
        coef = 2,
        ave.log.cpm = ave,
        block.size = 7,
        legacy = TRUE
    )
    fit <- edgeR::glmQLFit(
        dense,
        design = fixture$design,
        dispersion = disp$tagwise.dispersion,
        offset = prep$offset,
        AveLogCPM = ave,
        legacy = TRUE
    )
    expected <- edgeR::glmQLFTest(fit, coef = 2)

    expect_equal(observed$table, expected$table, tolerance = 1e-12)
    expect_equal(observed$coefficients, expected$coefficients, tolerance = 1e-12)
    expect_equal(observed$unshrunk.coefficients, expected$unshrunk.coefficients, tolerance = 1e-12)
    expect_equal(observed$df.total, expected$df.total, tolerance = 1e-12, ignore_attr = TRUE)
    expect_equal(observed$s2.post, fit$s2.post, tolerance = 1e-12, ignore_attr = TRUE)
    expect_equal(observed$df.residual.zeros, fit$df.residual.zeros, tolerance = 0, ignore_attr = TRUE)
})

test_that("bp_glm_qlf matches edgeR legacy contrast tests", {
    fixture <- make_counts(ngenes = 110, nsamples = 9, seed = 42)
    group <- factor(rep(c("A", "B", "C"), each = 3))
    design <- stats::model.matrix(~ 0 + group)
    colnames(design) <- levels(group)
    contrast <- c(A = -1, B = 1, C = 0)
    bp <- make_bp_matrix(fixture$counts)
    prep <- bp_prepare_dge(bp, design = design, group = group, block.size = 9)
    dense <- fixture$counts[prep$keep, , drop = FALSE]
    disp <- edgeR::estimateDisp(dense, design = design, offset = prep$offset)
    ave <- edgeR::aveLogCPM(
        dense,
        lib.size = prep$effective.lib.size,
        dispersion = disp$common.dispersion
    )

    observed <- bp_glm_qlf(
        prep$counts,
        design = design,
        dispersion = disp$tagwise.dispersion,
        offset = prep$offset,
        contrast = contrast,
        ave.log.cpm = ave,
        block.size = 8,
        legacy = TRUE
    )
    fit <- edgeR::glmQLFit(
        dense,
        design = design,
        dispersion = disp$tagwise.dispersion,
        offset = prep$offset,
        AveLogCPM = ave,
        legacy = TRUE
    )
    expected <- edgeR::glmQLFTest(fit, contrast = contrast)

    expect_equal(observed$table, expected$table, tolerance = 1e-12)
    expect_equal(observed$coefficients, expected$coefficients, tolerance = 1e-12)
    expect_identical(observed$comparison, expected$comparison)
})

test_that("bp_glm_qlf matches edgeR modern glmQLFit and glmQLFTest", {
    fixture <- make_counts(ngenes = 120, nsamples = 8, seed = 46)
    bp <- make_bp_matrix(fixture$counts)
    prep <- bp_prepare_dge(bp, design = fixture$design, group = fixture$group, block.size = 10)
    dense <- fixture$counts[prep$keep, , drop = FALSE]
    disp <- edgeR::estimateDisp(dense, design = fixture$design, offset = prep$offset)
    ave <- edgeR::aveLogCPM(
        dense,
        lib.size = prep$effective.lib.size,
        dispersion = disp$common.dispersion
    )

    observed <- bp_glm_qlf(
        prep$counts,
        design = fixture$design,
        dispersion = disp$tagwise.dispersion,
        offset = prep$offset,
        coef = 2,
        ave.log.cpm = ave,
        block.size = 9,
        legacy = FALSE
    )
    fit <- edgeR::glmQLFit(
        dense,
        design = fixture$design,
        dispersion = disp$tagwise.dispersion,
        offset = prep$offset,
        AveLogCPM = ave,
        legacy = FALSE
    )
    fit$AveLogCPM <- ave
    expected <- edgeR::glmQLFTest(fit, coef = 2)

    expect_equal(observed$table, expected$table, tolerance = 1e-10)
    expect_equal(observed$coefficients, expected$coefficients, tolerance = 1e-10)
    expect_equal(observed$unshrunk.coefficients, expected$unshrunk.coefficients, tolerance = 1e-10)
    expect_equal(observed$average.ql.dispersion, fit$average.ql.dispersion, tolerance = 1e-10)
    expect_equal(observed$df.residual.adj, fit$df.residual.adj, tolerance = 1e-10, ignore_attr = TRUE)
    expect_equal(observed$deviance.adj, fit$deviance.adj, tolerance = 1e-10, ignore_attr = TRUE)
    expect_equal(observed$df.total, expected$df.total, tolerance = 1e-10, ignore_attr = TRUE)
    expect_equal(observed$s2.post, fit$s2.post, tolerance = 1e-10, ignore_attr = TRUE)
})

test_that("bp_glm_qlf matches edgeR modern contrast tests", {
    fixture <- make_counts(ngenes = 125, nsamples = 9, seed = 47)
    group <- factor(rep(c("A", "B", "C"), each = 3))
    design <- stats::model.matrix(~ 0 + group)
    colnames(design) <- levels(group)
    contrast <- c(A = -1, B = 1, C = 0)
    bp <- make_bp_matrix(fixture$counts)
    prep <- bp_prepare_dge(bp, design = design, group = group, block.size = 11)
    dense <- fixture$counts[prep$keep, , drop = FALSE]
    disp <- edgeR::estimateDisp(dense, design = design, offset = prep$offset)
    ave <- edgeR::aveLogCPM(
        dense,
        lib.size = prep$effective.lib.size,
        dispersion = disp$common.dispersion
    )

    observed <- bp_glm_qlf(
        prep$counts,
        design = design,
        dispersion = disp$tagwise.dispersion,
        offset = prep$offset,
        contrast = contrast,
        ave.log.cpm = ave,
        block.size = 10,
        legacy = FALSE
    )
    fit <- edgeR::glmQLFit(
        dense,
        design = design,
        dispersion = disp$tagwise.dispersion,
        offset = prep$offset,
        AveLogCPM = ave,
        legacy = FALSE
    )
    fit$AveLogCPM <- ave
    expected <- edgeR::glmQLFTest(fit, contrast = contrast)

    expect_equal(observed$table, expected$table, tolerance = 1e-10)
    expect_equal(observed$coefficients, expected$coefficients, tolerance = 1e-10)
    expect_identical(observed$comparison, expected$comparison)
})

test_that("bp_glm_qlf modern pass fusion supports multi-coefficient tests", {
    fixture <- make_counts(ngenes = 130, nsamples = 9, seed = 470)
    group <- factor(rep(c("A", "B", "C"), each = 3))
    design <- stats::model.matrix(~ group)
    bp <- make_bp_matrix(fixture$counts)
    prep <- bp_prepare_dge(bp, design = design, group = group, block.size = 11)
    dense <- fixture$counts[prep$keep, , drop = FALSE]
    disp <- edgeR::estimateDisp(dense, design = design, offset = prep$offset)
    ave <- edgeR::aveLogCPM(
        dense,
        lib.size = prep$effective.lib.size,
        dispersion = disp$common.dispersion
    )

    observed <- bp_glm_qlf(
        prep$counts,
        design = design,
        dispersion = disp$tagwise.dispersion,
        offset = prep$offset,
        coef = 2:3,
        ave.log.cpm = ave,
        block.size = 10,
        legacy = FALSE,
        keep.fitted = TRUE
    )
    fit <- edgeR::glmQLFit(
        dense,
        design = design,
        dispersion = disp$tagwise.dispersion,
        offset = prep$offset,
        AveLogCPM = ave,
        legacy = FALSE
    )
    fit$AveLogCPM <- ave
    expected <- edgeR::glmQLFTest(fit, coef = 2:3)

    expect_equal(observed$table, expected$table, tolerance = 1e-10)
    expect_equal(observed$coefficients, expected$coefficients, tolerance = 1e-10)
    expect_equal(observed$fitted.values, fit$fitted.values, tolerance = 1e-10)
    expect_equal(observed$df.test, expected$df.test, tolerance = 1e-10, ignore_attr = TRUE)
    expect_equal(observed$df.total, expected$df.total, tolerance = 1e-10, ignore_attr = TRUE)
    expect_identical(observed$comparison, expected$comparison)
})

test_that("bp_glm_qlf modern fit reuse supports general designs and library sizes", {
    fixture <- make_counts(ngenes = 105, nsamples = 10, seed = 471)
    group <- rep(c(0, 1), each = 5)
    covariate <- seq(-1, 1, length.out = 10)
    design <- cbind(Intercept = 1, group = group, covariate = covariate)
    lib.size <- colSums(fixture$counts)
    ave <- edgeR::aveLogCPM(
        fixture$counts,
        lib.size = lib.size,
        dispersion = 0.12
    )
    bp <- make_bp_matrix(fixture$counts)

    observed <- bp_glm_qlf(
        bp,
        design = design,
        dispersion = 0.12,
        lib.size = lib.size,
        coef = 2:3,
        ave.log.cpm = ave,
        block.size = 9,
        legacy = FALSE
    )
    fit <- edgeR::glmQLFit(
        fixture$counts,
        design = design,
        dispersion = 0.12,
        lib.size = lib.size,
        AveLogCPM = ave,
        legacy = FALSE
    )
    fit$AveLogCPM <- ave
    expected <- edgeR::glmQLFTest(fit, coef = 2:3)

    expect_equal(observed$table, expected$table, tolerance = 1e-10)
    expect_equal(observed$coefficients, expected$coefficients, tolerance = 1e-10)
    expect_equal(observed$average.ql.dispersion, fit$average.ql.dispersion, tolerance = 1e-10)
    expect_equal(observed$s2.post, fit$s2.post, tolerance = 1e-10, ignore_attr = TRUE)
    expect_identical(observed$comparison, expected$comparison)
})

test_that("bp_glm_qlf rejects unsupported QL options", {
    fixture <- make_counts(ngenes = 40, nsamples = 6, seed = 43)
    bp <- make_bp_matrix(fixture$counts)

    expect_error(
        bp_glm_qlf(
            bp,
            design = fixture$design,
            dispersion = 0.1,
            lib.size = colSums(fixture$counts),
            robust = TRUE
        ),
        "`robust = TRUE` is not yet supported"
    )
})
