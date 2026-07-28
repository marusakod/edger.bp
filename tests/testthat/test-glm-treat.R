test_that("bp_glm_treat matches edgeR legacy glmQLFit and glmTreat", {
    fixture <- make_counts(ngenes = 100, nsamples = 6, seed = 61)
    bp <- make_bp_matrix(fixture$counts)
    prep <- bp_prepare_dge(bp, design = fixture$design, group = fixture$group, block.size = 10)
    dense <- fixture$counts[prep$keep, , drop = FALSE]
    disp <- edgeR::estimateDisp(dense, design = fixture$design, offset = prep$offset)
    ave <- edgeR::aveLogCPM(
        dense,
        lib.size = prep$effective.lib.size,
        dispersion = disp$common.dispersion
    )

    observed <- bp_glm_treat(
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
    expected <- edgeR::glmTreat(fit, coef = 2)

    expect_equal(observed$table, expected$table, tolerance = 1e-12)
    expect_equal(observed$coefficients, expected$coefficients, tolerance = 1e-12)
    expect_equal(observed$unshrunk.coefficients, expected$unshrunk.coefficients, tolerance = 1e-12)
    expect_equal(observed$s2.post, fit$s2.post, tolerance = 1e-12, ignore_attr = TRUE)
    expect_equal(observed$df.residual.zeros, fit$df.residual.zeros, tolerance = 0, ignore_attr = TRUE)
    expect_identical(observed$comparison, expected$comparison)
})

test_that("bp_glm_treat matches edgeR legacy contrast tests", {
    fixture <- make_counts(ngenes = 110, nsamples = 9, seed = 62)
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

    observed <- bp_glm_treat(
        prep$counts,
        design = design,
        dispersion = disp$tagwise.dispersion,
        offset = prep$offset,
        contrast = contrast,
        ave.log.cpm = ave,
        block.size = 8,
        lfc = log2(1.5),
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
    expected <- edgeR::glmTreat(fit, contrast = contrast, lfc = log2(1.5))

    expect_equal(observed$table, expected$table, tolerance = 1e-12)
    expect_equal(observed$coefficients, expected$coefficients, tolerance = 1e-12)
    expect_identical(observed$comparison, expected$comparison)
})

test_that("bp_glm_treat matches edgeR modern glmQLFit and glmTreat", {
    fixture <- make_counts(ngenes = 100, nsamples = 6, seed = 67)
    bp <- make_bp_matrix(fixture$counts)
    prep <- bp_prepare_dge(bp, design = fixture$design, group = fixture$group, block.size = 10)
    dense <- fixture$counts[prep$keep, , drop = FALSE]
    disp <- edgeR::estimateDisp(dense, design = fixture$design, offset = prep$offset)
    ave <- edgeR::aveLogCPM(
        dense,
        lib.size = prep$effective.lib.size,
        dispersion = disp$common.dispersion
    )

    observed <- bp_glm_treat(
        prep$counts,
        design = fixture$design,
        dispersion = disp$tagwise.dispersion,
        offset = prep$offset,
        coef = 2,
        ave.log.cpm = ave,
        block.size = 7
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
    expected <- edgeR::glmTreat(fit, coef = 2)

    expect_equal(observed$table, expected$table, tolerance = 1e-10)
    expect_equal(observed$coefficients, expected$coefficients, tolerance = 1e-10)
    expect_equal(observed$unshrunk.coefficients, expected$unshrunk.coefficients, tolerance = 1e-10)
    expect_equal(observed$s2.post, fit$s2.post, tolerance = 1e-10, ignore_attr = TRUE)
    expect_equal(observed$df.residual.adj, fit$df.residual.adj, tolerance = 1e-10, ignore_attr = TRUE)
    expect_equal(observed$average.ql.dispersion, fit$average.ql.dispersion, tolerance = 1e-10)
    expect_identical(observed$comparison, expected$comparison)
})

test_that("bp_glm_treat matches edgeR modern contrast tests", {
    fixture <- make_counts(ngenes = 110, nsamples = 9, seed = 68)
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

    observed <- bp_glm_treat(
        prep$counts,
        design = design,
        dispersion = disp$tagwise.dispersion,
        offset = prep$offset,
        contrast = contrast,
        ave.log.cpm = ave,
        block.size = 8,
        lfc = log2(1.5)
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
    expected <- edgeR::glmTreat(fit, contrast = contrast, lfc = log2(1.5))

    expect_equal(observed$table, expected$table, tolerance = 1e-10)
    expect_equal(observed$coefficients, expected$coefficients, tolerance = 1e-10)
    expect_identical(observed$comparison, expected$comparison)
})

test_that("bp_glm_treat with zero lfc delegates to default QLF testing", {
    fixture <- make_counts(ngenes = 80, nsamples = 6, seed = 63)
    bp <- make_bp_matrix(fixture$counts)
    prep <- bp_prepare_dge(bp, design = fixture$design, group = fixture$group, block.size = 8)

    observed <- bp_glm_treat(
        prep$counts,
        design = fixture$design,
        dispersion = 0.1,
        offset = prep$offset,
        coef = 2,
        lfc = 0,
        block.size = 8
    )
    expected <- bp_glm_qlf(
        prep$counts,
        design = fixture$design,
        dispersion = 0.1,
        offset = prep$offset,
        coef = 2,
        block.size = 8,
        legacy = FALSE
    )

    expect_equal(observed$table, expected$table, tolerance = 1e-12)
    expect_s3_class(observed, "BPCellsDGEQLF")
})

test_that("bp_glm_treat rejects unsupported options", {
    fixture <- make_counts(ngenes = 40, nsamples = 6, seed = 64)
    bp <- make_bp_matrix(fixture$counts)

    expect_error(
        bp_glm_treat(
            bp,
            design = fixture$design,
            dispersion = 0.1,
            lib.size = colSums(fixture$counts),
            lfc = -1
        ),
        "`lfc` must be a non-negative numeric scalar"
    )
    expect_error(
        bp_glm_treat(
            bp,
            design = fixture$design,
            dispersion = 0.1,
            lib.size = colSums(fixture$counts),
            robust = TRUE
        ),
        "`robust = TRUE` is not yet supported"
    )
})

test_that("Treat p-values are clamped after numerical integration", {
    integratepnorm <- function(a, b) rep(-1e-12, length(a))

    observed <- edger.bp:::.treat_p_value(
        z.left = c(10, 10),
        z.right = c(20, 20),
        null = "interval",
        integratepnorm = integratepnorm
    )

    expect_equal(observed, c(0, 0))
})
