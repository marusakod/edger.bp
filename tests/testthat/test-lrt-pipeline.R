test_that("bp_lrt matches canonical dense DGEList LRT workflow", {
    fixture <- make_counts(ngenes = 120, nsamples = 6, seed = 20)
    bp <- make_bp_matrix(fixture$counts)

    observed <- bp_lrt(
        bp,
        design = fixture$design,
        coef = 2,
        block.size = 13
    )

    dge <- edgeR::DGEList(fixture$counts)
    keep <- edgeR::filterByExpr(dge, design = fixture$design)
    dge <- dge[keep, , keep.lib.sizes = FALSE]
    dge <- edgeR::normLibSizes(dge)
    dge <- edgeR::estimateDisp(dge, design = fixture$design)
    fit <- edgeR::glmFit(dge, design = fixture$design)
    expected <- edgeR::glmLRT(fit, coef = 2)

    expect_identical(observed$prepare$keep, keep)
    expect_equal(
        observed$dispersion$common.dispersion,
        dge$common.dispersion,
        tolerance = 1e-12
    )
    expect_equal(
        observed$dispersion$trended.dispersion,
        dge$trended.dispersion,
        tolerance = 1e-12
    )
    expect_equal(
        observed$dispersion$tagwise.dispersion,
        dge$tagwise.dispersion,
        tolerance = 1e-12
    )
    expect_equal(observed$lrt$table, expected$table, tolerance = 1e-12)
    expect_equal(observed$lrt$coefficients, expected$coefficients, tolerance = 1e-12)
    expect_equal(observed$lrt$deviance, expected$deviance, tolerance = 1e-12, ignore_attr = TRUE)
    expect_equal(observed$lrt$AveLogCPM, dge$AveLogCPM, tolerance = 1e-12, ignore_attr = TRUE)
})

test_that("bp_lrt supports supplied dispersion without dense estimation", {
    fixture <- make_counts(ngenes = 100, nsamples = 6, seed = 21)
    bp <- make_bp_matrix(fixture$counts)
    dispersion <- 0.09

    observed <- bp_lrt(
        bp,
        design = fixture$design,
        dispersion = dispersion,
        coef = 2,
        block.size = 10
    )

    prep <- bp_prepare_dge(bp, design = fixture$design, block.size = 10)
    expected <- bp_glm_lrt(
        prep$counts,
        design = fixture$design,
        dispersion = dispersion,
        offset = prep$offset,
        coef = 2,
        block.size = 10
    )

    expect_null(observed$dispersion)
    expect_equal(observed$lrt$table, expected$table, tolerance = 1e-12)
    expect_equal(observed$lrt$coefficients, expected$coefficients, tolerance = 1e-12)
})

test_that("bp_lrt passes dispersion arguments to selected estimator", {
    fixture <- make_counts(ngenes = 90, nsamples = 6, seed = 22)
    bp <- make_bp_matrix(fixture$counts)

    observed <- bp_lrt(
        bp,
        design = fixture$design,
        coef = 2,
        block.size = 9,
        dispersion.args = list(tagwise = FALSE)
    )

    expect_null(observed$dispersion$tagwise.dispersion)
    expect_equal(
        observed$lrt$dispersion,
        observed$dispersion$trended.dispersion,
        tolerance = 0,
        ignore_attr = TRUE
    )
})

test_that("bp_lrt streamed default matches dense edgeR with an unbalanced design", {
    fixture <- make_counts(ngenes = 130, nsamples = 8, seed = 33)
    group <- factor(c("A", "A", "A", "B", "B", "B", "B", "B"))
    design <- stats::model.matrix(~ group)
    bp <- make_bp_matrix(fixture$counts)

    observed <- bp_lrt(
        bp,
        design = design,
        group = group,
        coef = 2,
        block.size = 11
    )

    dge <- edgeR::DGEList(fixture$counts)
    keep <- edgeR::filterByExpr(dge, design = design, group = group)
    dge <- dge[keep, , keep.lib.sizes = FALSE]
    dge <- edgeR::normLibSizes(dge)
    dge <- edgeR::estimateDisp(dge, design = design)
    fit <- edgeR::glmFit(dge, design = design)
    expected <- edgeR::glmLRT(fit, coef = 2)

    expect_identical(observed$prepare$keep, keep)
    expect_equal(observed$dispersion$common.dispersion, dge$common.dispersion, tolerance = 1e-12)
    expect_equal(observed$dispersion$trended.dispersion, dge$trended.dispersion, tolerance = 1e-12)
    expect_equal(observed$dispersion$tagwise.dispersion, dge$tagwise.dispersion, tolerance = 1e-12)
    expect_equal(observed$lrt$table, expected$table, tolerance = 1e-12)
    expect_equal(observed$lrt$coefficients, expected$coefficients, tolerance = 1e-12)
    expect_equal(observed$lrt$deviance, expected$deviance, tolerance = 1e-12, ignore_attr = TRUE)
})

test_that("bp_lrt streamed default matches dense edgeR with contrast tests", {
    fixture <- make_counts(ngenes = 140, nsamples = 9, seed = 34)
    group <- factor(rep(c("A", "B", "C"), each = 3))
    design <- stats::model.matrix(~ 0 + group)
    colnames(design) <- levels(group)
    contrast <- c(A = -1, B = 1, C = 0)
    bp <- make_bp_matrix(fixture$counts)

    observed <- bp_lrt(
        bp,
        design = design,
        group = group,
        contrast = contrast,
        block.size = 10
    )

    dge <- edgeR::DGEList(fixture$counts)
    keep <- edgeR::filterByExpr(dge, design = design, group = group)
    dge <- dge[keep, , keep.lib.sizes = FALSE]
    dge <- edgeR::normLibSizes(dge)
    dge <- edgeR::estimateDisp(dge, design = design)
    fit <- edgeR::glmFit(dge, design = design)
    expected <- edgeR::glmLRT(fit, contrast = contrast)

    expect_identical(observed$prepare$keep, keep)
    expect_equal(observed$dispersion$tagwise.dispersion, dge$tagwise.dispersion, tolerance = 1e-12)
    expect_equal(observed$lrt$table, expected$table, tolerance = 1e-12)
    expect_equal(observed$lrt$coefficients, expected$coefficients, tolerance = 1e-12)
    expect_equal(observed$lrt$comparison, expected$comparison)
})

test_that("bp_lrt streamed default handles varied zero-fit patterns", {
    fixture <- make_counts(ngenes = 150, nsamples = 6, seed = 35)
    counts <- fixture$counts
    counts[4:8, 1:3] <- 0
    counts[9:13, 4:6] <- 0
    counts[14:18, c(1, 4)] <- 0
    counts[19:23, c(2, 5)] <- 0
    counts[24:28, c(3, 6)] <- 0
    bp <- make_bp_matrix(counts)

    observed <- bp_lrt(
        bp,
        design = fixture$design,
        group = fixture$group,
        coef = 2,
        block.size = 8
    )

    dge <- edgeR::DGEList(counts)
    keep <- edgeR::filterByExpr(dge, design = fixture$design, group = fixture$group)
    dge <- dge[keep, , keep.lib.sizes = FALSE]
    dge <- edgeR::normLibSizes(dge)
    dge <- edgeR::estimateDisp(dge, design = fixture$design)
    fit <- edgeR::glmFit(dge, design = fixture$design)
    expected <- edgeR::glmLRT(fit, coef = 2)

    expect_identical(observed$prepare$keep, keep)
    expect_equal(observed$dispersion$common.dispersion, dge$common.dispersion, tolerance = 1e-12)
    expect_equal(observed$dispersion$tagwise.dispersion, dge$tagwise.dispersion, tolerance = 1e-12)
    expect_equal(observed$lrt$table, expected$table, tolerance = 1e-12)
    expect_equal(observed$lrt$AveLogCPM, dge$AveLogCPM, tolerance = 1e-12, ignore_attr = TRUE)
})

test_that("bp_lrt streamed results are invariant to block size", {
    fixture <- make_counts(ngenes = 120, nsamples = 6, seed = 36)
    bp <- make_bp_matrix(fixture$counts)

    small.blocks <- bp_lrt(
        bp,
        design = fixture$design,
        group = fixture$group,
        coef = 2,
        block.size = 5
    )
    large.blocks <- bp_lrt(
        bp,
        design = fixture$design,
        group = fixture$group,
        coef = 2,
        block.size = 29
    )

    expect_equal(small.blocks$dispersion$common.dispersion, large.blocks$dispersion$common.dispersion, tolerance = 1e-12)
    expect_equal(small.blocks$dispersion$tagwise.dispersion, large.blocks$dispersion$tagwise.dispersion, tolerance = 1e-12)
    expect_equal(small.blocks$lrt$table, large.blocks$lrt$table, tolerance = 1e-12)
    expect_equal(small.blocks$lrt$coefficients, large.blocks$lrt$coefficients, tolerance = 1e-12)
})

test_that("bp_lrt stream and dense dispersion methods agree on oracle data", {
    fixture <- make_counts(ngenes = 100, nsamples = 6, seed = 37)
    bp <- make_bp_matrix(fixture$counts)

    stream <- bp_lrt(
        bp,
        design = fixture$design,
        group = fixture$group,
        coef = 2,
        block.size = 9,
        dispersion.method = "stream"
    )
    dense <- bp_lrt(
        bp,
        design = fixture$design,
        group = fixture$group,
        coef = 2,
        block.size = 9,
        dispersion.method = "dense"
    )

    expect_equal(stream$prepare$keep, dense$prepare$keep, tolerance = 0)
    expect_equal(stream$dispersion$common.dispersion, dense$dispersion$common.dispersion, tolerance = 1e-12)
    expect_equal(stream$dispersion$trended.dispersion, dense$dispersion$trended.dispersion, tolerance = 1e-12)
    expect_equal(stream$dispersion$tagwise.dispersion, dense$dispersion$tagwise.dispersion, tolerance = 1e-12)
    expect_equal(stream$lrt$table, dense$lrt$table, tolerance = 1e-12)
    expect_equal(stream$lrt$coefficients, dense$lrt$coefficients, tolerance = 1e-12)
})

test_that("bp_lrt propagates dense estimation guard", {
    fixture <- make_counts(ngenes = 20, nsamples = 4, seed = 23)
    bp <- make_bp_matrix(fixture$counts)
    design <- fixture$design[seq_len(4), , drop = FALSE]

    expect_error(
        bp_lrt(
            bp,
            design = design,
            dispersion.method = "dense",
            max.elements = 10
        ),
        "dense dispersion estimation is only a reference/small-data path"
    )
})
