test_that("bp_glm_lrt matches edgeR with scalar dispersion and vector offset", {
    fixture <- make_counts(seed = 9)
    bp <- make_bp_matrix(fixture$counts)
    prep <- bp_prepare_dge(bp, design = fixture$design, block.size = 10)
    dense <- fixture$counts[prep$keep, , drop = FALSE]
    dispersion <- 0.12

    observed <- bp_glm_lrt(
        prep$counts,
        design = fixture$design,
        dispersion = dispersion,
        offset = prep$offset,
        coef = 2,
        block.size = 7,
        keep.fitted = TRUE
    )

    fit <- edgeR::glmFit(
        dense,
        design = fixture$design,
        dispersion = dispersion,
        offset = prep$offset
    )
    expected <- edgeR::glmLRT(fit, coef = 2)

    expect_equal(observed$table, expected$table, tolerance = 1e-12)
    expect_equal(observed$coefficients, expected$coefficients, tolerance = 1e-12)
    expect_equal(observed$unshrunk.coefficients, expected$unshrunk.coefficients, tolerance = 1e-12)
    expect_equal(observed$deviance, expected$deviance, tolerance = 1e-12, ignore_attr = TRUE)
    expect_equal(observed$df.residual, expected$df.residual, tolerance = 0, ignore_attr = TRUE)
    expect_equal(observed$dispersion, rep(dispersion, nrow(dense)), tolerance = 0, ignore_attr = TRUE)
    expect_equal(observed$fitted.values, expected$fitted.values, tolerance = 1e-12)
    expect_identical(observed$comparison, expected$comparison)
    expect_equal(observed$df.test, expected$df.test, tolerance = 0, ignore_attr = TRUE)
})

test_that("bp_glm_lrt matches edgeR with gene-level dispersion and row-constant matrix offset", {
    fixture <- make_counts(seed = 10)
    bp <- make_bp_matrix(fixture$counts)
    prep <- bp_prepare_dge(bp, design = fixture$design, block.size = 9)
    dense <- fixture$counts[prep$keep, , drop = FALSE]
    dispersion <- seq(0.04, 0.2, length.out = nrow(dense))
    offset <- matrix(
        prep$offset,
        nrow = nrow(dense),
        ncol = ncol(dense),
        byrow = TRUE,
        dimnames = dimnames(dense)
    )
    contrast <- c(0, 1)

    observed <- bp_glm_lrt(
        prep$counts,
        design = fixture$design,
        dispersion = dispersion,
        offset = offset,
        contrast = contrast,
        block.size = 5
    )

    fit <- edgeR::glmFit(
        dense,
        design = fixture$design,
        dispersion = dispersion,
        offset = prep$offset
    )
    expected <- edgeR::glmLRT(fit, contrast = contrast)

    expect_equal(observed$table, expected$table, tolerance = 1e-12)
    expect_equal(observed$coefficients, expected$coefficients, tolerance = 1e-12)
    expect_equal(observed$deviance, expected$deviance, tolerance = 1e-12, ignore_attr = TRUE)
    expect_equal(observed$df.residual, expected$df.residual, tolerance = 0, ignore_attr = TRUE)
    expect_equal(observed$dispersion, dispersion, tolerance = 0, ignore_attr = TRUE)
    expect_identical(observed$comparison, expected$comparison)
})

test_that("bp_glm_lrt rejects gene-specific offset matrices for now", {
    fixture <- make_counts(seed = 13)
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
        bp_glm_lrt(
            prep$counts,
            design = fixture$design,
            dispersion = 0.1,
            offset = offset
        ),
        "Gene-specific offset matrices are not yet supported"
    )
})

test_that("bp_glm_lrt can use supplied library sizes when offset is absent", {
    fixture <- make_counts(seed = 11)
    bp <- make_bp_matrix(fixture$counts)
    prep <- bp_prepare_dge(bp, design = fixture$design, block.size = 11)
    dense <- fixture$counts[prep$keep, , drop = FALSE]
    dispersion <- 0.08

    observed <- bp_glm_lrt(
        prep$counts,
        design = fixture$design,
        dispersion = dispersion,
        lib.size = prep$effective.lib.size,
        coef = 2,
        block.size = 6
    )

    fit <- edgeR::glmFit(
        dense,
        design = fixture$design,
        dispersion = dispersion,
        lib.size = prep$effective.lib.size
    )
    expected <- edgeR::glmLRT(fit, coef = 2)

    expect_equal(observed$table, expected$table, tolerance = 1e-12)
    expect_equal(observed$coefficients, expected$coefficients, tolerance = 1e-12)
})

test_that("bp_glm_lrt requires offset or global library sizes", {
    fixture <- make_counts(seed = 12)
    bp <- make_bp_matrix(fixture$counts)
    prep <- bp_prepare_dge(bp, design = fixture$design, block.size = 10)

    expect_error(
        bp_glm_lrt(
            prep$counts,
            design = fixture$design,
            dispersion = 0.1,
            block.size = 5
        ),
        "Supply `offset` or `lib.size`"
    )
})
