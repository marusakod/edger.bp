test_that("bp_prepare_dge matches the dense edgeR filtering and normalization workflow", {
    fixture <- make_counts(seed = 7)
    bp <- make_bp_matrix(fixture$counts)

    observed <- bp_prepare_dge(
        bp,
        design = fixture$design,
        block.size = 8
    )

    dge <- edgeR::DGEList(fixture$counts)
    keep <- edgeR::filterByExpr(dge, design = fixture$design)
    dge <- dge[keep, , keep.lib.sizes = FALSE]
    dge <- edgeR::normLibSizes(dge, method = "TMM")

    expect_identical(observed$keep, keep)
    expect_equal(unname(observed$lib.size), dge$samples$lib.size, tolerance = 0)
    expect_equal(unname(observed$norm.factors), dge$samples$norm.factors, tolerance = 1e-12)
    expect_equal(
        unname(observed$effective.lib.size),
        dge$samples$lib.size * dge$samples$norm.factors,
        tolerance = 1e-12
    )
    expect_equal(
        unname(observed$offset),
        log(dge$samples$lib.size * dge$samples$norm.factors),
        tolerance = 1e-12
    )
    expect_equal(dim(observed$counts), dim(dge$counts))
    expect_identical(rownames(observed$counts), rownames(dge$counts))
    expect_identical(colnames(observed$counts), colnames(dge$counts))
})

test_that("bp_prepare_dge can skip filtering", {
    fixture <- make_counts(seed = 8)
    bp <- make_bp_matrix(fixture$counts)

    observed <- bp_prepare_dge(bp, filter = FALSE)
    expected <- edgeR::normLibSizes(fixture$counts, method = "TMM")

    expect_true(all(observed$keep))
    expect_equal(unname(observed$lib.size), as.numeric(colSums(fixture$counts)), tolerance = 0)
    expect_equal(unname(observed$norm.factors), unname(expected), tolerance = 1e-12)
})

test_that("bp_prepare_dge subsets row-ordered BPCells matrices by position", {
    fixture <- make_counts(seed = 81)
    bp_col <- make_bp_matrix(fixture$counts)
    bp_row <- bp_ensure_storage_order(
        bp_col,
        order = "row",
        outdir = tempfile("bp_row_", tmpdir = tempdir())
    )

    observed <- bp_prepare_dge(
        bp_row,
        design = fixture$design,
        block.size = 8
    )

    dge <- edgeR::DGEList(fixture$counts)
    keep <- edgeR::filterByExpr(dge, design = fixture$design)
    dge <- dge[keep, , keep.lib.sizes = FALSE]
    dge <- edgeR::normLibSizes(dge, method = "TMM")

    expect_identical(observed$keep, keep)
    expect_equal(unname(observed$lib.size), dge$samples$lib.size, tolerance = 0)
    expect_identical(rownames(observed$counts), rownames(dge$counts))
})
