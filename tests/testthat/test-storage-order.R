test_that("bp_ensure_storage_order writes row-order companion when needed", {
    fixture <- make_counts(ngenes = 25, nsamples = 5, seed = 30)
    bp <- make_bp_matrix(fixture$counts)
    outdir <- tempfile("bp_row_", tmpdir = tempdir())

    observed <- bp_ensure_storage_order(bp, order = "row", outdir = outdir)

    expect_identical(BPCells::storage_order(bp), "col")
    expect_identical(BPCells::storage_order(observed), "row")
    expect_true(dir.exists(outdir))
    expect_equal(as.matrix(observed), fixture$counts, tolerance = 0)
})

test_that("bp_ensure_storage_order returns existing matrix when order already matches", {
    fixture <- make_counts(ngenes = 20, nsamples = 4, seed = 31)
    bp <- make_bp_matrix(fixture$counts)

    observed <- bp_ensure_storage_order(bp, order = "col")

    expect_identical(observed, bp)
})

test_that("bp_ensure_storage_order requires a fresh output directory", {
    fixture <- make_counts(ngenes = 20, nsamples = 4, seed = 32)
    bp <- make_bp_matrix(fixture$counts)

    expect_error(
        bp_ensure_storage_order(bp, order = "row"),
        "`outdir` is required"
    )
    expect_error(
        bp_ensure_storage_order(bp, order = "row", outdir = tempdir()),
        "`outdir` already exists"
    )
})
