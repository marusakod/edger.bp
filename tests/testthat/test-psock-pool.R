test_that("PSOCK reuses workers, reports errors, and does not export whole frames", {
    path <- tempfile("psock-counts-")
    y <- matrix(seq_len(120), 20, dimnames = list(paste0("g", 1:20), paste0("s", 1:6)))
    BPCells::write_matrix_dir(Matrix::Matrix(y, sparse = TRUE), path)
    pool <- bp_start_psock_pool(path, workers = 2, verbose = FALSE)
    on.exit(bp_stop_psock_pool(pool), add = TRUE)
    counts <- pool$counts
    unrelated <- new.env() # Must never be included in a serialized callback.
    fn <- function(rows) base::rowSums(base::as.matrix(counts[rows, , drop = FALSE]))
    packet <- .psock_pack_function(fn, pool)
    expect_false("unrelated" %in% names(packet$bindings))
    expect_lt(length(serialize(packet, NULL)), 20000)
    expect_identical(bp_with_psock_pool(pool,
        .block_lapply(1:3, 2, function(i, increment) i + increment, increment = 2L)),
        as.list(3:5))
    ids <- bp_with_psock_pool(pool, bp_block_lapply(1:20, 2, function(i) Sys.getpid()))
    for (iteration in 1:5) {
        expect_identical(bp_with_psock_pool(pool,
            bp_block_lapply(1:20, 2, function(i) Sys.getpid())), ids)
        result <- bp_with_psock_pool(pool, bp_block_lapply(list(1:10, 11:20), 2, fn))
        expect_equal(unlist(result), rowSums(y))
    }
    failure <- tryCatch(bp_with_psock_pool(pool, bp_block_lapply(1:3, 2,
        function(i) { if (i == 2) stop("deliberate failure"); i }, context = "test pass")),
        error = identity)
    expect_s3_class(failure, "error")
    expect_match(conditionMessage(failure), "test pass", fixed = TRUE)
    expect_match(conditionMessage(failure), "deliberate failure", fixed = TRUE)
    expect_null(.bp_psock_state$pool)
    expect_error(bp_with_psock_pool(pool, bp_block_lapply(1:2, 1, identity)), "mismatch")
    bp_stop_psock_pool(pool)
    bp_stop_psock_pool(pool)
    expect_true(pool$closed)
    expect_error(bp_with_psock_pool(pool, 1), "closed")
})

test_that("PSOCK numerical results match serial dispersion, LRT, QLF and Treat", {
    set.seed(41)
    path <- tempfile("psock-fit-")
    raw <- matrix(rnbinom(180 * 12, mu = 30, size = 5), 180,
        dimnames = list(paste0("g", 1:180), paste0("s", 1:12)))
    BPCells::write_matrix_dir(Matrix::Matrix(raw, sparse = TRUE), path)
    pool <- bp_start_psock_pool(path, rows = 180:1, cols = 12:1, workers = 2, verbose = FALSE)
    on.exit(bp_stop_psock_pool(pool), add = TRUE)
    y <- pool$counts
    group <- factor(rep(c("a", "b", "c"), each = 4))
    design <- model.matrix(~ group)
    offset <- log(as.numeric(BPCells::colSums(y)))
    serial <- bp_estimate_disp_stream(y, design, offset = offset, block.size = 30, workers = 1)
    socket <- bp_with_psock_pool(pool,
        bp_estimate_disp_stream(y, design, offset = offset, block.size = 30, workers = 2))
    expect_equal(socket$tagwise.dispersion, serial$tagwise.dispersion, tolerance = 1e-10)
    for (fit in list(bp_glm_lrt, bp_glm_qlf, bp_glm_treat)) {
        expected <- fit(y, design, dispersion = serial$tagwise.dispersion,
            offset = offset, block.size = 30, workers = 1)
        actual <- bp_with_psock_pool(pool, fit(y, design,
            dispersion = serial$tagwise.dispersion, offset = offset, block.size = 30, workers = 2))
        expect_equal(actual$table, expected$table, tolerance = 1e-10)
    }
})
