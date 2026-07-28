make_bp_matrix <- function(counts) {
    sparse <- Matrix::Matrix(counts, sparse = TRUE)
    suppressWarnings(BPCells::write_matrix_dir(
        sparse,
        tempfile("bp_", tmpdir = tempdir())
    ))
}

make_counts <- function(ngenes = 80, nsamples = 6, seed = 1) {
    set.seed(seed)
    group <- factor(rep(c("A", "B"), length.out = nsamples))
    base.mean <- rgamma(ngenes, shape = 2, rate = 0.2)
    lib.scale <- seq(0.7, 1.4, length.out = nsamples)
    counts <- vapply(seq_len(nsamples), function(j) {
        stats::rnbinom(ngenes, mu = base.mean * lib.scale[j], size = 4)
    }, numeric(ngenes))
    counts[sample(length(counts), size = floor(length(counts) * 0.35))] <- 0
    counts[1:3, ] <- 0
    rownames(counts) <- paste0("gene", seq_len(ngenes))
    colnames(counts) <- paste0("sample", seq_len(nsamples))
    storage.mode(counts) <- "numeric"
    list(
        counts = counts,
        group = group,
        design = stats::model.matrix(~ group)
    )
}
