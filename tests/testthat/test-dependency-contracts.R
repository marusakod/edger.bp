test_that("required edgeR exports are available", {
    required <- c(
        "estimateDisp",
        "glmFit",
        "glmLRT",
        "glmTreat",
        "aveLogCPM",
        "makeCompressedMatrix",
        "designAsFactor"
    )

    expect_true(all(vapply(required, exists, logical(1), envir = asNamespace("edgeR"))))
})

test_that("required edgeR internals retain the interfaces used by edger.bp", {
    required.formals <- list(
        maximizeInterpolant = c("x", "y"),
        WLEB = c(
            "theta", "loglik", "prior.n", "covariate", "trend.method", "span",
            "overall", "trend", "individual", "m0", "m0.out"
        ),
        .comboGroups = "truths",
        adjustedProfileLik = c(
            "dispersion", "y", "design", "offset", "weights", "start", "get.coef"
        ),
        .residDF = c("zero", "design"),
        .compressDispersions = c("y", "dispersion"),
        .isBelowPoissonBound = "glmfit",
        .integratepnorm = c("a", "b")
    )

    for (name in names(required.formals)) {
        fun <- getFromNamespace(name, "edgeR")
        expect_true(
            all(required.formals[[name]] %in% names(formals(fun))),
            info = name
        )
    }

    native <- getFromNamespace(".cxx_compute_adj_vec", "edgeR")
    expect_s3_class(native, "NativeSymbolInfo")
})

test_that("required limma functions retain the interfaces used by edger.bp", {
    required.formals <- list(
        squeezeVar = c(
            "var", "df", "covariate", "robust", "winsor.tail.p", "legacy"
        ),
        zscoreT = c("x", "df"),
        contrastAsCoef = c("design", "contrast", "first")
    )

    for (name in names(required.formals)) {
        fun <- getFromNamespace(name, "limma")
        expect_true(
            all(required.formals[[name]] %in% names(formals(fun))),
            info = name
        )
    }
})

test_that("required BPCells matrix operations are available", {
    required <- c("colSums", "rowSums", "storage_order", "transpose_storage_order")

    expect_true(all(vapply(required, exists, logical(1), envir = asNamespace("BPCells"))))
})
