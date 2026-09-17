# Tests for univariate_elasticnet(), in particular its parallel branch.

test_that("the parallel branch runs", {

    # The branch built its list of responses from an undefined lowercase `x`:
    #     lapply(seq_len(ncol(Y)), function(i) x[,i])
    # so every call raised "object 'x' not found". compute_isotwas() calls this
    # inside run_method()'s tryCatch, which turned the error into the silent
    # disappearance of the whole `univariate` method whenever par = TRUE.
    skip_if_not_installed("glmnet")
    skip_on_os("windows")   # mclapply is fork-based

    d <- withr::with_seed(31, mk_xy(n = 60, p = 20, q = 2))

    expect_no_error(
        univariate_elasticnet(X = d$X, Y = d$Y, Omega = NULL,
                              nfolds = 3, verbose = FALSE,
                              par = TRUE, n.cores = 2, seed = 1))
})

test_that("the parallel branch returns the same shape as the serial one", {

    skip_if_not_installed("glmnet")
    skip_on_os("windows")

    d <- withr::with_seed(32, mk_xy(n = 60, p = 20, q = 2))

    ser <- univariate_elasticnet(X = d$X, Y = d$Y, Omega = NULL,
                                 nfolds = 3, verbose = FALSE,
                                 par = FALSE, seed = 7)
    par <- univariate_elasticnet(X = d$X, Y = d$Y, Omega = NULL,
                                 nfolds = 3, verbose = FALSE,
                                 par = TRUE, n.cores = 2, seed = 7)

    expect_s3_class(par, "isotwas_model")
    expect_identical(names(par$transcripts), names(ser$transcripts))
    expect_identical(par$n_transcripts, ser$n_transcripts)
    expect_identical(names(par$transcripts[[1]]$weights),
                     names(ser$transcripts[[1]]$weights))
})

test_that("CHARACTERIZATION: parallel and serial disagree at the same seed", {

    # Not a consequence of the fix above: mclapply() forks, and the package
    # never sets RNGkind("L'Ecuyer-CMRG"), so each child continues from the
    # same RNG state rather than a distinct stream. cv.glmnet therefore draws
    # different folds under par = TRUE than under the sequential draws of
    # apply(), and repeated parallel runs are not reproducible either.
    #
    # Recorded here so that making parallel reproducible is a visible edit.
    skip_if_not_installed("glmnet")
    skip_on_os("windows")

    d <- withr::with_seed(32, mk_xy(n = 60, p = 20, q = 2))

    ser <- univariate_elasticnet(X = d$X, Y = d$Y, Omega = NULL,
                                 nfolds = 3, verbose = FALSE,
                                 par = FALSE, seed = 7)
    par <- univariate_elasticnet(X = d$X, Y = d$Y, Omega = NULL,
                                 nfolds = 3, verbose = FALSE,
                                 par = TRUE, n.cores = 2, seed = 7)

    r2_ser <- vapply(ser$transcripts, `[[`, numeric(1), "r2")
    r2_par <- vapply(par$transcripts, `[[`, numeric(1), "r2")

    expect_false(isTRUE(all.equal(r2_ser, r2_par)))
})

test_that("transcript names come through from the columns of Y", {

    skip_if_not_installed("glmnet")

    d <- withr::with_seed(33, mk_xy(n = 50, p = 12, q = 3))

    mod <- univariate_elasticnet(X = d$X, Y = d$Y, Omega = NULL,
                                 nfolds = 3, verbose = FALSE,
                                 par = FALSE, seed = 1)

    expect_identical(names(mod$transcripts), colnames(d$Y))
})
