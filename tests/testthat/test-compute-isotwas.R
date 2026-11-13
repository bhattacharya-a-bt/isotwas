# Tests for compute_isotwas()'s bookkeeping: the comparison matrix it prunes,
# the Omega it hands to each method, and what it reports when a method fails.
#
# The heavy fitting paths are covered by the per-method test files; these keep
# to the plumbing, which is where the silent failures were.

test_that("pruning all-NA columns keeps the comparison a matrix", {

    # `r2_mat[, keep]` drops to a vector when exactly one method survives, and
    # downstream code indexes it by column. Only one method surviving is the
    # normal case for a single-transcript gene, where the method set is forced
    # to `univariate`.
    r2_mat <- matrix(c(0.1, 0.2, NA, NA), nrow = 2,
                     dimnames = list(c("t1", "t2"), c("univariate", "spls")))

    keep <- apply(r2_mat, 2, function(x) !all(is.na(x)))
    pruned <- r2_mat[, keep, drop = FALSE]

    expect_true(is.matrix(pruned))
    expect_identical(dim(pruned), c(2L, 1L))
    expect_identical(colnames(pruned), "univariate")
})

test_that("a huge object's precision matrices live under $icov", {

    # The univariate methods were passed `omega_list[[omega_nlambda]]` while
    # every multivariate method was passed `omega_list$icov[[omega_nlambda]]`.
    # omega_list is a huge object, so the former indexes an arbitrary component
    # of that object rather than a precision matrix.
    skip_if_not_installed("huge")

    set.seed(41)
    Y <- matrix(stats::rnorm(50 * 4), 50, 4)
    om <- compute_omega(Y = Y, Y.rep = Y, R = 1, id = paste0("S", seq_len(50)),
                        method = "mean", nlambda = 5)

    expect_true(is.matrix(om$icov[[5]]))
    expect_identical(dim(om$icov[[5]]), c(4L, 4L))

    # What the univariate calls used to receive: not a 4x4 precision matrix.
    expect_false(is.matrix(om[[5]]) && identical(dim(om[[5]]), c(4L, 4L)))
})

test_that("a failed method is reported in the result, not only to the console", {

    # run_method() catches every error so one method's failure cannot abort the
    # run. That is the right call, but the failure was recorded only in a local
    # and printed under verbose, so at verbose = FALSE a method could vanish
    # from the comparison with nothing to show it had ever been attempted.
    res <- create_isotwas_result(
        best_models = create_isotwas_model(
            method = "best",
            transcripts = list(t1 = mk_tx(0.3)),
            n_samples = 100, n_snps = 10),
        failed_methods = data.frame(method = "Univariate ensemble (univariate)",
                                    message = "object 'x' not found",
                                    stringsAsFactors = FALSE))

    expect_identical(nrow(res$failed_methods), 1L)
    expect_identical(res$failed_methods$method, "Univariate ensemble (univariate)")
    expect_match(res$failed_methods$message, "not found")
})

test_that("a clean run reports no failures", {

    res <- create_isotwas_result(
        best_models = create_isotwas_model(
            method = "best",
            transcripts = list(t1 = mk_tx(0.3)),
            n_samples = 100, n_snps = 10),
        failed_methods = data.frame(method = character(0),
                                    message = character(0),
                                    stringsAsFactors = FALSE))

    expect_identical(nrow(res$failed_methods), 0L)
})

test_that("printing a result lists the methods that failed", {

    res <- create_isotwas_result(
        best_models = create_isotwas_model(
            method = "best",
            transcripts = list(t1 = mk_tx(0.3)),
            n_samples = 100, n_snps = 10),
        failed_methods = data.frame(method = "Joinet (joinet)",
                                    message = "boom",
                                    stringsAsFactors = FALSE))

    expect_output(print(res), "Methods that failed")
    expect_output(print(res), "Joinet \\(joinet\\): boom")
})

test_that("floor_r2 reaches a fitting method and only changes negatives", {

    # Y is pure noise against X, so every transcript is unpredictable and the
    # clamp is what is being measured. The clamp goes with the "oos" score
    # under r2_scale = 'legacy', and that score is negative on noise.
    skip_if_not_installed("tibble")

    d <- withr::with_seed(81, mk_xy(n = 60, p = 12, q = 3, h2 = 0))

    clamped <- multivariate_sgl(X = d$X, Y = d$Y, nfolds = 3, nlambda = 5,
                                verbose = FALSE, seed = 42, r2_method = "oos",
                                floor_r2 = TRUE)
    honest <- multivariate_sgl(X = d$X, Y = d$Y, nfolds = 3, nlambda = 5,
                               verbose = FALSE, seed = 42, r2_method = "oos",
                               floor_r2 = FALSE)

    expect_true(all(clamped$summary$r2 >= 0))
    expect_true(any(honest$summary$r2 < 0))

    # Identical wherever the clamp did not bite, and identical weights
    # throughout: the clamp is a reporting decision, not a fitting one.
    expect_equal(clamped$summary$r2, pmax(0, honest$summary$r2))
    expect_equal(lapply(clamped$transcripts, `[[`, "weights"),
                 lapply(honest$transcripts, `[[`, "weights"))
})
