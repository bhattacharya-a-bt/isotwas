# Parallel fitting must give the same answer every time, and the same answer
# regardless of how many cores it was given. Tests use at most 2 cores, the
# limit R CMD check --as-cran enforces.

test_that(".mclapply_seeded is reproducible and independent of core count", {

    skip_on_os("windows")   # mclapply is fork-based

    draw <- function(i) stats::runif(1)

    one <- unlist(.mclapply_seeded(1:8, draw, mc.cores = 1, seed = 42))
    two <- unlist(.mclapply_seeded(1:8, draw, mc.cores = 2, seed = 42))
    again <- unlist(.mclapply_seeded(1:8, draw, mc.cores = 2, seed = 42))

    expect_equal(one, two)
    expect_equal(two, again)
})

test_that(".mclapply_seeded gives different answers for different seeds", {

    skip_on_os("windows")

    draw <- function(i) stats::runif(1)

    a <- unlist(.mclapply_seeded(1:5, draw, mc.cores = 1, seed = 1))
    b <- unlist(.mclapply_seeded(1:5, draw, mc.cores = 1, seed = 2))

    expect_false(isTRUE(all.equal(a, b)))
})

test_that(".mclapply_seeded leaves the caller's RNG stream alone", {

    skip_on_os("windows")

    set.seed(99)
    before <- .Random.seed
    invisible(.mclapply_seeded(1:4, function(i) stats::runif(1),
                               mc.cores = 2, seed = 7))

    expect_identical(.Random.seed, before)
})

test_that(".mclapply_seeded passes extra arguments through", {

    skip_on_os("windows")

    got <- .mclapply_seeded(1:3, function(x, k) x * k, k = 10,
                            mc.cores = 1, seed = 1)

    expect_equal(unlist(got), c(10, 20, 30))
})

test_that("univariate_elasticnet is reproducible in parallel", {

    skip_if_not_installed("glmnet")
    skip_on_os("windows")

    d <- withr::with_seed(91, mk_xy(n = 60, p = 15, q = 2))

    fit <- function(cores) {
        univariate_elasticnet(X = d$X, Y = d$Y, Omega = NULL, nfolds = 3,
                              verbose = FALSE, par = TRUE, n.cores = cores,
                              seed = 5)
    }

    # With par = TRUE, one core still goes through the seeded path.
    a <- fit(2)
    b <- fit(2)
    c1 <- fit(1)

    r2 <- function(m) vapply(m$transcripts, `[[`, numeric(1), "r2")
    expect_equal(r2(a), r2(b))
    expect_equal(r2(a), r2(c1))
    expect_equal(lapply(a$transcripts, `[[`, "weights"),
                 lapply(c1$transcripts, `[[`, "weights"))
})

test_that("multivariate_sgl is reproducible in parallel", {

    skip_if_not_installed("tibble")
    skip_on_os("windows")

    d <- withr::with_seed(92, mk_xy(n = 50, p = 12, q = 3))

    fit <- function(cores) {
        multivariate_sgl(X = d$X, Y = d$Y, nfolds = 3, nlambda = 5,
                         verbose = FALSE, seed = 3, par = TRUE, n.cores = cores)
    }

    # One core takes the serial path, so only repeatability is checked here;
    # core-count independence is covered by the .mclapply_seeded test above.
    expect_equal(fit(2)$summary$r2, fit(2)$summary$r2)
})

test_that("the seed a run used is recorded in the result", {

    # Without this, a run left to the default seed cannot be reproduced: the
    # seed was drawn internally and reported only under verbose.
    res <- create_isotwas_result(
        best_models = create_isotwas_model(
            method = "best",
            transcripts = list(t1 = mk_tx(0.3)),
            n_samples = 100, n_snps = 10),
        seed = 4242)

    expect_identical(res$seed, 4242)
})
