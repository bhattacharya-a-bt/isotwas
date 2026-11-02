# With an intercept, adding a constant to expression changes only the
# intercept. Without one, the weights absorb the mean.

weights_of <- function(mod) {
    lapply(mod$transcripts, function(t) setNames(t$weights$Weight,
                                                 t$weights$SNP))
}

test_that("univariate_elasticnet weights do not depend on the expression mean", {

    skip_if_not_installed("glmnet")

    d <- withr::with_seed(61, mk_xy(n = 80, p = 20, q = 2, h2 = 0.5))

    fit <- function(Y) univariate_elasticnet(X = d$X, Y = Y, Omega = NULL,
                                             nfolds = 3, verbose = FALSE,
                                             par = FALSE, seed = 1)
    centred <- fit(d$Y)
    shifted <- fit(d$Y + 7)

    expect_equal(weights_of(shifted), weights_of(centred), tolerance = 1e-6)
})

test_that("univariate_susie weights do not depend on the expression mean", {

    skip_if_not_installed("susieR")

    d <- withr::with_seed(62, mk_xy(n = 80, p = 20, q = 2, h2 = 0.5))

    fit <- function(Y) univariate_susie(X = d$X, Y = Y, Omega = NULL,
                                        nfolds = 3, verbose = FALSE,
                                        par = FALSE, seed = 1)
    centred <- fit(d$Y)
    shifted <- fit(d$Y + 7)

    expect_equal(weights_of(shifted), weights_of(centred), tolerance = 1e-6)
})
