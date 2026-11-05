# compute_isotwas gives every method one set of folds (.isotwas_folds). For
# the methods that already drew exactly those folds themselves, passing them
# must change nothing.

same_model <- function(a, b) {
  expect_identical(names(a$transcripts), names(b$transcripts))
  for (tx in names(a$transcripts)) {
    expect_identical(a$transcripts[[tx]]$r2, b$transcripts[[tx]]$r2)
    expect_identical(a$transcripts[[tx]]$weights, b$transcripts[[tx]]$weights)
    expect_identical(a$transcripts[[tx]]$predicted,
                     b$transcripts[[tx]]$predicted)
  }
}

test_that(".isotwas_folds is createFolds under the seed and keeps the RNG state", {
  f <- withr::with_seed(5, {
    before <- .Random.seed
    f <- .isotwas_folds(60, 5, 71)
    expect_identical(.Random.seed, before)
    f
  })
  ref <- withr::with_seed(71, caret::createFolds(1:60, k = 5, returnTrain = TRUE))
  expect_identical(f$train, ref)
  expect_identical(f$nfolds, 5L)
  lab <- withr::with_seed(71, caret::createFolds(1:60, k = 5, list = FALSE))
  expect_identical(f$foldid, as.integer(lab))
})

test_that(".as_folds accepts training sets or fold labels, and checks them", {
  lab <- rep(1:4, length.out = 20)
  a <- .as_folds(lab, 20)
  b <- .as_folds(lapply(1:4, function(k) which(lab != k)), 20)
  expect_identical(a$foldid, b$foldid)
  expect_identical(a$train, b$train)
  expect_identical(.as_folds(a, 20), a)
  expect_error(.as_folds(lab, 21), "labels")
  expect_error(.as_folds(c(lab[-1], 6L), 20), "labels")
  expect_error(.as_folds(list(1:10, 5:20), 20), "exactly one fold")
  expect_error(.as_folds(a, 25), "samples")
})

test_that("methods that already drew the shared folds are unchanged by them", {
  X <- sim_geno(150, 30, 72)
  Y <- sim_expr(X, 3, 73)
  f <- .isotwas_folds(nrow(X), 5, 74)

  same_model(
    univariate_blup(X, Y, diag(3), nfolds = 5, verbose = FALSE, seed = 74),
    univariate_blup(X, Y, diag(3), nfolds = 5, verbose = FALSE, seed = 74,
                    folds = f))
  same_model(
    univariate_susie(X, Y, diag(3), nfolds = 5, verbose = FALSE, seed = 74),
    univariate_susie(X, Y, diag(3), nfolds = 5, verbose = FALSE, seed = 74,
                     folds = f))
  same_model(
    multivariate_spls(X, Y, nfolds = 5, verbose = FALSE, par = FALSE,
                      n.cores = 1, seed = 74),
    multivariate_spls(X, Y, nfolds = 5, verbose = FALSE, par = FALSE,
                      n.cores = 1, seed = 74, folds = f))
  same_model(
    multivariate_sgl(X, Y, nlambda = 5, nfolds = 5, seed = 74),
    multivariate_sgl(X, Y, nlambda = 5, nfolds = 5, seed = 74, folds = f))
  same_model(
    multivariate_graph_reg(X, Y, nlambda = 4, nlambda_graph = 2, nfolds = 5,
                           seed = 74),
    multivariate_graph_reg(X, Y, nlambda = 4, nlambda_graph = 2, nfolds = 5,
                           seed = 74, folds = f))
  # mtlasso's lambda folds are the shared ones; its R2 folds are drawn next
  same_model(
    multivariate_mtlasso(X, Y, nlambda = 5, nfolds = 5, seed = 74),
    multivariate_mtlasso(X, Y, nlambda = 5, nfolds = 5, seed = 74, folds = f))
  same_model(
    compute_mrce(X, Y, nlambda = 4, Omega = diag(3), nfolds = 5,
                 tol.in = 1e-3, maxit.in = 50, verbose = FALSE, seed = 74),
    compute_mrce(X, Y, nlambda = 4, Omega = diag(3), nfolds = 5,
                 tol.in = 1e-3, maxit.in = 50, verbose = FALSE, seed = 74,
                 folds = f))
})

test_that("user folds reach the methods, in either form", {
  X <- sim_geno(150, 30, 75)
  Y <- sim_expr(X, 3, 76)
  lab <- rep(1:4, length.out = nrow(X))
  a <- multivariate_sgl(X, Y, nlambda = 5, nfolds = 5, seed = 77, folds = lab)
  b <- multivariate_sgl(X, Y, nlambda = 5, nfolds = 5, seed = 77,
                        folds = lapply(1:4, function(k) which(lab != k)))
  same_model(a, b)
  d <- multivariate_sgl(X, Y, nlambda = 5, nfolds = 5, seed = 77)
  expect_false(identical(a$transcripts[[1]]$predicted,
                         d$transcripts[[1]]$predicted))
})
