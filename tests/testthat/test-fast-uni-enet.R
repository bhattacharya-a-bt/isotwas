# univariate_elasticnet on compressed rows must match the original rows.

withr::local_options(isotwas.glmnet_compress = TRUE)

test_that("univariate_elasticnet is unchanged by compression", {
  X <- sim_geno(300, 50, 54)
  Y <- sim_expr(X, 3, 55)
  run <- function(par = FALSE) {
    univariate_elasticnet(X, Y, diag(3), nfolds = 5, verbose = FALSE,
                          par = par, n.cores = 2, seed = 56)
  }
  expect_same_model(run(), uncompressed(run()))
  skip_on_os("windows")
  expect_same_model(run(TRUE), uncompressed(run(TRUE)))
})
