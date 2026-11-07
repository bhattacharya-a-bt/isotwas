# multivariate_elasticnet on compressed rows must match the original rows.

withr::local_options(isotwas.glmnet_compress = TRUE)

test_that("multivariate_elasticnet is unchanged by compression", {
  X <- sim_geno(300, 50, 51)
  Y <- sim_expr(X, 3, 52)
  run <- function() multivariate_elasticnet(X, Y, diag(3), nfolds = 5,
                                            verbose = FALSE, seed = 53)
  expect_same_model(run(), uncompressed(run()))
})
