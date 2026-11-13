# The transcript-to-gene model regresses gene expression on the selected
# models' cross-validated predictions. glmnet cannot fit when none of them
# varies ("All used predictors have zero variance"), so the fit is skipped.

test_that(".fit_tx2gene skips the fit when no prediction varies", {
  set.seed(271)
  gene <- stats::rnorm(60)
  foldid <- rep(1:5, length.out = 60)
  flat <- cbind(t1 = rep(0.2, 60), t2 = rep(-1, 60))
  got <- .fit_tx2gene(flat, gene, foldid, c("t1", "t2"))
  expect_identical(got$Feature, c("t1", "t2"))
  expect_identical(got$Weight_tx2gene, c(0, 0))
  expect_true(all(is.na(got$R2)))

  # one varying column is enough to fit
  some <- cbind(t1 = gene + stats::rnorm(60), t2 = rep(-1, 60))
  got <- .fit_tx2gene(some, gene, foldid, c("t1", "t2"))
  expect_true(is.finite(got$R2[1]))
  expect_gt(got$Weight_tx2gene[1], 0)
})

test_that("compute_isotwas returns when no selected model varies", {
  # On pure noise sgl stays at the null end of its path, so every
  # cross-validated prediction is the expression mean. This stopped with
  # glmnet's error 7777 before the fit was skipped.
  skip_if_not_installed("huge")
  X <- sim_geno(100, 20, 272)
  Y <- withr::with_seed(273, matrix(stats::rnorm(300), 100, 3,
                                    dimnames = list(NULL, paste0("tx", 1:3))))
  res <- compute_isotwas(X, Y, gene_exp = rowSums(Y), Y.rep = Y, R = 1,
                         id = seq_len(100), omega_est = "mean",
                         omega_nlambda = 2, method = "sgl", run_all = FALSE,
                         sgl_nlambda = 4, verbose = FALSE, seed = 274)
  expect_true(all(is.na(res$best_models$best_methods$best_method)))
  expect_identical(res$tx2gene$Weight_tx2gene, c(0, 0, 0))
  expect_true(all(is.na(res$tx2gene$R2)))
})
