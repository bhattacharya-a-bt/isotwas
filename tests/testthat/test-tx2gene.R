# The transcript-to-gene model regresses gene expression on the selected
# models' cross-validated predictions with an elastic net whose alpha and
# lambda are chosen by cross-validation. glmnet cannot fit when none of the
# predictions varies ("All used predictors have zero variance"), so the fit
# is then skipped.

test_that(".fit_tx2gene returns the coefficients of the fit it scores", {
  # Reference: each alpha's cross-validated path on the same folds, the
  # alpha and lambda with the best adjusted R2 of lm(gene ~ held-out
  # prediction), and that fit's coefficients. The code this replaced chose
  # among relaxed-lasso gammas, then refitted with alpha set to the chosen
  # gamma, so its weights came from a different model than its R2.
  set.seed(275)
  n <- 200; q <- 4
  G <- matrix(stats::rnorm(n * q), n, q)
  P <- G + matrix(stats::rnorm(n * q), n, q)
  gene <- rowSums(G) + stats::rnorm(n)
  foldid <- rep(1:5, length.out = n)
  alphas <- c(0, .25, .5, .75, 1)

  best <- list(r2 = -Inf)
  for (a in alphas) {
    cv <- glmnet::cv.glmnet(P, gene, foldid = foldid, alpha = a, keep = TRUE)
    for (l in seq_along(cv$lambda)) {
      x <- cv$fit.preval[, l]
      if (stats::sd(x) == 0) next
      r2 <- summary(stats::lm(gene ~ x))$adj.r.sq
      if (r2 > best$r2) {
        best <- list(r2 = r2,
                     w = as.numeric(stats::coef(cv, s = cv$lambda[l]))[-1])
      }
    }
  }

  got <- .fit_tx2gene(P, gene, foldid, paste0("t", 1:q))
  expect_equal(got$Weight_tx2gene, best$w, tolerance = 1e-10)
  expect_equal(got$R2, rep(best$r2, q), tolerance = 1e-10)
})

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
