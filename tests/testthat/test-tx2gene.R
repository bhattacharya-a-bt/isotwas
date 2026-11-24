# The documented rule, written out with lm(): each alpha's cross-validated
# path on the same folds; among lambdas with df > 0 whose held-out
# predictions correlate positively with the gene, the one with the best
# adjusted R2 of lm(gene ~ held-out prediction); and its coefficients.
ref_tx2gene <- function(P, gene, foldid, alphas = c(0, .25, .5, .75, 1)) {
  best <- list(r2 = NA_real_, w = rep(0, ncol(P)))
  for (a in alphas) {
    cv <- glmnet::cv.glmnet(P, gene, foldid = foldid, alpha = a, keep = TRUE)
    for (l in seq_along(cv$lambda)) {
      x <- cv$fit.preval[, l]
      if (cv$glmnet.fit$df[l] == 0 || stats::sd(x) == 0 ||
          stats::cor(x, gene) <= 0) next
      r2 <- summary(stats::lm(gene ~ x))$adj.r.sq
      if (is.na(best$r2) || r2 > best$r2) {
        best <- list(r2 = r2,
                     w = as.numeric(stats::coef(cv, s = cv$lambda[l]))[-1])
      }
    }
  }
  best
}

sim_tx2gene <- function(seed, n = 200, q = 4, signal = 1) {
  withr::with_seed(seed, {
    G <- matrix(stats::rnorm(n * q), n, q)
    list(P = G + matrix(stats::rnorm(n * q), n, q),
         gene = signal * rowSums(G) + stats::rnorm(n),
         noise = stats::rnorm(n))
  })
}

test_that(".fit_tx2gene returns the coefficients of the fit it scores", {
  foldid <- rep(1:5, length.out = 200)
  d <- sim_tx2gene(275)
  best <- ref_tx2gene(d$P, d$gene, foldid)
  got <- .fit_tx2gene(d$P, d$gene, foldid, paste0("t", 1:4))
  expect_equal(got$Weight_tx2gene, best$w, tolerance = 1e-10)
  expect_equal(got$R2, rep(best$r2, 4), tolerance = 1e-10)

  # a column with a tenth of the gene's standard deviation clears the floor
  P <- cbind(d$P, 0.1 * stats::sd(d$gene) * d$noise / stats::sd(d$noise))
  best <- ref_tx2gene(P, d$gene, foldid)
  got <- .fit_tx2gene(P, d$gene, foldid, paste0("t", 1:5))
  expect_equal(got$Weight_tx2gene, best$w, tolerance = 1e-10)
  expect_equal(got$R2[1], best$r2, tolerance = 1e-10)
})

test_that(".fit_tx2gene gives a nearly constant prediction weight 0", {
  foldid <- rep(1:5, length.out = 200)
  d <- sim_tx2gene(279)
  without <- .fit_tx2gene(cbind(d$P, 5), d$gene, foldid, paste0("t", 1:5))
  for (s in c(1e-3, 1e-12)) {
    P <- cbind(d$P, 5 + s * d$noise)
    got <- .fit_tx2gene(P, d$gene, foldid, paste0("t", 1:5))
    expect_identical(got$Weight_tx2gene, without$Weight_tx2gene)
    expect_identical(got$R2, without$R2)
  }

  # without the floor its weight is of order 1 / sd
  P <- cbind(d$P, 5 + 1e-5 * d$noise)
  unfloored <- .fit_tx2gene(P, d$gene, foldid, paste0("t", 1:5),
                            min_sd_ratio = 0)
  expect_gt(abs(unfloored$Weight_tx2gene[5]), 1e3)
})

test_that(".fit_tx2gene's floor is relative to the gene's scale", {
  foldid <- rep(1:5, length.out = 200)
  d <- sim_tx2gene(281)
  P <- cbind(d$P, 5 + 1e-4 * d$noise)
  got <- .fit_tx2gene(P, d$gene, foldid, paste0("t", 1:5))
  expect_identical(got$Weight_tx2gene[5], 0)
  for (k in c(1e-4, 1e4)) {
    scaled <- .fit_tx2gene(k * P, k * d$gene, foldid, paste0("t", 1:5))
    expect_equal(scaled$Weight_tx2gene, got$Weight_tx2gene, tolerance = 1e-6)
    expect_equal(scaled$R2, got$R2, tolerance = 1e-6)
  }
})

test_that(".fit_tx2gene does not score held-out predictions that anti-correlate", {
  # Null-end lambdas predict each fold's training mean, which anti-correlates
  # with the held-out values. Scored by unsigned R2, they won on weak genes
  # with weights of about 1e-37.
  foldid <- rep(1:5, length.out = 200)
  for (seed in 281:284) {
    d <- sim_tx2gene(seed, signal = 0.05)
    best <- ref_tx2gene(d$P, d$gene, foldid)
    got <- .fit_tx2gene(d$P, d$gene, foldid, paste0("t", 1:4))
    expect_equal(got$Weight_tx2gene, best$w, tolerance = 1e-10)
    expect_equal(got$R2[1], best$r2, tolerance = 1e-10)
    if (!is.na(got$R2[1])) {
      expect_gt(max(abs(got$Weight_tx2gene)), 1e-3)
    }
  }
})

test_that(".fit_tx2gene leaves out predictions constant within a training fold", {
  # As from multivariate_sgl selecting SNPs in one fold only. With every
  # column like that, glmnet stopped with error 7777.
  set.seed(276)
  n <- 200
  foldid <- rep(1:5, length.out = n)
  gene <- stats::rnorm(n)
  one_fold <- matrix(0.3, n, 2, dimnames = list(NULL, c("t1", "t2")))
  one_fold[foldid == 3, ] <- 0.3 + stats::rnorm(2 * sum(foldid == 3))

  got <- .fit_tx2gene(one_fold, gene, foldid, c("t1", "t2"))
  expect_identical(got$Weight_tx2gene, c(0, 0))
  expect_true(all(is.na(got$R2)))

  # beside an ordinary prediction, the fit is the one without them
  t3 <- gene + stats::rnorm(n)
  got <- .fit_tx2gene(cbind(one_fold, t3), gene, foldid, c("t1", "t2", "t3"))
  alone <- .fit_tx2gene(cbind(0.3, 0.3, t3), gene, foldid, c("t1", "t2", "t3"))
  expect_gt(got$Weight_tx2gene[3], 0)
  expect_identical(got, alone)
})

test_that(".fit_tx2gene skips the fit when nothing can be fitted", {
  set.seed(271)
  gene <- stats::rnorm(60)
  foldid <- rep(1:5, length.out = 60)
  some <- cbind(t1 = gene + stats::rnorm(60), t2 = stats::rnorm(60))
  cases <- list(
    flat = list(cbind(rep(0.2, 60), rep(-1, 60)), gene),
    below_floor = list(cbind(0.2 + 1e-6 * stats::rnorm(60),
                             -1 + 1e-9 * stats::rnorm(60)), gene),
    constant_gene = list(some, rep(3, 60)))
  for (case in cases) {
    got <- .fit_tx2gene(case[[1]], case[[2]], foldid, c("t1", "t2"))
    expect_identical(got$Feature, c("t1", "t2"))
    expect_identical(got$Weight_tx2gene, c(0, 0))
    expect_true(all(is.na(got$R2)))
  }
})

test_that(".fit_tx2gene leaves out non-finite predictions, with a warning", {
  set.seed(271)
  gene <- stats::rnorm(60)
  foldid <- rep(1:5, length.out = 60)
  some <- cbind(t1 = gene + stats::rnorm(60),
                t2 = c(NA, stats::rnorm(59)), t3 = c(Inf, stats::rnorm(59)))
  expect_warning(got <- .fit_tx2gene(some, gene, foldid, paste0("t", 1:3)),
                 "t2, t3")
  expect_gt(got$Weight_tx2gene[1], 0)
  expect_identical(got$Weight_tx2gene[2:3], c(0, 0))
})

test_that(".fit_tx2gene fits a single transcript", {
  # glmnet needs two columns; this stopped with "x should be a matrix with 2
  # or more columns"
  foldid <- rep(1:5, length.out = 200)
  d <- sim_tx2gene(282, q = 1)
  best <- ref_tx2gene(cbind(d$P, 0), d$gene, foldid)
  got <- .fit_tx2gene(d$P, d$gene, foldid, "t1")
  expect_gt(got$Weight_tx2gene, 0)
  expect_equal(got$Weight_tx2gene, best$w[1], tolerance = 1e-10)
  expect_equal(got$R2, best$r2, tolerance = 1e-10)
})

test_that("compute_isotwas checks gene_exp before fitting", {
  X <- sim_geno(40, 10, 283)
  Y <- withr::with_seed(284, matrix(stats::rnorm(80), 40, 2))
  for (bad in list(c(NA, rowSums(Y)[-1]), rowSums(Y)[-1])) {
    expect_error(
      compute_isotwas(X, Y, gene_exp = bad, Y.rep = Y, R = 1,
                      id = seq_len(40), method = "univariate",
                      verbose = FALSE),
      "gene_exp must have one finite value per sample")
  }
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
