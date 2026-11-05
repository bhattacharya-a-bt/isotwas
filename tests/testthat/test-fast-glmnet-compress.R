# .glmnet_compressed() and .cv_glmnet() must reproduce glmnet() and
# cv.glmnet() on the original rows.

# Compression on, so these small problems take the compressed route
withr::local_options(isotwas.glmnet_compress = TRUE)

coef_diff <- function(a, b) {
  if (is.list(a)) return(max(mapply(coef_diff, a, b)))
  max(abs(as.matrix(a) - as.matrix(b)))
}

test_that("compressed rows give glmnet's path (gaussian and mgaussian)", {
  X <- sim_geno(300, 60, 31)
  X[, 2:3] <- X[, 1]   # duplicated SNPs: rank-deficient X'X
  X[, 4] <- 1          # constant SNP
  Y <- sim_expr(X, 3, 32)
  d <- .glmnet_design(.gram_subset(.gram_full(X, Y)))
  for (a in c(0, .5)) {
    ref <- glmnet::glmnet(X, Y[, 2], alpha = a)
    new <- .glmnet_compressed(d, 2, "gaussian", alpha = a)
    expect_equal(new$lambda, ref$lambda, tolerance = 1e-12)
    expect_lt(coef_diff(coef(new), coef(ref)), 1e-10)
    expect_equal(new$nulldev, ref$nulldev, tolerance = 1e-12)
    expect_equal(new$dev.ratio, ref$dev.ratio, tolerance = 1e-10)

    ref <- glmnet::glmnet(X, Y, family = "mgaussian", alpha = a,
                          standardize = FALSE)
    new <- .glmnet_compressed(d, family = "mgaussian", alpha = a,
                              standardize = FALSE)
    expect_equal(new$lambda, ref$lambda, tolerance = 1e-12)
    expect_lt(coef_diff(coef(new), coef(ref)), 1e-10)
    expect_equal(dimnames(predict(new, X[1:5, ], s = .1)),
                 dimnames(predict(ref, X[1:5, ], s = .1)))
  }
})

test_that("a single user lambda fit matches glmnet's", {
  X <- sim_geno(300, 60, 33)
  Y <- sim_expr(X, 2, 34)
  d <- .glmnet_design(.gram_subset(.gram_full(X, Y)))
  ref <- glmnet::glmnet(X, Y[, 1], alpha = .5, lambda = .05)
  new <- .glmnet_compressed(d, 1, "gaussian", alpha = .5, lambda = .05)
  expect_lt(coef_diff(coef(new), coef(ref)), 1e-10)
})

test_that("subset cross-products equal those of the subset's rows", {
  X <- sim_geno(200, 30, 35)
  Y <- sim_expr(X, 2, 36)
  drop <- seq(1, 200, by = 4)
  s <- .gram_subset(.gram_full(X, Y), drop)
  Xs <- scale(X[-drop, ], scale = FALSE)
  Ys <- scale(Y[-drop, ], scale = FALSE)
  expect_equal(s$n, 150)
  expect_equal(unname(s$xx), unname(crossprod(Xs)), tolerance = 1e-12)
  expect_equal(unname(s$xy), unname(crossprod(Xs, Ys)), tolerance = 1e-12)
  expect_equal(unname(s$xm), unname(colMeans(X[-drop, ])), tolerance = 1e-12)
})

test_that(".cv_glmnet reproduces cv.glmnet", {
  X <- sim_geno(400, 60, 37)
  Y <- sim_expr(X, 3, 38)
  foldid <- withr::with_seed(39, sample(rep(1:5, length = nrow(X))))
  for (a in c(0, .5, 1)) {
    ref <- glmnet::cv.glmnet(X, Y[, 1], foldid = foldid, alpha = a, keep = TRUE)
    new <- .cv_glmnet(X, Y[, 1], foldid, "gaussian", alpha = a)
    expect_equal(new$lambda, ref$lambda, tolerance = 1e-12)
    expect_equal(new$cvm, ref$cvm, tolerance = 1e-10)
    expect_equal(new$cvsd, ref$cvsd, tolerance = 1e-8)
    expect_equal(new$lambda.min, ref$lambda.min)
    expect_equal(new$fit.preval, ref$fit.preval, tolerance = 1e-10)
    expect_lt(coef_diff(coef(new, s = "lambda.min"),
                        coef(ref, s = "lambda.min")), 1e-10)

    ref <- glmnet::cv.glmnet(X, Y, foldid = foldid, alpha = a, keep = TRUE,
                             family = "mgaussian", standardize = FALSE)
    new <- .cv_glmnet(X, Y, foldid, "mgaussian", alpha = a,
                      standardize = FALSE)
    expect_equal(new$cvm, ref$cvm, tolerance = 1e-10)
    expect_equal(new$lambda.min, ref$lambda.min)
    expect_equal(new$fit.preval, ref$fit.preval, tolerance = 1e-10)
    expect_lt(coef_diff(coef(new, s = "lambda.min"),
                        coef(ref, s = "lambda.min")), 1e-10)
  }
})

test_that(".cv_glmnet passes the original rows to glmnet when compression is off", {
  X <- sim_geno(120, 150, 40)   # p > n
  Y <- sim_expr(X, 2, 41)
  foldid <- withr::with_seed(42, sample(rep(1:5, length = nrow(X))))
  ref <- glmnet::cv.glmnet(X, Y[, 1], foldid = foldid, alpha = .5, keep = TRUE)
  new <- .cv_glmnet(X, Y[, 1], foldid, "gaussian", alpha = .5)
  expect_identical(new$cvm, ref$cvm)
  expect_identical(new$fit.preval, ref$fit.preval)
  expect_identical(coef(new, s = "lambda.min"), coef(ref, s = "lambda.min"))

  X <- sim_geno(300, 40, 43)
  Y <- sim_expr(X, 2, 44)
  foldid <- withr::with_seed(45, sample(rep(1:5, length = nrow(X))))
  ref <- glmnet::cv.glmnet(X, Y, foldid = foldid, alpha = .5, keep = TRUE,
                           family = "mgaussian", standardize = FALSE)
  new <- withr::with_options(list(isotwas.glmnet_compress = FALSE),
    .cv_glmnet(X, Y, foldid, "mgaussian", alpha = .5, standardize = FALSE))
  expect_identical(new$cvm, ref$cvm)
  expect_identical(new$fit.preval, ref$fit.preval)
})

test_that("designs can compress some fits and leave others on the original rows", {
  X <- sim_geno(300, 40, 49)
  Y <- sim_expr(X, 2, 50)
  foldid <- withr::with_seed(51, sample(rep(1:5, length = nrow(X))))
  d <- .glmnet_cv_designs(X, Y, foldid)
  d$folds[c(2, 4)] <- list(NULL)
  ref <- glmnet::cv.glmnet(X, Y[, 2], foldid = foldid, alpha = .5, keep = TRUE)
  new <- .cv_glmnet(X, Y[, 2], foldid, "gaussian", alpha = .5, designs = d,
                    j = 2)
  expect_equal(new$cvm, ref$cvm, tolerance = 1e-10)
  expect_equal(new$fit.preval, ref$fit.preval, tolerance = 1e-10)
})

test_that("compression is chosen by the option and the rows it saves", {
  expect_true(.glmnet_compress_on(1000, 100, 1L))
  expect_true(.glmnet_compress_on(1000, 748, 1L))    # 750 rows = 0.75 n
  expect_false(.glmnet_compress_on(1000, 749, 1L))
  expect_false(.glmnet_compress_on(1000, 740, 10L))
  withr::with_options(list(isotwas.glmnet_compress_ratio = 0.5),
                      expect_false(.glmnet_compress_on(1000, 600, 1L)))
  withr::with_options(list(isotwas.glmnet_compress = FALSE),
                      expect_false(.glmnet_compress_on(1000, 100, 1L)))
  withr::with_options(list(isotwas.glmnet_compress = NULL),
                      expect_false(.glmnet_compress_on(1000, 100, 1L)))
})

test_that("shared designs serve each response of a gaussian fit", {
  X <- sim_geno(300, 40, 46)
  Y <- sim_expr(X, 3, 47)
  foldid <- withr::with_seed(48, sample(rep(1:5, length = nrow(X))))
  designs <- .glmnet_cv_designs(X, Y, foldid)
  for (j in 1:3) {
    ref <- glmnet::cv.glmnet(X, Y[, j], foldid = foldid, alpha = .5, keep = TRUE)
    new <- .cv_glmnet(X, Y[, j], foldid, "gaussian", alpha = .5,
                      designs = designs, j = j)
    expect_equal(new$cvm, ref$cvm, tolerance = 1e-10)
    expect_equal(new$fit.preval, ref$fit.preval, tolerance = 1e-10)
  }
})
