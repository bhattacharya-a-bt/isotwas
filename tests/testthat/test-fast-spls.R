# .cv_spls() and .spls_fit() must reproduce spls::cv.spls() and spls::spls().

test_that(".cv_spls and .spls_fit reproduce spls::cv.spls and spls::spls", {
  X <- sim_geno(90, 30, 18)
  Y <- sim_expr(X, 3, 19)
  ref <- withr::with_seed(20, {
    utils::capture.output(cv <- spls::cv.spls(X, Y, fold = 3, K = 1:3,
                                               eta = seq(.1, .9, .2),
                                               scale.x = FALSE, scale.y = FALSE,
                                               plot.it = FALSE))
    cv
  })
  new <- withr::with_seed(20, .cv_spls(X, Y, fold = 3, K = 1:3, eta = seq(.1, .9, .2)))
  expect_equal(new$mspemat, ref$mspemat, tolerance = 1e-10, ignore_attr = TRUE)
  expect_equal(new$K.opt, ref$K.opt)
  expect_equal(new$eta.opt, ref$eta.opt)

  fit_ref <- spls::spls(X, Y, K = 2, eta = .5, scale.x = FALSE, scale.y = FALSE)
  fit_new <- .spls_fit(.spls_prep(X, Y), K = 2, eta = .5)
  expect_equal(fit_new$betahat, fit_ref$betahat, tolerance = 1e-10,
               ignore_attr = TRUE)
  expect_equal(fit_new$A, fit_ref$A)
})
