# Helpers for test-fast-*.R files that compare a glmnet-based method with
# compression on and off.

# With options(isotwas.glmnet_compress = FALSE) the glmnet-based methods
# pass the original rows to glmnet exactly as cv.glmnet did (see
# test-fast-glmnet-compress.R), so comparing the two routes checks that
# compression leaves each method's output unchanged.

expect_same_model <- function(new, ref, tol = 1e-10) {
  expect_identical(names(new$transcripts), names(ref$transcripts))
  for (tx in names(ref$transcripts)) {
    a <- new$transcripts[[tx]]
    b <- ref$transcripts[[tx]]
    expect_equal(a$r2, b$r2, tolerance = tol)
    expect_equal(a$pvalue, b$pvalue, tolerance = 1e-8)
    expect_equal(a$predicted, b$predicted, tolerance = tol)
    wa <- setNames(a$weights$Weight, a$weights$SNP)
    wb <- setNames(b$weights$Weight, b$weights$SNP)
    snps <- union(names(wa), names(wb))
    wa <- wa[snps]; wa[is.na(wa)] <- 0
    wb <- wb[snps]; wb[is.na(wb)] <- 0
    expect_lt(max(abs(wa - wb)), tol)
  }
}

uncompressed <- function(expr) {
  withr::with_options(list(isotwas.glmnet_compress = FALSE), expr)
}
