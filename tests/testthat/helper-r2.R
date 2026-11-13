# Checks that a fitted model reports the cross-validated R2 of the stated
# definition (see .score_r2), for the tests of r2_scale.

#' Expect every transcript's R2 to be its prediction scored by `method`
#'
#' A model emptied for having no signal reports 0. At least one transcript
#' must have a real model, and on one of them the two definitions must
#' differ, so that the check could fail.
#'
#' @param model an isotwas_model
#' @param Y matrix, the expression it was fitted to, columns named as its
#'   transcripts
#' @param method character, "refit" or "oos"
#' @param floor logical, the score is clamped at 0
expect_scored_by <- function(model, Y, method, floor = FALSE) {
  other <- setdiff(c("refit", "oos"), method)
  n_scored <- 0L
  gap <- 0
  for (tx in names(model$transcripts)) {
    tm <- model$transcripts[[tx]]
    if (is.null(tm$weights) || nrow(tm$weights) == 0) {
      expect_identical(tm$r2, 0)
      next
    }
    score <- function(m) {
      r2 <- .score_r2(Y[, tx], tm$predicted, m, pvalue = FALSE)$R2
      if (floor) max(0, r2) else r2
    }
    expect_equal(tm$r2, score(method), tolerance = 1e-12)
    gap <- max(gap, abs(score(method) - score(other)))
    n_scored <- n_scored + 1L
  }
  expect_gt(n_scored, 0L)
  expect_gt(gap, 1e-6)
  invisible(n_scored)
}
