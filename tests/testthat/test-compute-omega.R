# Tests for compute_omega().
#
# As in test-utils.R, tests marked CHARACTERIZATION record current behaviour
# -- including behaviour that looks unintended -- so that changing it is a
# visible edit here rather than a silent drift.

# ---- what the 'replicates' branch actually computes ---------------------
#
# These need no dependencies: they check the linear algebra the branch
# performs, independently of the glasso that consumes it.

test_that("the ~id projection returns fitted values, i.e. per-sample means", {

    # compute_omega's 'replicates' branch assigns
    #     Y.resid = mm %*% solve(t(mm) %*% mm) %*% t(mm) %*% Y.rep
    # which is H %*% Y.rep. For the one-way design ~id, H %*% Y.rep replaces
    # every row with its sample's across-replicate mean -- so despite the name
    # Y.resid, the within-sample variation has been removed, not extracted.
    set.seed(21)
    d <- mk_reps(N = 6, G = 3, R = 4)

    mm <- stats::model.matrix(~ d$id)
    fitted <- mm %*% solve(t(mm) %*% mm) %*% t(mm) %*% d$Y.rep

    grp <- match(d$id, unique(d$id))
    group_means <- (rowsum(d$Y.rep, grp) / as.vector(table(grp)))[grp, , drop = FALSE]

    expect_equal(unname(fitted), unname(group_means))
})

test_that("at R = 1 the ~id projection is the identity", {

    # With one row per sample the design is square and full rank, so H = I and
    # the 'replicates' branch receives exactly Y -- making it identical to
    # method = 'mean'.
    set.seed(22)
    Y <- matrix(stats::rnorm(8 * 3), 8, 3)
    id <- paste0("S", seq_len(8))

    mm <- stats::model.matrix(~id)
    fitted <- mm %*% solve(t(mm) %*% mm) %*% t(mm) %*% Y

    expect_equal(unname(fitted), unname(Y))
})

test_that("the true residual is the group-mean-subtracted matrix", {

    # Stated here as the reference the two independent implementations agree
    # on, for a balanced and an unbalanced design.
    set.seed(23)
    d <- mk_reps(N = 6, G = 3, R = 4)
    expect_equal(unname(resid_hat(d$Y.rep, d$id)),
                 unname(resid_qr(d$Y.rep, d$id)))

    u <- mk_reps(unbalanced = TRUE)
    expect_equal(unname(resid_hat(u$Y.rep, u$id)),
                 unname(resid_qr(u$Y.rep, u$id)))
})

# ---- the function itself ------------------------------------------------

test_that("method = 'mean' returns a huge object with nlambda precision matrices", {

    skip_if_not_installed("huge")

    set.seed(24)
    Y <- matrix(stats::rnorm(50 * 4), 50, 4)

    om <- compute_omega(Y = Y, Y.rep = Y, R = 1, id = paste0("S", seq_len(50)),
                        method = "mean", nlambda = 5)

    expect_length(om$icov, 5)
    expect_identical(dim(om$icov[[1]]), c(4L, 4L))
})

test_that("at R = 1, 'replicates' is identical to 'mean'", {

    # The consequence of H = I above, at the level of the returned object. The
    # isoTWAS manuscript ran method = 'replicates' with Y.rep = Y and R = 1, so
    # its precision matrices were those of method = 'mean'.
    skip_if_not_installed("huge")

    set.seed(25)
    Y <- matrix(stats::rnorm(50 * 4), 50, 4)
    id <- paste0("S", seq_len(50))

    reps <- compute_omega(Y = Y, Y.rep = Y, R = 1, id = id,
                          method = "replicates", nlambda = 5)
    mean <- compute_omega(Y = Y, Y.rep = Y, R = 1, id = id,
                          method = "mean", nlambda = 5)

    expect_equal(reps$icov, mean$icov)
})

test_that("the documented default for `method` is usable", {

    # The signature is method = c('replicates','mean'). Without match.arg(),
    # the branch test `if (method == 'mean')` received a length-2 condition,
    # which has been an error since R 4.2 -- so the exported function could not
    # be called at its own default.
    skip_if_not_installed("huge")

    set.seed(26)
    Y <- matrix(stats::rnorm(20 * 3), 20, 3)

    expect_no_error(
        compute_omega(Y = Y, Y.rep = Y, R = 1, id = paste0("S", seq_len(20)),
                      nlambda = 3))
})

test_that("an unrecognised `method` is rejected by name", {

    skip_if_not_installed("huge")

    set.seed(26)
    Y <- matrix(stats::rnorm(20 * 3), 20, 3)

    expect_error(
        compute_omega(Y = Y, Y.rep = Y, R = 1, id = paste0("S", seq_len(20)),
                      method = "nonsense", nlambda = 3),
        "should be one of")
})

test_that("an error inside the function does not leak the sink", {

    # Output is suppressed with sink(tempfile()). Released by a bare sink() at
    # the end of the body, any error in between skipped the release and left
    # the caller's console redirected for the rest of the session.
    skip_if_not_installed("huge")

    before <- sink.number()
    on.exit(while (sink.number() > before) sink(), add = TRUE)

    # Y.rep with the wrong number of rows for `id` makes model.matrix and the
    # matrix product disagree, so the failure happens inside the sunk region.
    set.seed(27)
    Y <- matrix(stats::rnorm(20 * 3), 20, 3)
    try(compute_omega(Y = Y, Y.rep = Y, R = 1, id = paste0("S", seq_len(19)),
                      method = "replicates", nlambda = 3),
        silent = TRUE)

    expect_identical(sink.number(), before)
})

# ---- the 'residuals' mode ----------------------------------------------

test_that(".group_project matches the explicit hat matrix, balanced and not", {

    # The whole justification for not forming H: group-mean subtraction is
    # the exact projection for the one-way design ~id, for equal and unequal
    # group sizes alike, at O(nrow * ncol) instead of O(nrow^2).
    set.seed(71)
    d <- mk_reps(N = 12, G = 4, R = 5)
    mm <- stats::model.matrix(~ d$id)
    H <- mm %*% solve(t(mm) %*% mm) %*% t(mm)

    expect_equal(unname(.group_project(d$Y.rep, d$id, "fitted")),
                 unname(H %*% d$Y.rep))
    expect_equal(unname(.group_project(d$Y.rep, d$id, "residual")),
                 unname(d$Y.rep - H %*% d$Y.rep))
    expect_equal(unname(.group_project(d$Y.rep, d$id, "residual")),
                 unname(resid_qr(d$Y.rep, d$id)))

    u <- mk_reps(unbalanced = TRUE)
    expect_equal(unname(.group_project(u$Y.rep, u$id, "residual")),
                 unname(resid_hat(u$Y.rep, u$id)))
})

test_that(".group_project rejects a mismatched id", {

    set.seed(72)
    d <- mk_reps(N = 6, G = 3, R = 3)

    expect_error(.group_project(d$Y.rep, d$id[-1], "residual"),
                 "must equal length")
})

test_that("anything constant within a sample cancels in the residual", {

    # Batch correction and covariate residualisation are applied per sample,
    # so they leave the within-sample deviations untouched.
    set.seed(73)
    d <- mk_reps(N = 10, G = 4, R = 6)
    offset <- matrix(stats::rnorm(10 * 4, sd = 5), 10, 4)
    grp <- match(d$id, unique(d$id))
    shifted <- d$Y.rep + offset[grp, , drop = FALSE]

    expect_equal(.group_project(shifted, d$id, "residual"),
                 .group_project(d$Y.rep, d$id, "residual"))
})

test_that("'residuals' differs from 'replicates' when replicates vary", {

    skip_if_not_installed("huge")

    set.seed(74)
    d <- mk_reps(N = 30, G = 4, R = 4, sd = 0.5)

    reps <- compute_omega(Y = d$Y, Y.rep = d$Y.rep, R = 4, id = d$id,
                          method = "replicates", nlambda = 4)
    resid <- compute_omega(Y = d$Y, Y.rep = d$Y.rep, R = 4, id = d$id,
                           method = "residuals", nlambda = 4)

    expect_false(isTRUE(all.equal(reps$icov, resid$icov)))
})

test_that("'residuals' refuses a design with one replicate per sample", {

    # The residual is identically zero there, so the within-sample covariance
    # is unidentifiable. Silently degrading to 'mean' would hide that.
    skip_if_not_installed("huge")

    set.seed(75)
    Y <- matrix(stats::rnorm(20 * 3), 20, 3)

    expect_error(
        compute_omega(Y = Y, Y.rep = Y, R = 1, id = paste0("S", seq_len(20)),
                      method = "residuals", nlambda = 3),
        ">= 2 replicates")
})

test_that("at R = 1 the residual is identically zero", {

    set.seed(76)
    Y <- matrix(stats::rnorm(15 * 3), 15, 3)

    expect_equal(unname(.group_project(Y, paste0("S", seq_len(15)), "residual")),
                 matrix(0, 15, 3))
})
