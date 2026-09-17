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

test_that("CHARACTERIZATION: at R = 1, 'replicates' is identical to 'mean'", {

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

test_that("CHARACTERIZATION: the documented default for `method` errors", {

    # The signature is method = c('replicates','mean') but there is no
    # match.arg(), so the branch test `if (method == 'mean')` receives a
    # length-2 condition. That has been an error since R 4.2, which makes the
    # exported function unusable at its own default.
    skip_if_not_installed("huge")

    set.seed(26)
    Y <- matrix(stats::rnorm(20 * 3), 20, 3)

    # The call leaks a sink on the way out (see below), so release it here.
    on.exit(while (sink.number() > 0) sink(), add = TRUE)

    expect_error(
        compute_omega(Y = Y, Y.rep = Y, R = 1, id = paste0("S", seq_len(20))),
        "condition has length")
})

test_that("CHARACTERIZATION: an error inside the function leaks the sink", {

    # Output is suppressed with sink(tempfile()) and released by a bare sink()
    # at the end of the body. Any error in between skips the release, leaving
    # the caller's console redirected for the rest of the session.
    skip_if_not_installed("huge")

    set.seed(27)
    Y <- matrix(stats::rnorm(20 * 3), 20, 3)

    before <- sink.number()
    on.exit(while (sink.number() > before) sink(), add = TRUE)

    try(compute_omega(Y = Y, Y.rep = Y, R = 1, id = paste0("S", seq_len(20))),
        silent = TRUE)

    expect_gt(sink.number(), before)
})
