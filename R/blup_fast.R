#' Ridge BLUP for many responses sharing one design
#'
#' For each response y (a column of Y), univariate_blup fits
#' y = 1 b + X u + e with u ~ N(0, s2u I), e ~ N(0, s2e I) and the variances
#' estimated by REML: ridge regression BLUP, which is what
#' rrBLUP::mixed.solve(y, Z = X, K = diag(m)) computes (rrBLUP 4.6.3).
#' mixed.solve redoes its decompositions for every response, although they
#' depend only on X. Here .blup_prep() decomposes X once per design and
#' .blup_solve() finishes each response, in five steps.
#'
#' 1. Centre. Replace X by Xc = X - 1 xbar'. This only moves the intercept
#' (b becomes b + xbar'u), so u, the REML estimate and the fitted values do
#' not change. It makes Xc'1 = 0, which simplifies everything below.
#'
#' 2. Decompose. mixed.solve needs the spectrum of H0 = Z K Z' = Xc Xc' and
#' that of S H0 S, where S = I - 11'/n removes the intercept. Because
#' Xc'1 = 0, S H0 S = H0, so both are the nonzero eigenvalues d2 of Xc Xc'
#' plus zeros. .blup_prep() takes d2 from the smaller of Xc'Xc (n > m,
#' eigenvectors V) and Xc Xc' (n <= m, eigenvectors U), dropping eigenvalues
#' within rounding error of zero.
#'
#' 3. Estimate the variances. With lambda = s2e / s2u, mixed.solve's REML
#' objective depends on y only through omega = U'y, its projections on the
#' eigenvectors for d2, and tail = ||Sy||^2 - sum(omega^2), the part in the
#' n - 1 - length(d2) zero directions. .blup_solve() minimises that
#' objective with mixed.solve's optimize() call and bounds. (In the n > m
#' branch, U'y = diag(1 / d) V'Xc'y.)
#'
#' 4. Solve. H = H0 + lambda I has 1 as an eigenvector, so the GLS intercept
#' is mean(y), and the marker effects u = Xc'(Xc Xc' + lambda I)^-1 y are
#' \itemize{
#'   \item n > m: u = V diag(1 / (d2 + lambda)) V'Xc'y;
#'   \item n <= m: u = Xc'U diag(1 / (d2 + lambda)) U'y.
#' }
#'
#' 5. Return to the scale of X. u is the same for X and Xc; the intercept
#' for X is mean(y) - xbar'u, so X %*% u + intercept is the fitted value. A
#' monomorphic SNP gets a weight of exactly 0.
#'
#' Fitted values match mixed.solve to about 1e-6 relative, and weights to the
#' precision of optimize()'s search (about 1e-4 relative at worst). In its
#' "cholesky" branch (n > m + 1) mixed.solve adds 1e-6 to the diagonal of K;
#' that rescales H0 by 1 + 1e-6, which REML absorbs into lambda, so only the
#' search path differs.
#'
#' Up to isotwas 2.0.0 the prior was K = X'X / (m - 1) on uncentred dosages;
#' NEWS.md says why it was replaced.
#'
#' @noRd
NULL

#' Ridge BLUP weights and intercepts for every column of Y
#'
#' A response with missing values is fitted on its observed rows, as
#' mixed.solve drops them; responses with the same missing rows share one
#' decomposition. A response with no observed values gets weights 0 and an
#' NA intercept.
#'
#' @param Y n x G matrix of responses
#' @param X n x m matrix of SNP dosages, no missing values
#' @return list(weights = m x G matrix of marker effects, on the scale of X;
#'   intercept = length-G vector)
#' @noRd
.blup_fit <- function(Y, X) {
    Y <- as.matrix(Y)
    X <- as.matrix(X)
    if (anyNA(X)) {
        stop("univariate_blup: X has missing values; impute them first")
    }
    if (!anyNA(Y)) return(.blup_solve(.blup_prep(X), Y))

    out <- list(weights = matrix(0, ncol(X), ncol(Y)),
                intercept = rep(NA_real_, ncol(Y)))
    miss <- is.na(Y)
    pattern <- apply(miss, 2, function(z) paste(which(z), collapse = ","))
    for (cols in split(seq_len(ncol(Y)), pattern)) {
        ok <- !miss[, cols[1]]
        if (!any(ok)) next
        fit <- .blup_solve(.blup_prep(X[ok, , drop = FALSE]),
                           Y[ok, cols, drop = FALSE])
        out$weights[, cols] <- fit$weights
        out$intercept[cols] <- fit$intercept
    }
    out
}

#' @param X n x m design (SNP dosages), no missing values
#' @return a decomposition for .blup_solve()
#' @noRd
.blup_prep <- function(X) {
    X <- as.matrix(X)
    n <- nrow(X)
    m <- ncol(X)
    if (n < 1) stop("univariate_blup: no samples to fit")
    xbar <- colMeans(X)
    Xc <- sweep(X, 2, xbar)
    # monomorphic SNPs: centred to exactly 0, and given weight exactly 0 in
    # .blup_solve (in the X'X branch their eigenvector entries are only
    # zero to rounding error)
    mono <- colSums(X != rep(X[1, ], each = n)) == 0
    Xc[, mono] <- 0
    branch <- if (n > m) "XtX" else "XXt"
    e <- eigen(if (branch == "XtX") crossprod(Xc) else tcrossprod(Xc),
               symmetric = TRUE)
    # eigenvalues within rounding error of zero are the directions of H0
    # that mixed.solve gives theta = 0, or as good as
    tol <- max(abs(e$values)) * max(n, m) * .Machine$double.eps * 8
    keep <- which(e$values > tol)
    # Xc has rank n - 1 at most: in the XX' branch the direction 1 must
    # never be kept, even if rounding lifts it above tol
    if (branch == "XXt") keep <- utils::head(keep, n - 1)
    list(branch = branch, n = n, m = m, xbar = xbar, Xc = Xc, mono = mono,
         d2 = e$values[keep], E = e$vectors[, keep, drop = FALSE])
}

#' @param prep output of .blup_prep()
#' @param Y n x G matrix of responses, no missing values
#' @param bounds REML search interval for lambda, as mixed.solve
#' @return list(weights = m x G matrix, intercept = length-G vector), as
#'   .blup_fit()
#' @noRd
.blup_solve <- function(prep, Y, bounds = c(1e-09, 1e+09)) {
    Y <- as.matrix(Y)
    n.p <- prep$n - 1
    d2 <- prep$d2
    n_tail <- n.p - length(d2)
    # mixed.solve's f.REML, with the n_tail theta = 0 directions, whose
    # squared projections sum to tail, summed into one term
    f.REML <- function(lambda, n.p, theta, omega.sq, tail, n_tail) {
        n.p * log(sum(omega.sq / (theta + lambda)) + tail / lambda) +
            sum(log(theta + lambda)) + n_tail * log(lambda)
    }

    ybar <- colMeans(Y)
    Sy <- sweep(Y, 2, ybar)
    ss_Sy <- colSums(Sy^2)
    if (prep$branch == "XtX") {
        # V'Xc'y = diag(d) U'y
        t_y <- crossprod(prep$E, crossprod(prep$Xc, Sy))
        omega.sq <- t_y^2 / d2
    } else {
        t_y <- crossprod(prep$E, Sy)
        omega.sq <- t_y^2
    }

    u <- matrix(0, prep$m, ncol(Y))
    for (k in seq_len(ncol(Y))) {
        # no genetic variation, or a constant response: u = 0
        if (length(d2) == 0 || ss_Sy[k] == 0) next
        tail <- if (n_tail > 0) max(ss_Sy[k] - sum(omega.sq[, k]), 0) else 0
        lambda <- stats::optimize(f.REML, interval = bounds, n.p, d2,
                                  omega.sq[, k], tail, n_tail)$minimum
        r <- t_y[, k] / (d2 + lambda)
        u[, k] <- if (prep$branch == "XtX") {
            prep$E %*% r
        } else {
            crossprod(prep$Xc, prep$E %*% r)
        }
    }
    u[prep$mono, ] <- 0
    list(weights = u, intercept = ybar - drop(crossprod(prep$xbar, u)))
}
