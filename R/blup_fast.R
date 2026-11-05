#' rrBLUP::mixed.solve for many responses sharing one design
#'
#' univariate_blup calls rrBLUP::mixed.solve(y, Z = X, K = X'X / (p - 1))
#' once per transcript, with the default intercept as the only fixed effect
#' and REML estimation. Everything mixed.solve decomposes depends on Z and K
#' only, so .blup_prep() does that work once per design and .blup_solve()
#' finishes each response: the REML search, beta and u.
#'
#' mixed.solve (rrBLUP 4.6.3) needs two spectra: phi, of H0 = Z K Z', and
#' theta with the projections omega = Q'y, of S H0 S where S = I - 11'/n.
#' Because K is a function of Z'Z, both follow from one eigendecomposition of
#' the smaller Gram matrix and one of a diagonal-minus-rank-one matrix:
#' \itemize{
#'   \item n > m + 1 (mixed.solve's "cholesky" branch, where K carries a
#'     1e-6 ridge): Z'Z = R diag(mu) R', so K = R diag(Lambda) R' with
#'     Lambda = mu / (m - 1) + 1e-6, and H0 has nonzero eigenvalues
#'     a = Lambda mu with eigenvectors Z R Lambda^(1/2) a^(-1/2). Its centred
#'     version S H0 S has the nonzero spectrum of diag(a) - v v' / n with
#'     v = Lambda^(1/2) R' Z'1.
#'   \item n <= m + 1 (the "eigen" branch, K without the ridge): ZZ' =
#'     U diag(nu) U', so H0 = U diag(nu^2 / (m - 1)) U', and S H0 S has the
#'     spectrum of diag(phi) - w w' / n with w = phi^(1/2) U'1.
#' }
#' H^-1 = (H0 + lambda I)^-1 is applied in the same bases. The REML
#' objective, its bounds and the optimize() call are mixed.solve's; the
#' objective is evaluated on the nonzero part of theta, with the rest (exact
#' zeros in mixed.solve, or eigenvalues below rounding error) summed into
#' one term. The estimates match mixed.solve to the precision of its own
#' REML search.
#'
#' @noRd
NULL

#' Marker effects of univariate_blup for every column of Y
#'
#' Equal to sapply over columns of rrBLUP::mixed.solve(Y[, g], Z = Z,
#' K = crossprod(Z) / (ncol(Z) - 1))$u, which it falls back to when there
#' are missing values.
#'
#' @return ncol(Z) x ncol(Y) matrix
#' @noRd
.blup_u <- function(Y, Z) {
    Y <- as.matrix(Y)
    Z <- as.matrix(Z)
    if (anyNA(Y) || anyNA(Z)) {
        K <- crossprod(Z) / (ncol(Z) - 1)
        return(sapply(seq_len(ncol(Y)), function(g)
            as.numeric(rrBLUP::mixed.solve(Y[, g], Z = Z, K = K)$u)))
    }
    .blup_solve(.blup_prep(Z), Y)
}

#' Eigenvalues and eigenvectors of diag(a) - v v' / n
#'
#' Eigenvalues within rounding error of zero are dropped: they stand for the
#' directions of S H0 S that mixed.solve gives theta = 0, or as good as.
#' @noRd
.blup_centred_spectrum <- function(a, v, n) {
    T <- -tcrossprod(v) / n
    diag(T) <- diag(T) + a
    e <- eigen(T, symmetric = TRUE)
    keep <- e$values > max(abs(e$values)) * length(a) * .Machine$double.eps * 8
    list(theta = e$values[keep], V = e$vectors[, keep, drop = FALSE])
}

#' @param Z n x m design (SNP dosages), no missing values
#' @return a decomposition for .blup_solve()
#' @noRd
.blup_prep <- function(Z) {
    Z <- as.matrix(Z)
    n <- nrow(Z)
    m <- ncol(Z)
    s <- colSums(Z)
    if (n > m + 1) {
        e <- eigen(crossprod(Z), symmetric = TRUE)
        R <- e$vectors
        mu <- e$values
        Lambda <- mu / (m - 1) + 1e-06
        a <- Lambda * mu
        Rs <- drop(crossprod(R, s))
        cs <- .blup_centred_spectrum(a, sqrt(Lambda) * Rs, n)
        list(branch = "cholesky", n = n, m = m, Z = Z, R = R,
             Lambda = Lambda, phi = a, Rs = Rs,
             theta = cs$theta, V2 = cs$V)
    } else {
        e <- eigen(tcrossprod(Z), symmetric = TRUE)
        U <- e$vectors
        nu <- e$values
        phi <- nu^2 / (m - 1)
        U1 <- colSums(U)
        cs <- .blup_centred_spectrum(phi, sqrt(pmax(phi, 0)) * U1, n)
        list(branch = "eigen", n = n, m = m, Z = Z, U = U,
             nu = nu, phi = phi, Rs = U1,
             theta = cs$theta, V2 = cs$V)
    }
}

#' @param prep output of .blup_prep()
#' @param Y n x G matrix of responses, no missing values
#' @param bounds REML search interval, as mixed.solve
#' @return m x G matrix of marker effects u
#' @noRd
.blup_solve <- function(prep, Y, bounds = c(1e-09, 1e+09)) {
    Y <- as.matrix(Y)
    n <- prep$n
    m <- prep$m
    n.p <- n - 1
    theta <- prep$theta
    n_tail <- n.p - length(theta)
    phi <- prep$phi
    # mixed.solve's f.REML, with the theta = 0 directions summed into one term
    f.REML <- function(lambda, n.p, theta, omega.sq, tail) {
        n.p * log(sum(omega.sq / (theta + lambda)) + tail / lambda) +
            sum(log(theta + lambda)) + n_tail * log(lambda)
    }

    Sy <- sweep(Y, 2, colMeans(Y))
    if (prep$branch == "cholesky") {
        # rotated Z'y, and Z'Sy in the basis of the centred spectrum
        t_y <- crossprod(prep$R, crossprod(prep$Z, Y))
        t_Sy <- crossprod(prep$R, crossprod(prep$Z, Sy))
        g <- crossprod(prep$V2, sqrt(prep$Lambda) * t_Sy)
    } else {
        t_y <- crossprod(prep$U, Y)
        g <- crossprod(prep$V2, sqrt(pmax(phi, 0)) * crossprod(prep$U, Sy))
    }
    ss_Sy <- colSums(Sy^2)
    one_y <- colSums(Y)

    out <- matrix(0, m, ncol(Y))
    for (k in seq_len(ncol(Y))) {
        omega.sq <- g[, k]^2 / theta
        tail <- max(ss_Sy[k] - sum(omega.sq), 0)
        lambda <- stats::optimize(f.REML, interval = bounds, n.p, theta,
                                  omega.sq, tail)$minimum
        if (prep$branch == "cholesky") {
            # 1'H^-1 v = (1'v - Rs' diag(Lambda / (lambda + phi)) R'Z'v) / lambda
            # Z'H^-1 v = R diag(1 / (lambda + phi)) R'Z'v
            wgt <- prep$Lambda / (lambda + phi)
            Rs <- prep$Rs
            oneHy <- (one_y[k] - sum(Rs * wgt * t_y[, k])) / lambda
            oneH1 <- (n - sum(Rs * wgt * Rs)) / lambda
            beta <- oneHy / oneH1
            r <- (t_y[, k] - beta * Rs) / (lambda + phi)
            # u = K Z'w = R diag(Lambda) R'Z'w
            out[, k] <- prep$R %*% (prep$Lambda * r)
        } else {
            # H^-1 = U diag(1 / (phi + lambda)) U'
            U1 <- prep$Rs
            dinv <- 1 / (phi + lambda)
            beta <- sum(U1 * dinv * t_y[, k]) / sum(U1^2 * dinv)
            r <- (t_y[, k] - beta * U1) * dinv
            # u = K Z'w = Z'(ZZ') w / (m - 1) = Z'U diag(nu) U'w / (m - 1)
            out[, k] <- crossprod(prep$Z, prep$U %*% (prep$nu * r)) / (m - 1)
        }
    }
    out
}
