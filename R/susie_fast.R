#' susieR::susie as univariate_susie calls it, for many responses
#'
#' univariate_susie fits susieR::susie(X, y, estimate_residual_variance =
#' TRUE, estimate_prior_variance = FALSE, scaled_prior_variance = 0.1) with
#' the other arguments at their defaults (L = min(10, p), intercept and
#' standardize TRUE, max_iter = 100, tol = 1e-4, convergence on the ELBO).
#' .susie_prep() computes what depends on X alone, once per training set;
#' .susie_fit() runs IBSS for one response.
#'
#' The algorithm is susieR 0.16.6's, step for step: the same order of
#' single-effect updates, the same log Bayes factors, posterior moments, KL
#' terms, ELBO, convergence test (0 <= ELBO change < tol), residual variance
#' update and bounds, and the same trimming of null effects. What changes is
#' how the products with the standardised X are formed:
#' \itemize{
#'   \item X'R for the residual R is taken from X'y and a cached X'X b_l per
#'     effect when n >= p (a p x p product instead of two n x p ones), or
#'     from a cached X b_l per effect otherwise;
#'   \item the KL term's R'X E[b] is formed as (X'R)'E[b], and the ELBO's
#'     ||X b_l||^2 from the same caches.
#' }
#' The estimates agree with susieR to rounding error, as long as no ELBO
#' change falls within rounding error of the tolerance.
#'
#' @noRd
NULL

#' @param X n x p matrix, no missing values
#' @param gram logical, work from the standardised Gram matrix
#' @param XtX crossprod(X), if already known
#' @noRd
.susie_prep <- function(X, gram = nrow(X) >= ncol(X) && .use_gram(ncol(X)),
                        XtX = NULL) {
    X <- as.matrix(X)
    storage.mode(X) <- "double"
    n <- nrow(X)
    p <- ncol(X)
    # susieR's compute_colstats with center = scale = TRUE; X reaches susie
    # as a sparse Matrix in univariate_susie, hence the sparse-matrix sd
    cm <- colMeans(X)
    csd_raw <- sqrt((colMeans(X^2) - cm^2) * n / (n - 1))
    csd <- csd_raw
    csd[csd == 0] <- 1
    d <- n * cm^2 + (n - 1) * csd_raw^2
    d <- (d - n * cm^2) / csd^2
    G <- NULL
    if (gram) {
        G <- (if (is.null(XtX)) crossprod(X) else XtX) - n * tcrossprod(cm)
        G <- G / csd
        G <- t(t(G) / csd)
    }
    list(X = X, n = n, p = p, cm = cm, csd = csd, d = d, G = G)
}

#' @param prep output of .susie_prep()
#' @param y numeric response
#' @return list(intercept, beta, converged, niter) with
#'   c(intercept, beta) equal to coef() of the susieR fit
#' @noRd
.susie_fit <- function(prep, y, scaled_prior_variance = 0.1,
                       max_iter = 100, tol = 1e-4, prior_tol = 1e-9) {
    X <- prep$X
    n <- prep$n
    p <- prep$p
    cm <- prep$cm
    csd <- prep$csd
    d <- prep$d
    G <- prep$G
    L <- min(10, p)
    y <- as.numeric(y)
    if (anyNA(y)) stop("Input y must not contain missing values.")

    lowerbound <- stats::var(y) / 1e4
    mean_y <- mean(y)
    y <- y - mean_y
    var_y <- stats::var(y)
    if (!(var_y > 0)) {
        stop("Residual variance sigma2 must be positive (is your var(Y) zero?).")
    }
    sigma2 <- var_y
    V <- scaled_prior_variance * var_y
    lpi <- log(rep(1 / p, p) + sqrt(.Machine$double.eps))

    # standardised products: Xs = t((t(X) - cm) / csd)
    Xty <- drop(crossprod(X, y)) / csd - cm / csd * sum(y)
    Xb <- function(b) drop(X %*% (b / csd)) - sum(cm * b / csd)

    alpha <- matrix(1 / p, L, p)
    mu <- matrix(0, L, p)
    mu2 <- matrix(0, L, p)
    KL <- rep(NA_real_, L)
    Eb <- matrix(0, L, p)                      # alpha * mu
    if (!is.null(G)) {
        Gb <- matrix(0, p, L)                  # Xs'Xs b_l
        GB <- numeric(p)                       # Xs'Xs sum_l b_l
    } else {
        Xbl <- matrix(0, n, L)                 # Xs b_l
        Xr <- numeric(n)
    }
    yty <- sum(y^2)

    prev_elbo <- -Inf
    converged <- FALSE
    for (iter in seq_len(max_iter)) {
        for (l in seq_len(L)) {
            if (!is.null(G)) {
                XtR <- Xty - (GB - Gb[, l])
            } else {
                R <- y - (Xr - Xbl[, l])
                XtR <- drop(crossprod(X, R)) / csd - cm / csd * sum(R)
            }
            betahat <- XtR / d
            shat2 <- sigma2 / d

            # gaussian_ser_lbf and lbf_stabilization
            shat2_safe <- pmax(shat2, .Machine$double.eps)
            lbf <- -0.5 * log(1 + V / shat2_safe) +
                0.5 * betahat^2 * V / (shat2_safe * (V + shat2_safe))
            lbf[!is.finite(betahat) | !is.finite(shat2)] <- 0
            lpo <- lbf + lpi
            inf_idx <- is.infinite(shat2)
            lpo[inf_idx] <- lpi[inf_idx]
            # compute_posterior_weights
            mx <- max(lpo)
            w <- exp(lpo - mx)
            sw <- sum(w)
            a_l <- w / sw
            lbf_model <- log(sw) + mx
            # gaussian_ser_moments
            post_var <- V * shat2 / (V + shat2)
            post_mean <- post_var / shat2 * betahat
            no_info <- !is.finite(betahat) | !is.finite(shat2)
            post_var[no_info] <- 0
            post_mean[no_info] <- 0

            alpha[l, ] <- a_l
            mu[l, ] <- post_mean
            mu2[l, ] <- post_var + post_mean^2
            eb <- a_l * post_mean
            Eb[l, ] <- eb
            # compute_kl.individual, with the R'R and log 2 pi sigma2 terms
            # that cancel between its two halves removed
            KL[l] <- -lbf_model +
                (2 * sum(XtR * eb) - sum(d * a_l * mu2[l, ])) / (2 * sigma2)

            if (!is.null(G)) {
                gb <- drop(G %*% eb)
                GB <- GB - Gb[, l] + gb
                Gb[, l] <- gb
            } else {
                xb <- Xb(eb)
                Xr <- Xr - Xbl[, l] + xb
                Xbl[, l] <- xb
            }
        }

        # get_ER2.individual
        B <- colSums(Eb)
        Eb2_d <- sum((alpha * mu2) %*% d)
        if (!is.null(G)) {
            rss <- yty - 2 * sum(Xty * B) + sum(B * GB)
            xb2 <- sum(Eb * t(Gb))
        } else {
            rss <- sum((y - Xr)^2)
            xb2 <- sum(Xbl^2)
        }
        ER2 <- rss + Eb2_d - xb2
        elbo <- -n / 2 * log(2 * pi * sigma2) - ER2 / (2 * sigma2) - sum(KL)
        if (is.infinite(elbo)) stop("get_objective() produced an infinite ELBO value")

        if (iter > 1) {
            diff <- elbo - prev_elbo
            converged <- !is.na(diff) && diff >= 0 && diff < tol
        }
        prev_elbo <- elbo
        if (converged) break

        # update_model_variance.default
        sigma2 <- ER2 / n
        if (sigma2 < 0) stop("est_residual_variance() failed: the estimated value is negative")
        sigma2 <- min(max(sigma2, lowerbound), Inf)
    }
    if (!converged) {
        warning("IBSS algorithm did not converge in ", max_iter, " iterations!",
                call. = FALSE)
    }

    # trim_null_effects: V is fixed and equal for every effect
    if (V < prior_tol) {
        alpha[] <- 1 / p
        mu[] <- 0
    }
    b <- colSums(alpha * mu) / csd
    list(intercept = mean_y - sum(cm * b),
         beta = b,
         converged = converged,
         niter = iter)
}

#' Coefficients (without intercept) of univariate_susie's fit for each column
#' of Y
#' @return p x ncol(Y) matrix
#' @noRd
.susie_coef <- function(X, Y, XtX = NULL) {
    prep <- .susie_prep(X, XtX = XtX)
    Y <- as.matrix(Y)
    vapply(seq_len(ncol(Y)), function(g) .susie_fit(prep, Y[, g])$beta,
           numeric(ncol(prep$X)))
}
