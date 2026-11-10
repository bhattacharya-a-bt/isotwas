# Go/no-go benchmark for compressing glmnet's input rows (not part of the
# package; benchmarks/ is in .Rbuildignore).
#
# Gaussian and mgaussian elastic-net paths depend on the data only through
# X'X, X'Y, Y'Y and n once X and Y are centred and X is standardised. So the
# n rows can be replaced by the p + q rows of the R factor of [X, Y], rescaled
# by sqrt(k / n), and glmnet returns the same path (see compress() below). This script measures, at
# one fold-sized subset (0.8 n rows):
#
#   t_plain  glmnet on the original rows, as isotwas calls it today
#   t_comp   glmnet on the compressed rows (gaussian: p + 1 rows per response)
#   t_qr     building the compressed rows (one pivoted QR of X, then Q'Y)
#   t_chol   the same rows from X'X by pivoted Cholesky
#
# and checks the two paths agree. It then projects the CV cost per gene for
# univariate_enet and multi_enet (6 subsets: the full data and 5 folds) with
# today's per-transcript / per-alpha fold draws ("legacy", one factorisation
# per distinct subset) and with one shared fold set ("shared", 6 in total),
# and for joinet (61 subsets, each shared across outputs either way).
#
# Usage:
#   Rscript benchmarks/bench_glmnet_compress.R [n] [p,p,...] [q] [reps]
# Defaults: n = 2000, p = 300,1000,2000,6000, q = 8, reps = 3.
# Run with single-threaded BLAS to match training jobs.

args <- commandArgs(trailingOnly = TRUE)
n    <- if (length(args) >= 1) as.integer(args[1]) else 2000L
ps   <- if (length(args) >= 2) as.integer(strsplit(args[2], ",")[[1]]) else
          c(300L, 1000L, 2000L, 6000L)
q    <- if (length(args) >= 3) as.integer(args[3]) else 8L
reps <- if (length(args) >= 4) as.integer(args[4]) else 3L

suppressPackageStartupMessages(library(glmnet))
cat(sprintf("glmnet %s | %s | n = %d, q = %d, reps = %d\n",
            packageVersion("glmnet"), R.version.string, n, q, reps))
# BLAS sanity check: a 2000 x 1000 crossprod is ~2 GFLOP. Single-threaded
# OpenBLAS manages ~10-40 GFLOP/s, reference BLAS ~1-3; well under 1 means
# oversubscribed threads (set OPENBLAS_NUM_THREADS=1 inside the allocation).
cat("BLAS:", extSoftVersion()[["BLAS"]], "| LAPACK:", La_library(), "\n")
cat("OPENBLAS_NUM_THREADS =", Sys.getenv("OPENBLAS_NUM_THREADS", "<unset>"),
    "| LD_PRELOAD =", Sys.getenv("LD_PRELOAD", "<unset>"), "\n")
local({
  set.seed(1); A <- matrix(rnorm(2000 * 1000), 2000)
  t <- system.time(crossprod(A))[["elapsed"]]
  cat(sprintf("crossprod 2000 x 1000: %.3fs = %.1f GFLOP/s\n\n",
              t, 2000 * 1000^2 / max(t, 1e-3) / 1e9))
})

# Genotypes with local LD: two AR(1) latent haplotypes thresholded at MAF.
sim_geno <- function(n, p, rho = 0.9) {
  maf <- runif(p, 0.05, 0.5)
  thr <- qnorm(1 - maf)
  hap <- function() {
    Z <- matrix(rnorm(n * p), n, p)
    for (j in seq_len(p)[-1]) Z[, j] <- rho * Z[, j - 1] + sqrt(1 - rho^2) * Z[, j]
    sweep(Z, 2, thr, ">") * 1
  }
  hap() + hap()
}

# Compressed rows from one pivoted QR of the centred X: X = Q R (columns
# unpivoted), so X'X = R'R for R scaled or not, and every response needs only
# Q'y. Rows for response j: [R, (Q'y_j)[1:p]] plus one row carrying the norm of
# the residual part (Q'y_j)[-(1:p)]. Rows for mgaussian: [R, Q'Y[1:p, ]] over
# [0, R_e], with R_e the q x q factor of the residual part. Both reproduce
# X'X, X'Y and Y'Y exactly; rows are rescaled by sqrt(k / n) so glmnet's 1/k
# weighting equals the original 1/n.
compress <- function(X, Y) {
  n <- nrow(X); p <- ncol(X)
  xm <- colMeans(X); Xc <- sweep(X, 2, xm)
  xs <- sqrt(colSums(Xc^2) / n); xs[xs == 0] <- 1
  qrx <- qr(Xc, LAPACK = TRUE)
  R <- qr.R(qrx)[, order(qrx$pivot), drop = FALSE]
  QtY <- qr.qty(qrx, sweep(Y, 2, colMeans(Y)))
  E <- QtY[-(1:p), , drop = FALSE]
  list(n = n, p = p, xs = xs, R = R, QtY = QtY[1:p, , drop = FALSE],
       enorm = sqrt(colSums(E^2)), Re = qr.R(qr(E)))
}

# Same rows from the Gram matrix: pivoted Cholesky of X'X (rank r, trailing
# rows dropped) gives an r x p R with R'R = X'X; Q'Y = R^{-T} X'Y on the
# pivoted r columns; the residual Gram Y'Y - (Q'Y)'(Q'Y) gives the last rows.
# Cheaper than the pivoted QR (BLAS-3 crossprod + p^3/3), and fold Grams can be
# downdated from the full one, at the cost of working with X'X (as glmnet's
# own covariance mode does).
compress_chol <- function(X, Y) {
  n <- nrow(X); p <- ncol(X)
  Xc <- sweep(X, 2, colMeans(X)); Yc <- sweep(Y, 2, colMeans(Y))
  xs <- sqrt(colSums(Xc^2) / n); xs[xs == 0] <- 1
  C <- suppressWarnings(chol(crossprod(Xc), pivot = TRUE))
  r <- attr(C, "rank"); piv <- attr(C, "pivot")
  C <- C[seq_len(r), , drop = FALSE]
  QtY <- backsolve(C[, seq_len(r), drop = FALSE],
                   crossprod(Xc[, piv[seq_len(r)], drop = FALSE], Yc),
                   transpose = TRUE)
  E <- crossprod(Yc) - crossprod(QtY)
  Ce <- suppressWarnings(chol(E, pivot = TRUE))
  if (attr(Ce, "rank") < ncol(E)) Ce[-seq_len(attr(Ce, "rank")), ] <- 0
  list(n = n, p = p, xs = xs, R = C[, order(piv), drop = FALSE], QtY = QtY,
       enorm = sqrt(pmax(diag(E), 0)),
       Re = Ce[, order(attr(Ce, "pivot")), drop = FALSE])
}

gaussian_rows <- function(cx, j, standardize = TRUE) {
  R <- if (standardize) sweep(cx$R, 2, cx$xs, "/") else cx$R
  s <- sqrt((nrow(R) + 1) / cx$n)
  list(X = rbind(R, 0) * s, y = c(cx$QtY[, j], cx$enorm[j]) * s)
}

mgaussian_rows <- function(cx, standardize = FALSE) {
  R <- if (standardize) sweep(cx$R, 2, cx$xs, "/") else cx$R
  q <- ncol(cx$QtY)
  s <- sqrt((nrow(R) + q) / cx$n)
  list(X = rbind(R, matrix(0, q, cx$p)) * s, Y = rbind(cx$QtY, cx$Re) * s)
}

med_time <- function(expr) {
  expr <- substitute(expr); env <- parent.frame()
  median(replicate(reps, system.time(eval(expr, env))[["elapsed"]]))
}

path_diff <- function(B1, B2) {
  L <- min(ncol(B1), ncol(B2))
  max(abs(B1[, seq_len(L), drop = FALSE] - B2[, seq_len(L), drop = FALSE]))
}

rows <- list()
for (p in ps) {
  set.seed(p)
  nf <- round(0.8 * n)
  X <- sim_geno(nf, p)
  X <- X[, apply(X, 2, sd) > 0, drop = FALSE]
  p <- ncol(X)
  causal <- sample(p, min(p, 10))
  Y <- sapply(seq_len(q), function(j)
    drop(X[, causal] %*% rnorm(length(causal), 0, 0.2)) + rnorm(nf))
  can <- p + q < nf

  t_qr <- if (can) med_time(cx <- compress(X, Y)) else NA
  t_chol <- if (can) med_time(cc <- compress_chol(X, Y)) else NA

  for (a in c(0, 0.5, 1)) {
    t_plain <- med_time(f1 <- glmnet(X, Y[, 1], alpha = a))
    t_comp <- diff <- diff_chol <- NA
    if (can) {
      g <- gaussian_rows(cx, 1)
      t_comp <- med_time(f2 <- glmnet(g$X, g$y, alpha = a, intercept = FALSE,
                                      standardize = FALSE))
      diff <- path_diff(as.matrix(f1$beta), as.matrix(f2$beta) / cx$xs)
      g <- gaussian_rows(cc, 1)
      f3 <- glmnet(g$X, g$y, alpha = a, intercept = FALSE, standardize = FALSE)
      diff_chol <- path_diff(as.matrix(f1$beta), as.matrix(f3$beta) / cc$xs)
    }
    rows[[length(rows) + 1]] <- data.frame(
      p = p, family = "gaussian", alpha = a, t_plain = t_plain, t_comp = t_comp,
      t_qr = t_qr, t_chol = t_chol, maxdiff = diff, maxdiff_chol = diff_chol)

    t_plain <- med_time(f1 <- glmnet(X, Y, family = "mgaussian", alpha = a,
                                     standardize = FALSE))
    t_comp <- diff <- diff_chol <- NA
    if (can) {
      mg <- mgaussian_rows(cx)
      t_comp <- med_time(f2 <- glmnet(mg$X, mg$Y,
                                      family = "mgaussian", alpha = a,
                                      intercept = FALSE, standardize = FALSE))
      diff <- max(sapply(seq_len(q), function(j)
        path_diff(as.matrix(f1$beta[[j]]), as.matrix(f2$beta[[j]]))))
      mg <- mgaussian_rows(cc)
      f3 <- glmnet(mg$X, mg$Y, family = "mgaussian", alpha = a,
                   intercept = FALSE, standardize = FALSE)
      diff_chol <- max(sapply(seq_len(q), function(j)
        path_diff(as.matrix(f1$beta[[j]]), as.matrix(f3$beta[[j]]))))
    }
    rows[[length(rows) + 1]] <- data.frame(
      p = p, family = "mgaussian", alpha = a, t_plain = t_plain, t_comp = t_comp,
      t_qr = t_qr, t_chol = t_chol, maxdiff = diff, maxdiff_chol = diff_chol)
  }
  r <- do.call(rbind, rows[(length(rows) - 5):length(rows)])
  for (i in seq_len(nrow(r)))
    cat(sprintf("p=%5d %-9s alpha=%.1f  plain %7.3fs  comp %7.3fs  (%5.1fx)  qr %6.3fs  chol %6.3fs  maxdiff qr %.1e chol %.1e\n",
                r$p[i], r$family[i], r$alpha[i], r$t_plain[i], r$t_comp[i],
                r$t_plain[i] / r$t_comp[i], r$t_qr[i], r$t_chol[i],
                r$maxdiff[i], r$maxdiff_chol[i]))
}
res <- do.call(rbind, rows)

# Projected CV cost per gene (seconds), from the fold-sized timings above,
# using the cheaper of the two factorisations.
proj <- do.call(rbind, lapply(unique(res$p), function(p) {
  g <- res[res$p == p & res$family == "gaussian" & res$alpha == 0.5, ]
  m <- res[res$p == p & res$family == "mgaussian", ]
  tg <- g$t_plain; tgc <- g$t_comp; tq <- min(g$t_qr, g$t_chol)
  tm <- sum(m$t_plain); tmc <- sum(m$t_comp); tqu <- tq
  data.frame(
    p = p,
    uni_plain = 6 * q * tg,
    uni_legacy = 6 * q * tgc + (5 * q + 1) * tq,
    uni_shared = 6 * q * tgc + 6 * tq,
    multi_plain = 6 * tm,
    multi_legacy = 6 * tmc + 16 * tqu,
    multi_shared = 6 * tmc + 6 * tqu,
    joinet_plain = 61 * q * tg,
    joinet_comp = 61 * q * tgc + 61 * tq)
}))
cat("\nProjected CV seconds per gene (q =", q, "):\n")
print(format(proj, digits = 3), row.names = FALSE)
cat("\nSpeed-up vs plain:\n")
print(format(data.frame(
  p = proj$p,
  uni_legacy = proj$uni_plain / proj$uni_legacy,
  uni_shared = proj$uni_plain / proj$uni_shared,
  multi_legacy = proj$multi_plain / proj$multi_legacy,
  multi_shared = proj$multi_plain / proj$multi_shared,
  joinet = proj$joinet_plain / proj$joinet_comp), digits = 3), row.names = FALSE)

cat("\nTSV\n")
write.table(res, stdout(), sep = "\t", quote = FALSE, row.names = FALSE)
