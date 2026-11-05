# Simulated dosages and expression for the tests that check a faster solver
# against the one it replaces (test-fast-*.R).

sim_geno <- function(n, p, seed) {
  withr::with_seed(seed, {
    X <- sapply(stats::runif(p, .05, .5), function(f) stats::rbinom(n, 2, f))
    colnames(X) <- paste0("rs", seq_len(p))
    X
  })
}

sim_expr <- function(X, q, seed, n_causal = 5) {
  withr::with_seed(seed, {
    B <- matrix(0, ncol(X), q)
    B[sample(ncol(X), n_causal), ] <- stats::rnorm(n_causal * q, 0, .4)
    Y <- X %*% B + matrix(stats::rnorm(nrow(X) * q), nrow(X), q)
    colnames(Y) <- paste0("tx", seq_len(q))
    Y
  })
}
