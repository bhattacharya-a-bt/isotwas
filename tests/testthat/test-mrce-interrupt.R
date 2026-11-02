# The MRCE kernel is a .C() call; without an interrupt check, neither Ctrl-C
# nor setTimeLimit() can stop it. Takes about 3 s.

test_that("the blasso solver can be interrupted by setTimeLimit", {

    skip_on_cran()
    skip_if_not(is.loaded("blasso"), "compiled code not loaded")

    set.seed(101)
    p <- 260L; q <- 24L; n <- 400L
    S <- crossprod(matrix(stats::rnorm(n * p), n, p)) / n
    diag(S) <- diag(S) + 0.1
    M <- matrix(stats::rnorm(p * q), p, q)

    started <- Sys.time()
    setTimeLimit(elapsed = 3, transient = TRUE)
    on.exit(setTimeLimit(), add = TRUE)

    # A tolerance this tight with a penalty this small will not converge inside
    # maxit, so the loop runs until something stops it.
    expect_error(
        .C("blasso",
           B = as.double(matrix(0, p, q)), S = as.double(S), M = as.double(M),
           Om = as.double(diag(q)), soft = as.double(-M),
           p = p, q = q, n = n,
           lam = as.double(matrix(1e-8, p, q)),
           tol = as.double(1e-14), maxit = 100000L,
           totalit = integer(1), objective = as.double(0)))

    setTimeLimit()
    elapsed <- as.numeric(difftime(Sys.time(), started, units = "secs"))
    expect_lt(elapsed, 15)
})
