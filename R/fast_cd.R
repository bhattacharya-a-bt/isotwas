#' Compiled row-wise coordinate descent
#'
#' Helpers around the C solver in src/rowcd.c, which multivariate_sgl,
#' multivariate_mtlasso and multivariate_graph_reg share. The solver makes the
#' same updates in the same order as the R reference solvers (.fit_sgl_r,
#' .fit_mtl_l21_r, .fit_graph_reg_r) and stops by the same rule, so its
#' estimates match theirs to rounding error.
#'
#' @name fast_cd
#' @keywords internal
#' @noRd
NULL

#' Sufficient statistics for one training set
#'
#' @param X,Y matrices already centred the way the method centres them
#' @param idx integer, training rows
#' @param G_full,XtY_full crossprod(X) and crossprod(X, Y) on all rows, or
#'   NULL. The training-set versions are formed by subtracting the held-out
#'   rows, which costs a fraction of a fresh crossprod.
#'
#' @return list(G, XtY, X, Y, d, n): G/XtY in covariance mode, X/Y otherwise
#' @noRd
.cd_data <- function(X, Y, idx = NULL, G_full = NULL, XtY_full = NULL) {
    if (is.null(idx)) {
        Xt <- X
        Yt <- Y
        G <- G_full
        XtY <- XtY_full
    } else {
        Xt <- X[idx, , drop = FALSE]
        Yt <- Y[idx, , drop = FALSE]
        if (!is.null(G_full)) {
            Xo <- X[-idx, , drop = FALSE]
            G <- G_full - crossprod(Xo)
            XtY <- XtY_full - crossprod(Xo, Y[-idx, , drop = FALSE])
        } else {
            G <- NULL
            XtY <- NULL
        }
    }
    n <- nrow(Xt)
    list(G = G,
         XtY = XtY,
         X = if (is.null(G)) Xt else NULL,
         Y = if (is.null(G)) Yt else NULL,
         d = colSums(Xt^2) / n,
         n = n)
}

#' X %*% B over the rows of B that are not all zero
#'
#' The skipped rows add exact zeros, so the product is X %*% B.
#' @noRd
.xb_sparse <- function(X, B) {
    nz <- which(rowSums(B != 0) > 0)
    X[, nz, drop = FALSE] %*% B[nz, , drop = FALSE]
}

#' Run the compiled solver on prepared sufficient statistics
#'
#' @param cd output of .cd_data()
#' @param type 0 sparse group lasso, 1 L21, 2 elastic net, 3 graph-regularised
#' @param pars numeric c(l1, l2, lg)
#' @param B0 p x q warm start, or NULL for zeros
#' @param V,ev eigenvectors and eigenvalues of the Laplacian (type 3)
#'
#' @return p x q coefficient matrix
#' @noRd
.cd_solve <- function(cd, type, pars, B0 = NULL, max_iter, tol,
                      V = NULL, ev = NULL) {
    p <- length(cd$d)
    q <- if (!is.null(cd$XtY)) ncol(cd$XtY) else ncol(cd$Y)
    if (is.null(B0)) B0 <- matrix(0, p, q)
    storage.mode(B0) <- "double"
    .Call("isotwas_rowcd",
          cd$G, cd$XtY, cd$X, cd$Y, B0,
          as.double(cd$d), as.double(cd$n), as.integer(type),
          as.double(pars),
          if (is.null(V)) NULL else as.double(V),
          if (is.null(ev)) NULL else as.double(ev),
          as.integer(max_iter), as.double(tol),
          PACKAGE = "isotwas")$B
}
