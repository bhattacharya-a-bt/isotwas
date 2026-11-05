/*
 * Row-wise block coordinate descent shared by multivariate_sgl,
 * multivariate_mtlasso and multivariate_graph_reg.
 *
 * These reproduce the R solvers .fit_sgl_r, .fit_mtl_l21_r and
 * .fit_graph_reg_r update for update: the same cyclic order over SNPs, the
 * same row update and the same stopping rule (largest absolute coefficient
 * change in a sweep below tol). Only the way the partial-residual gradient
 * is obtained differs, so results agree to rounding error.
 *
 * Two ways to get the gradient X_j' R / n for row j:
 *   covariance mode  keeps C = X'Y - X'X B (p x q) and updates it with one
 *                    column of X'X when a row changes: O(q) per row visited,
 *                    O(p q) per row that changes.
 *   naive mode       keeps the residual R = Y - X B (n x q): O(n q) per row
 *                    visited. Used when X'X would be too large to hold.
 *
 * Row updates, with g the partial-residual gradient (X_j' R_{-j} / n),
 * d = X_j' X_j / n:
 *   type 0  sparse group lasso   s = soft(g, l1); new = max(0, 1 - l2 /
 *                                (d ||s|| + 1e-10)) s / d
 *   type 1  L21 multi-task lasso new = (1 - l1 / ||g||) g / d if ||g|| > l1
 *   type 2  elastic net          new = soft(g, l1) / (d + l2)
 *   type 3  graph-regularised    u = (d + l2) I + lg L)^{-1} g, solved with
 *                                L = V diag(ev) V'; new = soft(u, l1 / (d +
 *                                l2 + lg))
 */

#include <R.h>
#include <Rinternals.h>
#include <R_ext/Utils.h>
#include <math.h>
#include <string.h>

static inline double soft1(double x, double t)
{
    double a = fabs(x) - t;
    if (a > 0.0) return (x > 0.0) ? a : ((x < 0.0) ? -a : 0.0);
    return 0.0;
}

/* Compute the new row in out[] from gradient g[] (length q). */
static void row_update(int type, int q, const double *g, double d,
                       double l1, double l2, double lg,
                       const double *V, const double *ev,
                       double *work, double *out)
{
    int k, m;
    if (type == 0) {
        double ns = 0.0;
        for (k = 0; k < q; k++) {
            out[k] = soft1(g[k], l1);
            ns += out[k] * out[k];
        }
        ns = sqrt(ns);
        if (ns > 0.0) {
            double sh = 1.0 - l2 / (d * ns + 1e-10);
            if (sh < 0.0) sh = 0.0;
            for (k = 0; k < q; k++) out[k] = sh * out[k] / d;
        } else {
            for (k = 0; k < q; k++) out[k] = 0.0;
        }
    } else if (type == 1) {
        double gn = 0.0;
        for (k = 0; k < q; k++) gn += g[k] * g[k];
        gn = sqrt(gn);
        if (gn > l1) {
            double sh = 1.0 - l1 / gn;
            for (k = 0; k < q; k++) out[k] = sh * g[k] / d;
        } else {
            for (k = 0; k < q; k++) out[k] = 0.0;
        }
    } else if (type == 2) {
        double den = d + l2;
        for (k = 0; k < q; k++) out[k] = soft1(g[k], l1) / den;
    } else {
        /* work[0..q-1] = V' g / (d + l2 + lg ev) */
        for (m = 0; m < q; m++) {
            double s = 0.0;
            const double *vm = V + (size_t)m * q;
            for (k = 0; k < q; k++) s += vm[k] * g[k];
            work[m] = s / (d + l2 + lg * ev[m]);
        }
        double thr = l1 / (d + l2 + lg);
        for (k = 0; k < q; k++) {
            double u = 0.0;
            for (m = 0; m < q; m++) u += V[k + (size_t)m * q] * work[m];
            out[k] = soft1(u, thr);
        }
    }
}

/*
 * isotwas_rowcd(G, XtY, X, Y, B0, n, type, pars, V, ev, maxit, tol)
 *   G, XtY  : X'X (p x p) and X'Y (p x q), unscaled; covariance mode
 *   X, Y    : n x p and n x q; naive mode (used when G is NULL)
 *   B0      : p x q warm start
 *   d       : length-p X_j'X_j / n, the curvature the R solvers use
 *   pars    : c(l1, l2, lg)
 * Returns list(B, iter).
 */
SEXP isotwas_rowcd(SEXP sG, SEXP sXtY, SEXP sX, SEXP sY, SEXP sB0, SEXP sD,
                   SEXP sN, SEXP sType, SEXP sPars, SEXP sV, SEXP sEv,
                   SEXP sMaxit, SEXP sTol)
{
    const int cov = !isNull(sG);
    const int *dim = INTEGER(getAttrib(sB0, R_DimSymbol));
    const int pp = dim[0], q = dim[1];
    const double n = asReal(sN);
    const int type = asInteger(sType);
    const double l1 = REAL(sPars)[0], l2 = REAL(sPars)[1], lg = REAL(sPars)[2];
    const int maxit = asInteger(sMaxit);
    const double tol = asReal(sTol);
    const double *d = REAL(sD);
    const double *V = (type == 3) ? REAL(sV) : NULL;
    const double *ev = (type == 3) ? REAL(sEv) : NULL;

    SEXP sB = PROTECT(duplicate(sB0));
    double *B = REAL(sB);
    double *g = (double *) R_alloc(q, sizeof(double));
    double *bn = (double *) R_alloc(q, sizeof(double));
    double *work = (double *) R_alloc(q, sizeof(double));
    double *delta = (double *) R_alloc(q, sizeof(double));

    int nobs = 0;
    double *C = NULL, *Rr = NULL;
    const double *G = NULL, *X = NULL;
    size_t P = (size_t) pp;

    if (cov) {
        G = REAL(sG);
        C = (double *) R_alloc(P * q, sizeof(double));
        memcpy(C, REAL(sXtY), sizeof(double) * P * q);
        for (int j = 0; j < pp; j++) {
            for (int k = 0; k < q; k++) {
                double b = B[j + P * k];
                if (b != 0.0) {
                    const double *gj = G + P * j;
                    double *ck = C + P * k;
                    for (int a = 0; a < pp; a++) ck[a] -= gj[a] * b;
                }
            }
        }
    } else {
        X = REAL(sX);
        nobs = nrows(sX);
        size_t N = (size_t) nobs;
        Rr = (double *) R_alloc(N * q, sizeof(double));
        memcpy(Rr, REAL(sY), sizeof(double) * N * q);
        for (int j = 0; j < pp; j++) {
            for (int k = 0; k < q; k++) {
                double b = B[j + P * k];
                if (b != 0.0) {
                    const double *xj = X + N * j;
                    double *rk = Rr + N * k;
                    for (int i = 0; i < nobs; i++) rk[i] -= xj[i] * b;
                }
            }
        }
    }

    int iter;
    for (iter = 1; iter <= maxit; iter++) {
        R_CheckUserInterrupt();
        double max_change = 0.0;
        for (int j = 0; j < pp; j++) {
            if (cov) {
                for (int k = 0; k < q; k++)
                    g[k] = C[j + P * k] / n + d[j] * B[j + P * k];
            } else {
                size_t N = (size_t) nobs;
                const double *xj = X + N * j;
                for (int k = 0; k < q; k++) {
                    const double *rk = Rr + N * k;
                    double s = 0.0;
                    for (int i = 0; i < nobs; i++) s += xj[i] * rk[i];
                    g[k] = s / n + d[j] * B[j + P * k];
                }
            }
            row_update(type, q, g, d[j], l1, l2, lg, V, ev, work, bn);
            int changed = 0;
            for (int k = 0; k < q; k++) {
                double old = B[j + P * k];
                delta[k] = bn[k] - old;
                if (delta[k] != 0.0) changed = 1;
                double ac = fabs(delta[k]);
                if (ac > max_change || ISNAN(ac)) max_change = ac;
                B[j + P * k] = bn[k];
            }
            if (!changed) continue;
            if (cov) {
                const double *gj = G + P * j;
                for (int k = 0; k < q; k++) {
                    double dk = delta[k];
                    if (dk == 0.0) continue;
                    double *ck = C + P * k;
                    for (int a = 0; a < pp; a++) ck[a] -= gj[a] * dk;
                }
            } else {
                size_t N = (size_t) nobs;
                const double *xj = X + N * j;
                for (int k = 0; k < q; k++) {
                    double dk = delta[k];
                    if (dk == 0.0) continue;
                    double *rk = Rr + N * k;
                    for (int i = 0; i < nobs; i++) rk[i] -= xj[i] * dk;
                }
            }
        }
        if (max_change < tol) break;
    }
    if (iter > maxit) iter = maxit;

    SEXP out = PROTECT(allocVector(VECSXP, 2));
    SEXP nms = PROTECT(allocVector(STRSXP, 2));
    SET_VECTOR_ELT(out, 0, sB);
    SET_VECTOR_ELT(out, 1, ScalarInteger(iter));
    SET_STRING_ELT(nms, 0, mkChar("B"));
    SET_STRING_ELT(nms, 1, mkChar("iter"));
    setAttrib(out, R_NamesSymbol, nms);
    UNPROTECT(3);
    return out;
}
