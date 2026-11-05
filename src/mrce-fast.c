/*
 * Faster equivalent of blasso() in mrce-solveb.c, for the case rblasso()
 * computes the starting soft-threshold values itself (soft = NULL).
 *
 * blasso keeps soft (p x q) current: every change to one coefficient B[r,c]
 * updates all p*q entries, so a sweep costs O(p q) per changed coefficient.
 * This version keeps U = S B (p x q) instead, which one changed coefficient
 * touches in a single column (O(p)), and forms the soft value for (r, c)
 * when it is visited (O(q)):
 *
 *   soft[r,c] = M[r,c] - sum_b U[r,b] Om[c,b] + S[r,r] Om[c,c] B[r,c]
 *
 * that is, M - S B Om with the coefficient's own term added back: the
 * partial-residual score whose soft-threshold minimises the objective in
 * B[r,c]. blasso tracks the same value, since it never updates an entry's
 * soft value when that entry itself changes and rblasso() starts it from
 * M - S B0 Om + diag(S) diag(Om)' * B0.
 *
 * The objective change is formed the same way:
 *   sum over (a,b) != (r,c) of B[a,b] S[a,r] Om[b,c] d
 *     = d (sum_b U[r,b] Om[b,c] - S[r,r] B[r,c] Om[c,c])
 *
 * Visit order, thresholding, the objective bookkeeping and the stopping rule
 * are blasso's, so the two agree to rounding error. S must be symmetric (it is
 * crossprod of centred X).
 */

#include <R.h>
#include <Rinternals.h>
#include <R_ext/Utils.h>
#include <math.h>
#include <string.h>

SEXP isotwas_blasso_fast(SEXP sB0, SEXP sS, SEXP sM, SEXP sOm, SEXP sLam,
                         SEXP sN, SEXP sTol, SEXP sMaxit)
{
    const int *dim = INTEGER(getAttrib(sB0, R_DimSymbol));
    const int p = dim[0], q = dim[1];
    const size_t P = (size_t) p, Q = (size_t) q;
    const double *S = REAL(sS), *M = REAL(sM), *Om = REAL(sOm);
    const double *lam = REAL(sLam);
    const double n = asReal(sN), tol = asReal(sTol);
    const int maxit = asInteger(sMaxit);

    SEXP sB = PROTECT(duplicate(sB0));
    double *B = REAL(sB);
    double *U = (double *) R_alloc(P * Q, sizeof(double));
    memset(U, 0, sizeof(double) * P * Q);
    for (int b = 0; b < q; b++) {
        double *ub = U + P * b;
        for (int a = 0; a < p; a++) {
            double v = B[a + P * b];
            if (v == 0.0) continue;
            const double *sa = S + P * a;   /* column a of S */
            for (int i = 0; i < p; i++) ub[i] += sa[i] * v;
        }
    }

    double obj = 0.0;
    double max_diff = tol + 1.0;
    int kit = 0;
    while ((max_diff > tol) && (kit < maxit)) {
        R_CheckUserInterrupt();
        kit++;
        max_diff = 0.0;
        for (int c = 0; c < q; c++) {
            const double omcc = Om[c + Q * c];
            for (int r = 0; r < p; r++) {
                const size_t rc = r + P * c;
                const double srr = S[r + P * r];
                double us = 0.0;
                for (int b = 0; b < q; b++) us += U[r + P * b] * Om[c + Q * b];
                const double soft = M[rc] - us + srr * omcc * B[rc];

                double H = soft, AH = fabs(H), tmp = AH - lam[rc], bnew = 0.0;
                if (tmp > 0.0) {
                    if (H > 0.0) bnew = tmp;
                    else if (H < 0.0) bnew = -tmp;
                    else bnew = 0.0;
                }
                bnew = bnew / (srr * omcc);
                if (bnew != B[rc]) {
                    const double old = B[rc];
                    const double thisdiff = old - bnew;
                    double uo = 0.0;
                    for (int b = 0; b < q; b++) uo += U[r + P * b] * Om[b + Q * c];
                    double this_obj_diff = thisdiff * (uo - srr * old * omcc);
                    this_obj_diff += thisdiff * (-M[rc] - soft + srr * omcc * (old + bnew));
                    this_obj_diff += 2.0 * lam[rc] * (fabs(old) - fabs(bnew));
                    this_obj_diff = this_obj_diff / n;
                    B[rc] = bnew;
                    double *uc = U + P * c;
                    const double *sr = S + P * r;   /* column r of S */
                    for (int i = 0; i < p; i++) uc[i] -= sr[i] * thisdiff;
                    obj -= this_obj_diff;
                    if (this_obj_diff > max_diff) max_diff = this_obj_diff;
                }
            }
        }
    }

    SEXP out = PROTECT(allocVector(VECSXP, 3));
    SEXP nms = PROTECT(allocVector(STRSXP, 3));
    SET_VECTOR_ELT(out, 0, sB);
    SET_VECTOR_ELT(out, 1, ScalarReal(obj));
    SET_VECTOR_ELT(out, 2, ScalarInteger(kit));
    SET_STRING_ELT(nms, 0, mkChar("B"));
    SET_STRING_ELT(nms, 1, mkChar("obj"));
    SET_STRING_ELT(nms, 2, mkChar("iterations"));
    setAttrib(out, R_NamesSymbol, nms);
    UNPROTECT(3);
    return out;
}
