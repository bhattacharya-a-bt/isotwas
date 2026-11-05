#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>
#include <Rinternals.h>

/* .C calls */
extern void blasso(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);

/* .Call calls: one block per source file, here and in CallEntries, so that
   each file's routines can be added or removed on their own */

/* src/rowcd.c */
extern SEXP isotwas_rowcd(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP,
                          SEXP, SEXP, SEXP, SEXP, SEXP);

/* src/mrce-fast.c */

static const R_CMethodDef CEntries[] = {
    {"blasso", (DL_FUNC) &blasso, 13},
    {NULL, NULL, 0}
};

static const R_CallMethodDef CallEntries[] = {
    /* src/rowcd.c */
    {"isotwas_rowcd", (DL_FUNC) &isotwas_rowcd, 13},
    /* src/mrce-fast.c */
    {NULL, NULL, 0}
};

void R_init_isotwas(DllInfo *dll)
{
    R_registerRoutines(dll, CEntries, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
