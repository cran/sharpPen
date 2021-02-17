/*
  *  Part of R package sharpPen
*  Copyright (C) 2020 D. Wang 
*
  *  Unlimited use and distribution (see LICENCE).
*/
  
  #include <stddef.h>
  #include <R.h>
  #include <Rinternals.h>
  #include <R_ext/Rdynload.h>
  
  void F77_SUB(blkest46)(double *x, double *y, int *n, int *q, int *qq,
                       int *nval, double *xj, double *yj, double *coef,
                       double *xmat, double *wk, double *qraux,
                       double *sigsqe, double *th44e, double *th46e);

static const R_FortranMethodDef FortEntries[] = {
  {"blkest46", (DL_FUNC) &F77_SUB(blkest46), 15},
 {NULL, NULL, 0}
};


void R_init_sharpPen(DllInfo *dll)
{
  R_registerRoutines(dll, NULL, NULL, FortEntries, NULL);
  R_useDynamicSymbols(dll, FALSE);
}

