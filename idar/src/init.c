#include <R_ext/RS.h>
#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

/* .Fortran calls */
extern void F77_NAME(mtb)(double *x1, double *y1, double *x2, double *y2, int *l1, int *l2, int *nombresp, int *nsp, double *r, int *nr, int *ltab, double *abu);

static const R_FortranMethodDef FortranEntries[] = {
    {"mtb", (DL_FUNC) &F77_NAME(mtb), 12},
    {NULL, NULL, 0}
};

void R_init_idar(DllInfo *dll)
{
    R_registerRoutines(dll, NULL, NULL, FortranEntries, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
