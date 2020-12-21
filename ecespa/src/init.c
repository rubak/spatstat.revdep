#include <R_ext/RS.h>
#include <stdlib.h> 
#include <R_ext/Rdynload.h>

/* .Fortran calls */
extern void F77_NAME(syrjalala)(double *x, double *y, double *var1, double *var2, int *nd, int *nperm, double *cvm, double *ks);

static const R_FortranMethodDef FortranEntries[] = {
    {"syrjalala", (DL_FUNC) &F77_NAME(syrjalala), 8},
    {NULL, NULL, 0}
};

void R_init_ecespa(DllInfo *dll)
{
    R_registerRoutines(dll, NULL, NULL, FortranEntries, NULL);
    R_useDynamicSymbols(dll, FALSE);
}

