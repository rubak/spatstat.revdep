#include <R_ext/RS.h>
#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

/* .Fortran calls */
extern void F77_NAME(dixonloop)(int *k, double *N, double *R, double *Q, double *SP, int *VN, double *VarN, double *EN);

static const R_FortranMethodDef FortranEntries[] = {
    {"dixonloop", (DL_FUNC) &F77_NAME(dixonloop), 8},
    {NULL, NULL, 0}
};

void R_init_dixon(DllInfo *dll)
{
    R_registerRoutines(dll, NULL, NULL, FortranEntries, NULL);
    R_useDynamicSymbols(dll, FALSE);
}

