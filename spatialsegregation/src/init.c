#include <R.h>
#include <Rinternals.h>
#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

/* .External calls */
extern SEXP fun_c(SEXP);
extern SEXP graph_c(SEXP);

static const R_ExternalMethodDef ExternalEntries[] = {
  {"fun_c",   (DL_FUNC) &fun_c,   14},
  {"graph_c", (DL_FUNC) &graph_c,  7},
  {NULL, NULL, 0}
};

void R_init_spatialsegregation(DllInfo *dll)
{
  R_registerRoutines(dll, NULL, NULL, NULL, ExternalEntries);
  R_useDynamicSymbols(dll, FALSE);
}
