#include <R.h>
#include <Rinternals.h>
#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

  /* .External calls */
extern SEXP pairwise_distances_c(SEXP);
extern SEXP SGCS_clustfun_c(SEXP);
extern SEXP SGCS_clustfun_denominator_c(SEXP);
extern SEXP SGCS_confun_c(SEXP);
extern SEXP SGCS_getArcs_c(SEXP);
extern SEXP SGCS_Kfun_c(SEXP);
extern SEXP SGCS_morphoArea_c(SEXP);
extern SEXP SGCS_morphoEuler_c(SEXP);
extern SEXP SGCS_morphoLength_c(SEXP);
extern SEXP SGCS_Rfun_c(SEXP);
extern SEXP SGCS_Tfun_c(SEXP);
extern SEXP translation_weights_c(SEXP);

static const R_ExternalMethodDef ExternalEntries[] = {
  {"pairwise_distances_c",        (DL_FUNC) &pairwise_distances_c,        2},
  {"SGCS_clustfun_c",             (DL_FUNC) &SGCS_clustfun_c,             2},
  {"SGCS_clustfun_denominator_c", (DL_FUNC) &SGCS_clustfun_denominator_c, 2},
  {"SGCS_confun_c",               (DL_FUNC) &SGCS_confun_c,               4},
  {"SGCS_getArcs_c",              (DL_FUNC) &SGCS_getArcs_c,              2},
  {"SGCS_Kfun_c",                 (DL_FUNC) &SGCS_Kfun_c,                 2},
  {"SGCS_morphoArea_c",           (DL_FUNC) &SGCS_morphoArea_c,           2},
  {"SGCS_morphoEuler_c",          (DL_FUNC) &SGCS_morphoEuler_c,          2},
  {"SGCS_morphoLength_c",         (DL_FUNC) &SGCS_morphoLength_c,         2},
  {"SGCS_Rfun_c",                 (DL_FUNC) &SGCS_Rfun_c,                 2},
  {"SGCS_Tfun_c",                 (DL_FUNC) &SGCS_Tfun_c,                 2},
  {"translation_weights_c",       (DL_FUNC) &translation_weights_c,       1},
  {NULL, NULL, 0}
};

void R_init_SGCS(DllInfo *dll)
{
  R_registerRoutines(dll, NULL, NULL, NULL, ExternalEntries);
  R_useDynamicSymbols(dll, FALSE);
}

