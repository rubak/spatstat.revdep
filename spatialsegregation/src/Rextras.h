#include <R.h>
#include <Rinternals.h>
#include <Rdefines.h>
#include <vector>
SEXP getListElement(SEXP list, const char *str);

void VectsxpToVector(SEXP nodelistR, std::vector<std::vector<int> > &result);

SEXP vectorToSEXP(std::vector<std::vector<int> > nodelist);
