#include <R.h>
#include <Rinternals.h>
#include <vector>
#include "Pp.h"
extern "C" {
SEXP pairwise_distances_c(SEXP Args) {
  Args = CDR(Args);
  Pp * pp = new Pp(CAR(Args)); // read pp
  
  Args = CDR(Args);
  int tor = INTEGER(CAR(Args))[0]; // toroidal?
  bool torb = (bool) tor;
  pp -> setToroidal(&torb);
  
  int n = pp->size();
  
  std::vector<double> D((int) (n*(n-1) * 0.5));
  
  int i,j;
  int ind = 0;
  for(i=0; i < n-1; i++) {
    for(j=i+1; j < n; j++) {
      D.at(ind++) = pp -> getDistance(&i, &j);
    }
  }
  return vectorToSEXP(D);
}
}
