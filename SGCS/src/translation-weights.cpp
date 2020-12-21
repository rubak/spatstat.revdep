#include <R.h>
#include <Rinternals.h>
#include <vector>
#include "Pp.h"
extern "C" {
SEXP translation_weights_c(SEXP Args) {
  Args = CDR(Args);
  Pp * pp = new Pp(CAR(Args)); // read pp
  
  int n = pp->size();
  int dim = pp->getDim();
  
  std::vector<double> boxlen(dim);
  std::vector<double> w((int) (n*(n-1) * 0.5));
  
  int i,j,k;
  
  for(i=0; i < dim; i++) boxlen.at(i) = pp->getBoundingBoxExtent(i,1)  - pp->getBoundingBoxExtent(i,0);
    
  int ind = 0;
  double wij;
  for(i=0; i < n-1; i++) {
    for(j=i+1; j < n; j++) {
      wij = 1.0;
      for(k=0; k < dim; k++) wij *= boxlen.at(k) - fabs(pp->getCoord(&i,&k) - pp->getCoord(&j,&k));
      w.at(ind++) = wij;

    }
  }
  return vectorToSEXP(w);
}
}
