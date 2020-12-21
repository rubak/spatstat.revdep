#include <R.h>
#include <vector>
#include "Pp.h"
#include "Graph.h"
#include "Rextras.h"

#define DBG 0

extern "C" {
  SEXP SGCS_Kfun_c(SEXP Args)
{
  //start parsing the args
  Args = CDR(Args);
  Pp *pp = new Pp(CAR(Args)); // init pp
  Args = CDR(Args);
  double *rvec = REAL(CAR(Args)); // r vector
  int nrvec = length(CAR(Args));


  //// setup the main graph object
  double r0=0, prepr0=0;
  int gtype = 0, i0=0;
  Graph graph(pp, gtype, r0, prepr0, i0, DBG);
  
  std::vector<double > value(nrvec);
  // set old par so we start anew
  graph.oldpar = rvec[nrvec-1]-1;
  graph.dbg = 0;
  int i,j,k,l;
  double v;
  for(k=nrvec-1; k > -1; k--){
    graph.par = rvec[k];  
    graph.sg_calc();
    v  = 0.0;
    for(i=0; i < pp->size(); i++){
      for(j=0; j < graph.nodelist.at(i).size(); j++){
        l = graph.nodelist.at(i).at(j)-1;
        v += 1.0/pp->getWeight(&i,&l);
      }
    }
    graph.oldpar = rvec[k];
    value.at(k) = v;
    
  }
  return vectorToSEXP(value);
}
}
