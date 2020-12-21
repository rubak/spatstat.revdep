#include <R.h>
#include <vector>
#include "Pp.h"
#include "Graph.h"
#include "Rextras.h"

#define DBG 0

extern "C" {
  SEXP SGCS_Tfun_c(SEXP Args)
{
    //start parsing the args
    Args = CDR(Args);
    Pp *pp = new Pp(CAR(Args)); // init pp
    Args = CDR(Args);
    double *rvec = REAL(CAR(Args)); // r vector
    int nrvec = length(CAR(Args));
    
    // edge distances should be computed
    
    //// setup the main graph object
    double r0=0, prepr0=0;
    int gtype = 0, i0=0;
    Graph graph(pp, gtype, r0, prepr0, i0, DBG);
    
    std::vector<double > value(nrvec);
    // set old par so we start anew
    graph.oldpar = rvec[nrvec-1]-1;
    int i,j,k,m;
    int n1, n2;
    double v;
    int ok;
    for(k=nrvec-1; k > -1; k--){
      graph.par = rvec[k];  
      graph.sg_calc();
      // connections made, lets go:
      v  = 0.0;
      ok = 0;
      for(i=0; i < pp->size(); i++){
        if(rvec[k] < pp->getEdgeDistance(&i)){
          ok++;
          if(graph.nodelist.at(i).size() > 1){
            for(j=0; j < graph.nodelist.at(i).size()-1; j++){
              for(m=j+1; m < graph.nodelist.at(i).size(); m++){
                n1 = graph.nodelist.at(i).at(j)-1;
                n2 = graph.nodelist.at(i).at(m)-1;
                if(pp->getDistance(&n1, &n2) < rvec[k]) v += 1.0;
              }
            }
          }
        }
      }
      if(ok>0) v = v/(double)ok;
      graph.oldpar = rvec[k];
      value.at(k) = v;
    }
    return vectorToSEXP(value);
  }
}
