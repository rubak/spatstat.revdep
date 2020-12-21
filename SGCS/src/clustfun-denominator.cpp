#include <R.h>
#include <vector>
#include "Pp.h"
#include "Graph.h"
#include "Rextras.h"

#define DBG 0

extern "C" {
  SEXP SGCS_clustfun_denominator_c(SEXP Args)
{
    /// parsing the args ///
    Args = CDR(Args);
    Pp *pp = new Pp(CAR(Args)); // init pp
    Args = CDR(Args);
    double *rvec = REAL(CAR(Args)); // r vector
    int nrvec = length(CAR(Args));
    
    // edge distances should be computed already
    
    //// setup the main graph object
    double r0=0, prepr0=0;
    int gtype = 0, i0=0;
    Graph graph(pp, gtype, r0, prepr0, i0, DBG);
    
    // set old par so we start anew
    graph.oldpar = rvec[nrvec-1]-1;
    int iter, i,j,m;
    std::vector<double > value(nrvec);
    double v, cr, r;
    int ok, nn;
    for(iter=nrvec-1; iter > -1; iter--){
      r = rvec[iter];
      graph.par = r;
      graph.sg_calc();
      cr = 0.0;
      ok = 0;
      
      for(i=0; i < pp->size(); i++){
        if(r <= pp->getEdgeDistance(&i)){
          ok++;
          nn = graph.nodelist.at(i).size();  
          cr += 0.5 * nn * (nn-1);
        }
      }
      if(ok > 0) cr = cr / ( (double) ok );
      value.at(iter) = cr;
      graph.oldpar = r;
    }
    
    return vectorToSEXP(value);
  }
}


