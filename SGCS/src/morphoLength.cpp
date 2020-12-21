#include <R.h>
#include <vector>
#include "Pp.h"
#include "Graph.h"
#include "Rextras.h"
#include "arcs.h"

#define DBG 0

extern "C" {
  SEXP SGCS_morphoLength_c(SEXP Args){
    /// parsing the args ///
    Args = CDR(Args);
    Pp *pp = new Pp(CAR(Args)); // init pp
    Args = CDR(Args);
    double *rvec = REAL(CAR(Args)); // r vector
    int nrvec = length(CAR(Args));
    
    /// setup ///
    std::vector<double>  value(nrvec);
    std::vector<std::vector<double> > arcs;
    int i;
    
    /// The graph container ///
    double r0=0, prepr0=0;
    int gtype = 0, i0=0;
    Graph *graph = new Graph(pp, gtype, r0, prepr0, i0, DBG);
    /// main loop ///
    int k;
    double R, v;
    for(k=nrvec-1; k >-1; k-- ){
      R = rvec[k];
      graph->par = 2*R;
      graph->sg_calc();
      arcs = morphoArcs(graph);
      v = 0.0;
      if(arcs.size()>0){
        for(i=0; i< (int) arcs.size(); i++){
          v += arcs.at(i).at(2)-arcs.at(i).at(1);
        }
        v*= graph->par/2.0;
      }
      else v = 0;
      graph->oldpar = graph->par;
      value.at(k) = v;
    }
    return vectorToSEXP(value);
  } // eof function


}// C


