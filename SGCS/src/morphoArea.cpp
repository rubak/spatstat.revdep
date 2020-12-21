#include <R.h>
#include <vector>
#include "Pp.h"
#include "Graph.h"
#include "Components.h"
#include "Rextras.h"
#include "arcs.h"

#define DBG 0

extern "C" {
  SEXP SGCS_morphoArea_c(SEXP Args){
    /// parsing the args ///
    Args = CDR(Args);
    Pp *pp = new Pp(CAR(Args)); // init pp
    Args = CDR(Args);
    double *rvec = REAL(CAR(Args)); // r vector
    int nrvec = length(CAR(Args));
    
    /// setup ///
    std::vector<double>  value(nrvec);
    std::vector<std::vector<double> > arcs;
    int i, ci;
    double ten, tst;
    
    /// The graph container ///
    double r0=0, prepr0=0;
    int gtype = 0, i0=0;
    Graph *graph = new Graph(pp, gtype, r0, prepr0, i0, DBG);
    /// main loop ///
    int k;
    double R, v;
    for(k=nrvec-1; k >-1; k-- ){
      R = rvec[k];
      graph->par = R * 2.0;   // for the arc-calculator
      graph->sg_calc();
      arcs = morphoArcs(graph);
      v = 0.0;
      if(arcs.size()>0){
        for(i=0; i< (int) arcs.size(); i++){
          ten = arcs.at(i).at(2);
          tst = arcs.at(i).at(1);
          ci  = (int)arcs.at(i).at(0);
          v += R * R * (ten-tst);
          v += R * (pp->getX(&ci)*sin(ten)- pp->getY(&ci)*cos(ten) +
                    pp->getY(&ci)*cos(tst)- pp->getX(&ci)*sin(tst));
        }
        //v*= graph->par/4.0;
      }
      else v = NA_REAL;
      graph->oldpar = graph->par;
      value.at(k) = v;
    }
    return vectorToSEXP(value);
  } // eof function


}// C


