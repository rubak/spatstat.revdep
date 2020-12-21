#include <R.h>
#include <vector>
#include "Pp.h"
#include "Graph.h"
#include "Components.h"
#include "Rextras.h"
#include "arcs.h"

#define DBG 0

extern "C" {
  SEXP SGCS_getArcs_c(SEXP Args){
    /// parsing the args ///
      Args = CDR(Args);
      Pp *pp = new Pp(CAR(Args)); // init pp
      Args = CDR(Args);
      double *rvec = REAL(CAR(Args)); // ball radius
      /// setup ///
      std::vector<double> values;
      std::vector<std::vector<double> > arcs;
      int i, ci;
      double ten, tst;
      
      /// The graph container ///
        double r0=0, prepr0=0;
      int gtype = 0, i0=0;
      Graph *graph = new Graph(pp, gtype, r0, prepr0, i0, DBG);
      double R = rvec[0];
      graph->par = R * 2.0;   // for the arc-calculator
      graph->sg_calc();
      arcs = morphoArcs(graph);
      if(arcs.size()>0){
        for(i=0; i< (int) arcs.size(); i++){
          values.push_back(arcs.at(i).at(0));
          values.push_back(arcs.at(i).at(1));
          values.push_back(arcs.at(i).at(2));
        }
      }
      return vectorToSEXP(values);
  } // eof function
  
  
}// C


