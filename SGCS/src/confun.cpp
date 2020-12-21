#include <R.h>
#include <vector>
#include "Pp.h"
#include "Graph.h"
#include "Components.h"
#include "Rextras.h"

#define DBG 0

// two kernels
double k_box(double *r, double *d, double *h){
  if( fabs(*r-*d) < *h) return 1.0/(2.0*(*h));
  return 0.0;
}
// cumulation
double k_cumu(double *r, double *d, double *h){
  if(*d <= *r ) return 1.0;
  return 0.0;
}

extern "C" {
  SEXP SGCS_confun_c(SEXP Args)
{
    /// parsing the args ///
    Args = CDR(Args);
    Pp *pp = new Pp(CAR(Args)); // init pp
    Args = CDR(Args);
    double *rvec = REAL(CAR(Args)); // r vector
    int nrvec = length(CAR(Args));
    
    Args = CDR(Args);
    double *Rh = REAL(CAR(Args)); // R and h vector. If h=0, use cumulative
    
    double (*kernel)(double *r, double *d, double *h);
    kernel = &k_box;
    double h=0;
    bool cumu = false;
    if(Rh[1]==0) {
      kernel = &k_cumu;
      cumu = true;
    }
    else h = Rh[1];
    
    /// preGraph ///
    Args = CDR(Args);
    SEXP preGraph = CAR(Args); 
    
    // translation weights should be computed already
    
    //// setup the main graph object
    double r0=0, prepr0=0;
    int gtype = 0, i0=0;
    Graph graph(pp, gtype, r0, prepr0, i0, DBG);
    
    
    // The component network:
    if(!Rf_isNull(preGraph)){
      graph.setNodelist(preGraph);
    } 
    else{
      graph.par = Rh[0];
      graph.oldpar = graph.par;
      graph.sg_calc();
    }
    
    /// determine the components ///
    Components components;
    components.calculate(&graph);
    // and to save time...
    components.preComputeConnected();
    // set old par so we start anew
    graph.oldpar = rvec[nrvec-1]-1;
    int i,j,k;
    std::vector<double > value(nrvec);
    double v, d, w, yla, ala, r;
    int ok;
    /// main loop
    for(k=0; k < nrvec; k++){
      yla = 0.0;
      ala = 0.0;
      ok = 0;
      r = rvec[k];
      for(i=0; i < pp->size()-1; i++){
        for(j=i+1; j < pp->size(); j++){
          w = pp->getWeight(&i,&j);
          if(w>0){
            ok++;
            d = pp -> getDistance(&i,&j);
            v = (*kernel)(&r, &d, &h);
            yla += components.connected(&i,&j) * v / w;
            ala += v / w;
          }
        }
      }
      if(cumu) ala = 0.5;
      value.at(k) = yla/ala;
    }
    return vectorToSEXP(value);
  }
}


