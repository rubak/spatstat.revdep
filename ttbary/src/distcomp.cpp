#include <Rcpp.h>
#include "distcomp.h"

extern "C" {
  #include "auctionbfnew.h"
}

using namespace Rcpp;


// Here we specify functions that compute basic distances between points (including cutoffs)



// =====================================================
// R^2 with p=2
// =====================================================

// A performance test showed that the difference between dprime2 and dprimep with p=2
// is quite impressive (p=2, 25% faster)
//
double dprime2(double x0, double y0, double x1, double y1, double penp) {
  // double penp = pow(penalty,p);
  double res;
  
  if (ISNA(x1)) {
    if (ISNA(x0)) {
      res = 0.0; 
    } else {
      res = penp; 
    }
  } else {
    if (ISNA(x0)) {
      res = penp; 
    } else {
      res = (x0-x1) * (x0-x1) + (y0-y1) * (y0-y1);
      res = std::min(res, 2*penp); 
    }
  }
  return res;
}


// compute dprime crossdist to the 2nd power
// from a single point to a point pattern and return only a vector!
NumericVector cross_dprime2(double x, double y,
                            NumericVector etax, NumericVector etay, double penp) {
  int n = etax.length();
  // double penp = pow(penalty,p);
  NumericVector res(n);
  
  // j first since NumericMatrix is col major
  for (int i = 0; i < n; i++) {
    res(i) = dprime2(x, y, etax[i], etay[i], penp);
  }
  return res;
}


// compute dprime crossdist to the p-th power
NumericMatrix cross_dprime2(NumericVector xix, NumericVector xiy,
                            NumericVector etax, NumericVector etay, double penp) {
  int n = xix.length();
  if (etax.length() != n) {
    stop("cross_dprimep called with point patterns of different cardinalities");
  }  
  // double penp = pow(penalty,p);
  NumericMatrix res(n,n);
  
  // j first since NumericMatrix is col major
  for (int j = 0; j < n; j++) {
    for (int i = 0; i < n; i++) {
      res(i,j) = dprime2(xix[i], xiy[i], etax[j], etay[j], penp);
    }
  }
  return res;
}


IntegerVector closest_dprime2(double x, double y, NumericVector ppx, NumericVector ppy, double penp) {
  int ind = -1;
  int info = 1;
  
  double dist = R_PosInf;
  double mindist = R_PosInf;
  for (int i = 0; i < ppx.length(); i++) {
    dist = dprime2(x, y, ppx[i], ppy[i], penp);
    if (dist < mindist) {
      ind = i;
      mindist = dist;
    }
  }
  
  if (ISNA(ppx(ind))) {
    info = 0;
  } else if (mindist == 2*penp || ISNA(x)) {
    info = -1;
  }
  
  IntegerVector res(2);
  res[0] = ind;
  res[1] = info;
  return res;
}



// =====================================================
// R^2 with general p (typically != 2)
// =====================================================

// A performance test showed that the difference between dprime2 and dprimep with p=2
// is quite impressive (p=2, 25% faster). Probably some performance can be gained for
// doing sqrt instead of pow(...,1) in the case p=1
//
double dprimep(double x0, double y0,
               double x1, double y1, double p, double penp) {
  // double penp = pow(penalty,p);
  double res;
  
  if (ISNA(x1)) {
    if (ISNA(x0)) {
      res = 0.0; 
    } else {
      res = penp; 
    }
  } else {
    if (ISNA(x0)) {
      res = penp; 
    } else {
      res = (x0-x1) * (x0-x1) + (y0-y1) * (y0-y1);
      res = std::min(pow(res,p/2), 2*penp); 
    }
  }
  return res;
}


// compute dprime crossdist to the p-th power
// from a single point to a point pattern and return only a vector!
NumericVector cross_dprimep(double x, double y,
                            NumericVector etax, NumericVector etay, double p, double penp) {
  int n = etax.length();
  // double penp = pow(penalty,p);
  NumericVector res(n);
  
  // j first since NumericMatrix is col major
  for (int i = 0; i < n; i++) {
    res(i) = dprimep(x, y, etax[i], etay[i], p, penp);
  }
  return res;
}


// compute dprime crossdist to the p-th power
NumericMatrix cross_dprimep(NumericVector xix, NumericVector xiy,
                            NumericVector etax, NumericVector etay, double p, double penp) {
  int n = xix.length();
  if (etax.length() != n) {
    stop("cross_dprimep called with point patterns of different cardinalities");
  }  
  // double penp = pow(penalty,p);
  NumericMatrix res(n,n);
  
  // j first since NumericMatrix is col major
  for (int j = 0; j < n; j++) {
    for (int i = 0; i < n; i++) {
      res(i,j) = dprimep(xix[i], xiy[i], etax[j], etay[j], p, penp);
    }
  }
  return res;
}


IntegerVector closest_dprimep(double x, double y, NumericVector ppx, NumericVector ppy, double p, double penp) {
  int ind = -1;
  int info = 1;
  
  double dist = R_PosInf;
  double mindist = R_PosInf;
  for (int i = 0; i < ppx.length(); i++) {
    dist = dprimep(x, y, ppx(i), ppy(i), p, penp);
    if (dist < mindist) {
      ind = i;
      mindist = dist;
    }
  }
  
  if (ISNA(ppx(ind))) {
    info = 0;
  } else if (mindist == 2*penp || ISNA(x)) {
    info = -1;
  }
  
  IntegerVector res(2);
  res(0) = ind;
  res(1) = info;
  return res;
}









