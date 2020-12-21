#ifndef CENTERING_H
#define CENTERING_H

#include <Rcpp.h>
using namespace Rcpp;

// Functions for computing "centers" for clusters (on various spaces, for various p, and so on)
void optimClusterCenterEuclid2(NumericVector clustx, NumericVector clusty, double& centerx, double& centery);
  
#endif
