#include <Rcpp.h>
#include "clustercenter.h"

// Here we can specify functions that compute "centers" for clusters
// Currently it is not possible for the code to use information about happy
// points (and so on) with respect to the previous center, because the caller
// will not pass this information

// [[Rcpp::export]]
void optimClusterCenterEuclid2(NumericVector clustx, NumericVector clusty, double& centerx, double& centery) {
  centerx = mean(clustx);
  centery = mean(clusty);
  return;
}
