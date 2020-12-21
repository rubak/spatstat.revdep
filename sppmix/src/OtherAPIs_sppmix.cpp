//Written by Yuchen Wang, 2015
#include "sppmix.h"
//#include <Rcpp.h>
#include <mvtnormAPI.h>

using namespace Rcpp;

// [[Rcpp::export]]
double ApproxBivNormProb_sppmix(vec const& xlims,
                                vec const& ylims,vec const& mu,
                                mat const& sigma,int type)
{
  NumericVector lls(2),uls(2);
  lls(0) = (xlims(0) - mu(0)) / sqrt(sigma(0, 0));
  lls(1) = (ylims(0) - mu(1)) / sqrt(sigma(1, 1));

  uls(0) = (xlims(1) - mu(0)) / sqrt(sigma(0, 0));
  uls(1) = (ylims(1) - mu(1)) / sqrt(sigma(1, 1));

  NumericVector muxy(2);
  muxy(0)=0;
  muxy(1)=0;

  int n = 2, nu = 0, maxpts = 2000, inform;


  //INFIN INTEGER, array of integration limits flags:
  //if INFIN(I) < 0, Ith limits are (-infinity, infinity);
  //if INFIN(I) = 0, Ith limits are (-infinity, UPPER(I)];
  //if INFIN(I) = 1, Ith limits are [LOWER(I), infinity);
  //if INFIN(I) = 2, Ith limits are [LOWER(I), UPPER(I)].
  IntegerVector infin(2);
  infin(0) = type;
  infin(1) = type;

  double abseps=1.0/1000.0, releps=1.0/1000.0, error, value;
  int rnd=0;
  double corr=sigma(0, 1) / sqrt(sigma(0, 0) * sigma(1, 1));

  /* mvtnorm_C_mvtdst is defined in mvtnorm/inst/include/mvtnormAPI.h */
  mvtnorm_C_mvtdst(&n, &nu, lls.begin(), uls.begin(),
                   infin.begin(), &corr, muxy.begin(),
                   &maxpts, &abseps, &releps, &error,
                   &value, &inform, &rnd);
                   return value;
}
