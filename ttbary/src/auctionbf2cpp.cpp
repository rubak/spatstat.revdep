#include <Rcpp.h>

extern "C" {
#include "auctionbfnew.h"
}

using namespace Rcpp;



// wrapper for C function (seems to tedious to maintain separate registering routines)
// [[Rcpp::export]]
List auctionbf2cpp(IntegerMatrix d, int n, IntegerVector pers_to_obj, IntegerVector obj_to_pers,
                   NumericVector price, NumericVector profit, int neps, NumericVector epsvec) {
  int* dp = d.begin();
  int* pto = pers_to_obj.begin();
  int* otp = obj_to_pers.begin();
  double* pricep = price.begin();
  double* profitp = profit.begin();
  double* epsvecp = epsvec.begin();
  
  /*  Rcpp::Rcout << "d = " << d << std::endl;
   Rcpp::Rcout << "pers_to_obj = " << pers_to_obj << std::endl;
   Rcpp::Rcout << "obj_to_pers = " << obj_to_pers << std::endl;
   Rcpp::Rcout << "price = " << price << std::endl;
   Rcpp::Rcout << "profit = " << profit << std::endl;
   Rcpp::Rcout << "epsvec = " << epsvec << std::endl;
   Rcpp::Rcout << n << " " << neps << std::endl; */
  
  auctionbf2(dp,&n,pto,otp,pricep,profitp,&neps,epsvecp);
  
  List res = List::create(Named("pers_to_obj") = pers_to_obj ,
                          _["obj_to_pers"] = obj_to_pers);
  return res;
}
