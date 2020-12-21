#ifndef DOMAUX_H
#define DOMAUX_H

#include <Rcpp.h>
using namespace Rcpp;

IntegerVector whichNA(NumericVector x);
IntegerVector whichNA(IntegerVector x);
IntegerVector whichnotNA(NumericVector x);
IntegerVector whichnotNA(IntegerVector x);

IntegerVector which(IntegerVector x, int k);
IntegerVector which(NumericVector x, double a);
int whichTwice(NumericVector x, NumericVector y, double a, double b);
IntegerVector which(LogicalVector x);
  
void swap(IntegerVector& x, int i, int j);
void swap(NumericVector& x, int i, int j);
void swap(IntegerMatrix& x, int i, int j, int col);

NumericVector prepare_epsvec(double first, double last, double epsfac);
#endif
