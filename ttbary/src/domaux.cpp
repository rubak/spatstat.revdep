// Dom's little auxiliary functions for working with Rcpp
// Version 0.0, 16/03/2019

#include <Rcpp.h>
#include "domaux.h"

using namespace Rcpp;

// from a post by Dirk Eddelbuettel
IntegerVector whichNA(NumericVector x) {
  IntegerVector v = seq(0, x.size()-1);
  return v[is_na(x)];
}

IntegerVector whichNA(IntegerVector x) {
  IntegerVector v = seq(0, x.size()-1);
  return v[is_na(x)];
}

IntegerVector whichnotNA(NumericVector x) {
  IntegerVector v = seq(0, x.size()-1);
  return v[!is_na(x)];
}

IntegerVector whichnotNA(IntegerVector x) {
  IntegerVector v = seq(0, x.size()-1);
  return v[!is_na(x)];
}


// Return first index of x that is equal to k
IntegerVector which(IntegerVector x, int k) {
  IntegerVector v = seq(0, x.size()-1);
  IntegerVector w = v[x == k];
  return w;
} 

// Return indices of x that are equal to a
IntegerVector which(NumericVector x, double a) {
  IntegerVector v = seq(0, x.size()-1);
  IntegerVector w = v[x == a];
  return w;
} 

// Return first index i, where x[i]=a and y[i]=b 
int whichTwice(NumericVector x, NumericVector y, double a, double b) {
  IntegerVector v = seq(0, x.size()-1);
  IntegerVector w = v[(x == a) & (y == b)];
  if (w.length() == 0) {
    return -1;
  } else {
    return w(0);
  }
} 

// Return all indices where x is true
IntegerVector which(LogicalVector x) {
  IntegerVector v = seq(0, x.size()-1);
  return v[x];
}


// for Vectors that change in length but are 
// treated as fixed length in memory with an end point
// (typically used in the compination swap(x,i,n); n--;
// to erase an element)
void swap(IntegerVector& x, int i, int j) {
  int temp = x(i);
  x(i) = x(j);
  x(j) = temp;
  
  return;
}

// for Vectors that change in length but are 
// treated as fixed length in memory with an end point
// (typically used in the compination swap(x,i,n); n--;
// to erase an element)
void swap(NumericVector& x, int i, int j) {
  double temp = x(i);
  x(i) = x(j);
  x(j) = temp;
  
  return;
}

// matrix version of the above; somehow referencing
// the column did not work
void swap(IntegerMatrix& x, int i, int j, int col) {
  int temp = x(i,col);
  x(i,col) = x(j,col);
  x(j,col) = temp;
  
  return;
}



// Generate a decreasing vector of maximal length that starts with value < first,
// ends exactly in value 'last', and we get from component i-1 to component i
// by multiplication of epsfac>1
// 
// The main purpose of this function is to obtain a suitable vector for epsilon scaling
// in the auction algorithm. In the grand scheme we should have first = blowup/10.0 (=1e8 with default values)
// Cannot be Rcpp::export'ed since we still have an R function that does the same
NumericVector prepare_epsvec(double first, double last, double epsfac) {  
  int neps = ceil(log(first/last)/log(epsfac));
  // we precompute the length, otherwise like R-version
  NumericVector epsvec(neps);
  
  epsvec(neps-1) = last; // without .0 this is integer division
  for (int i = neps-2; i > -1; i--) {
    epsvec(i) = epsvec(i+1) * epsfac;
  }
  
  return(epsvec);
}
