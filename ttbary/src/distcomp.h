#ifndef DISTCOMP_H
#define DISTCOMP_H

#include <Rcpp.h>
using namespace Rcpp;

// Here we can specify functions that compute basic distances between points (including cutoffs)

// =====================================================
// R^2 with p=2
// =====================================================

double dprime2(double x0, double y0, double x1, double y1, double penp);
NumericVector cross_dprime2(double x, double y, NumericVector etax, NumericVector etay, double penp);
NumericMatrix cross_dprime2(NumericVector xix, NumericVector xiy,
                            NumericVector etax, NumericVector etay, double penp);
IntegerVector closest_dprime2(double x, double y, NumericVector ppx, NumericVector ppy, double penp);


// =====================================================
// R^2 with general p (typically != 2)
// =====================================================

double dprimep(double x0, double y0, double x1, double y1, double p, double penp);
NumericVector cross_dprimep(double x, double y, NumericVector etax, NumericVector etay,
                            double p, double penp);
NumericMatrix cross_dprimep(NumericVector xix, NumericVector xiy,
                            NumericVector etax, NumericVector etay, double p, double penp);
IntegerVector closest_dprimep(double x, double y, NumericVector ppx, NumericVector ppy,
                             double p, double penp);


#endif
