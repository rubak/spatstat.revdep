#ifndef MATCHING_H
#define MATCHING_H

#include <Rcpp.h>
using namespace Rcpp;


class MultiMatching {

  int n, k; // no. of points (in all patterns the same, including virtual points), no. of patterns
  NumericVector zetax, zetay;
  LogicalVector isvirtual;
  int nvirtual;  // order!
  const NumericMatrix ppmatx, ppmaty;
  const NumericVector allxvals, allyvals;
  const int nallxyvals;   // length of allxvals and allyvals (order!)

  IntegerMatrix perm, perminfo;  // perm(_,j) contains order of rows of ppmatx/y as matched to zetax
                                 // hence ppmatx(perm(i,j),j) for j = 1,...,k are (x-coords of) points of
                                 // the cluster with center zetax(i)
                                 // perminfo(i,j) = 0 if such a point is in aleph, -1 if the point is miserable
                                 // (distance >= 2^(1/p)*penalty away from zeta(i) OR zetax(i) is virtual and
                                 // point is not), +1 if point is happy (i.e. matched to zeta(i) at dist < 2^(1/p)*penalty)
  NumericMatrix pricemat; // pricemat(_,j) contains the prices that are calculated by auctionbf for the matching
                          // of zeta and the jth data point pattern (entry (i,j): price of i-th point/object of
                          // j-th data point pattern)
  NumericVector profitvec; // scratch; stores the profits that are calculated by auctionbf
  NumericVector obj_to_pers; // scratch; stores the permutations as seen from the second point pattern
  IntegerVector khappy; // the no. of positive values in the rows of perminfo, respectively
  NumericMatrix happyclusterx_kn, happyclustery_kn; // k times n matrices for faster access:
                                  // the i-th column contains from 0-th to (khappy(i)-1)st row
                                  // the x/y-coordinates of the happy points associated with zetax/y(i)
// PERFORMANCE: given how we fill happycluster and that we read only once from it per filling
//              it might be just as well (or better?) to work with the transposed one!!
  int navail;
  NumericVector allxavail, allyavail;
    // vectors of x/y-coordinates that we "currently want to sample from"
    // These vectors only have to be maintained for sampling when trying to add points in optimAdd.
    // They are only properly initialized at the beginning of optimAdd based on all miserable data points,
    // During optimAdd they might temporarily be different from the miserable points
    // (e.g. sampled points are immediately removed from allxavail, allyavail;
    // if the variable UPDATE_ALLXYAVAIL is not defined in multimatching.cpp, then allx/yavail are
    // not updated in optimAdd elsewhere (i.e. apart from removing the sampled), but miserable
    // points in perminfo are)

  double sumttdistp;  // stores the current total cost of the matching
  const double p, penp;  // penp = penalty^p, currently p=2, but penalty is never needed directly


public:
  // MultiMatching() {}; problematic because of the const members that have to be defined then
  MultiMatching(NumericVector zetax, NumericVector zetay, NumericMatrix ppmatx, NumericMatrix ppmaty,
                double penalty, double p);

  void printAll();
  void printSome();

  double cost();    // compute cost for current perm (perminfo not used!)
  double getCost(); // read cost from sumttdistp

  // The four main steps of the kmeansbary algorithm (Variant 1):
  void optimPerm();
  void optimPerm(NumericVector epsvec);  // passing vector for eps-scaling
  void optimBary(); // calls updateHappyClusterInfo() in the beginning
  int optimDelete(); // calls updateHappyClusterInfo() in the beginning, then checkDeletePoint()
  int optimAdd();    // calls updateAllAvailable() in the beginning, then checkAddPoint()
                     // (as well as sampleOneFromAvailable, buildSingleCluster,
                     // searchAndDeleteAvailable if UPDATE_ALLXYAVAIL is defined)

  // Variant 2 omits optimDelete and optimAdd
  // Does optimDelete once algo gets stuck, calls the following, then restarts
  int bringBackPts(int ndel);
  int getn();


private:
  void updateHappyClusterInfo();  // called from optimPerm; update all data members having to do with
                             // clusters: khappy, happycluster_kn, navail, allxavail, allyavail
  void updateAllAvailable();

  double doSingleMatch(int j, NumericVector epsvec); // called from optimPerm
  bool checkDeletePoint(int i); // called from optimDelete

  // The following are all called from optim Add:
  int sampleOneFromAvailable(double& sampx, double&sampy);
  bool checkAddPoint(int i, double& propx, double& propy);
  // NumericMatrix sampleFromAvailable(int npoints);   // probably no longer needed
  void buildSingleCluster(int i, double zx, double zy, IntegerVector& rowinperm,
                          NumericVector& clustx, NumericVector& clusty, IntegerVector& clustinfo);
    // build a new cluster around point (zx,zy) taking from each data pattern
    // the closest (more precisely cheapest) miserable point
  bool searchAndDeleteAvailable(int i, int j);
    // update allxavail, allyavail after adding a point
    // only needed if we define UPDATE_ALLXYAVAIL in .cpp file
};

#endif
