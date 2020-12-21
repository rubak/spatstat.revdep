#define UPDATE_ALLXYAVAIL
#define BLOWUP 1e9  // the blowup factor for the distances, should be <= MAX_INT
                    // the latter is on all but the most exotic systems 2,147,483,647
                    // (standard requires only 32,767)

#include <Rcpp.h>
#include "multimatching.h"
#include "domaux.h"
#include "distcomp.h"
#include "clustercenter.h"

extern "C" {
  #include "auctionbfnew.h"
}

using namespace Rcpp;


// constructor
MultiMatching::MultiMatching(NumericVector zetax, NumericVector zetay, NumericMatrix ppmatx, NumericMatrix ppmaty,
                             double penalty, double p) :
  n(zetax.length()), k(ppmatx.ncol()),
  zetax(zetax), zetay(zetay), isvirtual(is_na(zetax)), nvirtual(sum(isvirtual)),
  ppmatx(ppmatx), ppmaty(ppmaty),
  allxvals(na_omit(ppmatx)), allyvals(na_omit(ppmaty)), nallxyvals(allxvals.length()),
  perm(n,k), perminfo(n,k), pricemat(n,k),
  profitvec(n), obj_to_pers(n),
  khappy(n), happyclusterx_kn(k,n), happyclustery_kn(k,n),
  navail(0), allxavail(n*k), allyavail(n*k),
  sumttdistp(0.0), p(p), penp(pow(penalty,p)) {

  if (p != 2.0) {
    stop("Code currently only works for p=2");
  }
  //
  // perm, perminfo are properly initialized by optimPerm
  // khappy, happyclusterx/y_kn, navail, allxavail, allyavail are properly initialized
  // when first needed.
  // That is for khappy, happyclusterx/y_kn in the beginning of optimBary (then updated
  // at the beginning of optimDelete) and for navail, allxavail, allyavail only in the
  // beginning of optimAdd

  optimPerm(rep(1e8,1));
  // optimPerm(prepare_epsvec(1e8, 1.0/(n+1), 10.0));
  // optimPerm(rep(1e7,1));
  //
  // Rcout << "MultiMatching object constructed with optimal matching... Initial cost: " << sumttdistp << std::endl;
}


void MultiMatching::printAll() {
  Rcout << "=============================" << std::endl;
  Rcout << n << " points; " << k << "+1 patterns." << std::endl;
  Rcout << zetax << std::endl;
  Rcout << zetay << std::endl;
  Rcout << nvirtual << " virtual points in zeta." << std::endl;
  Rcout << isvirtual << std::endl << std::endl;
  Rcout << perm << std::endl;
  Rcout << perminfo << std::endl;
  Rcout << "happy: " << khappy << std::endl << std::endl;
  Rcout << happyclusterx_kn << std::endl;
  Rcout << happyclustery_kn << std::endl;
  Rcout << "available points: " << navail << std::endl;
  Rcout << allxavail << std::endl;
  Rcout << allyavail << std::endl << std::endl;
  Rcout << "total cost: " << sumttdistp << std::endl;
  Rcout << "=============================" << std::endl << std::endl;
  return;
}


void MultiMatching::printSome() {
  Rcout << std::endl << zetax << std::endl;
  Rcout << zetay << std::endl;
  Rcout << isvirtual << std::endl << std::endl;
  Rcout << perm << std::endl;
  Rcout << perminfo << std::endl;
  Rcout << "total cost: " << sumttdistp << std::endl << std::endl;
  return;
}


double MultiMatching::cost() {

  double res = 0.0;

  for (int j=0; j<k; j++) {
    for (int i=0; i<n; i++) {
      res += dprime2(zetax[i], zetay[i], ppmatx(perm(i,j),j), ppmaty(perm(i,j),j), penp);
    }
  }

  return(res);
}


double MultiMatching::getCost() {
  return(sumttdistp);
}


// optimize matchings with all data point pattern and update cluster information
void MultiMatching::optimPerm() {
  NumericVector epsvec = prepare_epsvec(1e8, 1.0/(n+1), 10.0);  // 1.0 needed!
  NumericVector dist(k);

  for (int j = 0; j < k; j++) { // do columnwise ttdistp-optimal matching
    // Rcout << j << ", " << std::endl;
    dist(j) = doSingleMatch(j, epsvec);
  }

  sumttdistp = sum(dist);
  return;
}


// same as optimPerm, but with the possibility to pass individual epsvecs for the
// epsilon-scaling in the auction algorithm
void MultiMatching::optimPerm(NumericVector epsvec) {

  NumericVector dist(k);

  for (int j = 0; j < k; j++) { // do columnwise ttdistp-optimal matching
    dist(j) = doSingleMatch(j, epsvec);
  }

  sumttdistp = sum(dist);
  return;
}


void MultiMatching::optimBary() {

  // updateHappyClusterInfo fills khappy, whichhappy_kn
  // needs up-to-date perminfo (and perm anyway)
  updateHappyClusterInfo();

  for (int i = 0; i < n; i++) {
    if (!isvirtual(i)) {
      if (khappy(i) > 0) {  // otherwise point is deleted in optimDelete/checkDeletePoint
                            // anyway (no need to repeat the cleanup code here)
                            // if (2*khappy(i) <= k) might also be an option ...
                            // ("small" suboptimality, saves us from updating updating Happy again)
        NumericVector clustx = happyclusterx_kn(_,i);
        NumericVector clusty = happyclustery_kn(_,i);
        optimClusterCenterEuclid2(clustx[seq(0, khappy(i)-1)], clusty[seq(0, khappy(i)-1)], zetax(i), zetay(i));

        // update perminfo for subsequent updateHappyClusterInfo (in optimDelete)
        for (int j = 0; j < k; j++) {
          if (perminfo(i,j) != 0) {  // zero stays zero and does not contribute to anything happy
            if (dprime2(zetax(i), zetay(i), ppmatx(perm(i,j),j),
                        ppmaty(perm(i,j),j), penp) == 2*penp) {
                // is the new center beyond maximal distance?
              perminfo(i,j) = -1;
            } else {
              perminfo(i,j) = 1;
            }
          }
        }  // for j

      }  // if (khappy(i) > 0)
    }  // if (!isvirtual(i)) {
  }  // for i

  sumttdistp = cost();
  return;
}


int MultiMatching::optimDelete() {

  int ndel = 0; // counts number of deletions

  // update of perminfo happens at the end of optimBary
  updateHappyClusterInfo();
  // Here update of khappy and happyclusterx/y_kn: Unfortunately both changes from
  // happy to miserable as well as "accidental" changes from miserable to happy are possible,
  // albeit maybe not too frequent... Anyway, if we don't update, at the very least the
  // 2*khappy(i) <= k check is not reliable anymore, and as far as I can see the cost could
  // even increase; the remaining checks in checkDeletePoint only seem to get somewhat
  // conservative, so maybe not updating can be an option if we leave out 2*khappy(i) <= k??
  // OR, MORE EFFICIENTLY, we don't update barycenter if before already 2*khappy(i) <= k
  // (then we have the behaviour that was intended before, BUT this time we reach the
  // correct formal behaviour that from after optimBary to after optimDelete total cost
  // does not decrease)

  for (int i=0; i<n; i++) {
    if (!isvirtual(i)) {

      bool del = checkDeletePoint(i);

      if (del) {
        ndel++;
        // Rcout << "Point " << i << " deleted" << std::endl;
      }
    }
  }

  sumttdistp = cost();
  return ndel;
}


int MultiMatching::optimAdd() {

  int nadd = 0; // counts number of additions

  if (nvirtual > 0) {
    updateAllAvailable();

    IntegerVector whichvirtual = which(isvirtual);  // indices of virtual points
    int nsample = std::min(nvirtual, navail);  // number of points to sample
    IntegerVector sam = sample(whichvirtual, nsample, false);
      // random order even if we sample everything. Seems safer, but not clear if necessary
      // Note: this determines only in which row (of perm and so on) we fill the point

    for (int isam = 0; isam < nsample; isam++) {
      int i = sam(isam);
      double propx = 0.0;
      double propy = 0.0;

      if (navail == 0) break;  // if UPDATE_ALLXYAVAIL is definied we typically remove more than just
                               // the sampled points from allx/yavail; also (other) points that
                               // switch cluster and go from miserable to happy

      int availind = sampleOneFromAvailable(propx, propy);
      // the point is written to propx, propy
      navail--;
      swap(allxavail, availind, navail);
      swap(allyavail, availind, navail);
      // deleted sampled point from allxavail/allyvail
      // (DO NOT move this to after checkAddPoint: if we use the code for "UPDATE_ALLXYAVAIL defined",
      // the sampled point would otherwise be deleted by checkAddPoint and (worse)
      // availind may not be the correct index any more due to this and/or other deletions

      bool add = checkAddPoint(i, propx, propy);

      if (add) {
        nadd++;
        // Rcout << "Point " << i << " added" << std::endl;
      }
    }
  }

  sumttdistp = cost();
  return nadd;
}


int MultiMatching::bringBackPts(int ndel){
    if (nvirtual-ndel <=1){ //less than two points are deleted in this iteration
        return(-1);
      }
      else{ //(first) half of the deleted points are brought back
        int rev = floor(0.5*(nvirtual-ndel));
        int counter=0;
        for (int i=0; i<n; i++){
            if (isvirtual(i) && counter < rev){
                isvirtual(i) = FALSE;
                counter++;
                nvirtual--;
            }
        }
        return(nvirtual);
      }
}

int MultiMatching::getn(){
    return n;
}




// ========================
// private member functions
// ========================

// update khappy, whichhappy_kn
void MultiMatching::updateHappyClusterInfo() {

  khappy.fill(0);

  for (int j = 0; j < k; j++) {
    for (int i = 0; i < n; i++) {
      if (perminfo(i,j) == 1) {  // happy point
        happyclusterx_kn(khappy(i),i) = ppmatx(perm(i,j),j);
        happyclustery_kn(khappy(i),i) = ppmaty(perm(i,j),j);
        khappy(i)++;
      }
    }
  }

  return;
}


void MultiMatching::updateAllAvailable() {

  navail = 0;

  for (int j = 0; j < k; j++) {
    for (int i = 0; i < n; i++) {
      if (perminfo(i,j) == -1) {  // miserable point
        allxavail(navail) = ppmatx(perm(i,j),j);
        allyavail(navail) = ppmaty(perm(i,j),j);
        navail++;
      }
    }
  }

  return;
}


// dedicated version of ttdist_euclidp for greater flexibility
// j = index of point pattern to be matched
// epsvec = epsilon-scaling schedule
// output is directly written to perm, perminfo and the scratch
// variable profitvec
double MultiMatching::doSingleMatch(int j, NumericVector epsvec) {
  // Special case 1: two empty point patterns
  // Cannot happen so we just catch it for now and can actually skip this later on
  if (n == 0) {
    stop("optimSingleMatchEuclidp called with zero point objects");
  }

  // dfix already taken to the p
  NumericMatrix dfix = cross_dprime2(zetax, zetay, ppmatx(_,j), ppmaty(_,j), penp);
  double dfixmax = max(dfix);

  // Special case 2: all distances 0 (e.g. all NA entries)
  if (dfixmax == 0.0) {
    perm(_,j) = seq(0,n-1);
    perminfo(_,j) = rep(0,n-1);
    for (int i = 0; i < n; i++) {
      // price(_,j) = rep(0,n) not worth it; any matching is fine; even more for profit
      // which is just a scratch vector
      if (ISNA(ppmatx(perm(i,j),j))) {
        perminfo(i,j) = 0;
      } else if (ISNA(zetax[i])) {
        perminfo(i,j) = -1;
      }
    }
    return 0.0;
  }

  IntegerMatrix dwork(n,n);
  // dwork will be a working copy of dfix that is rescaled and rounded to int.

  // j first since NumericMatrix is col major:
  for (int j = 0; j < n; j++) {
    for (int i = 0; i < n; i++) {
      dwork(i,j) = round(dfix(i,j)/dfixmax * BLOWUP);
    }
  }

  // NumericVector epsvec = prepare_epsvec(dwork,n);
  int neps = epsvec.length();
  if (neps < 1) {
    stop("epsvec has length < 1");
  }

  dwork = max(dwork)-dwork;
    // auctionbf uses a "desire matrix"

  // reset matching (auctionbf2 terminates once everything is matched):
  perm(_,j) = rep(-1,n);
  IntegerVector obj_to_pers(n,-1);

  int* dworkp = dwork.begin();
  int* pto = perm.begin() + j * n;  // points to start of column j
  int* otp = obj_to_pers.begin();
  double* pricep = pricemat.begin() + j * n;  // points to start of column j
  double* profitp = profitvec.begin();
  double* epsvecp = epsvec.begin();

  // https://herbsutter.com/2008/04/07/cringe-not-vectors-are-guaranteed-to-be-contiguous/
  // says newer C++ standards guarantee contiguity for std::vector;
  // assuming Rcpp Vectors are modelled to some extent alike

  /* The following shows that permcol = perm(_,j) makes a copy (the memory
     addresses of pto and temp are not the same).
     Question is, why is the code that passes permcol = perm(_,j) not at all
     slower (but maybe even a bit faster)?
     Rcout << pto << std::endl;
     IntegerVector permcol = perm(_,j);
     int* temp = permcol.begin();
     Rcout << temp << std::endl;  */

  /* IntegerVector perminfocol(n);
     NumericVector price(n);
     NumericVector profit(n);
     for (int i = 0; i < n; i++) {
       permcol(i) = *(pto + i);
       price(i) = *(pricep + i);
       profit(i) = *(profitp + i);
    }                              */

  /* Rcpp::Rcout << "dwork = " << dwork << std::endl;
  Rcpp::Rcout << "pers_to_obj = " << permcol << std::endl;
  Rcpp::Rcout << "obj_to_pers = " << obj_to_pers << std::endl;
  Rcpp::Rcout << "price = " << price << std::endl;
  Rcpp::Rcout << "profit = " << profit << std::endl;
  Rcpp::Rcout << "epsvec = " << epsvec << std::endl;
  Rcpp::Rcout << n << " " << neps << std::endl; */

  // stop("Designed stop");

  auctionbf2(dworkp,&n,pto,otp,pricep,profitp,&neps,epsvecp);

  // compute ttdist
  NumericVector distvec(n);
  for (int i = 0; i < n; i++) {
    distvec(i) = dfix(i,perm(i,j));
  }
  double ttdist = sum(distvec);  // p-th power already included in dfix

  // update perminfo
  perminfo(_,j) = rep(1,n);
  for (int i = 0; i < n; i++) {
    if (ISNA(ppmatx(perm(i,j),j))) {
      perminfo(i,j) = 0;
    } else if (distvec(i) == 2*penp || ISNA(zetax[i])) {
      // if the distance between points happens to be exactly 2*penp, we
      // register it as teleport rather than as move
      perminfo(i,j) = -1;
    }
  }

  return ttdist;
}


// Find the optimal cluster for (zx, zy), taking from each data pattern the -1 point closest to (zx, zy)
// (or the 0 at row i if there are no -1 points)
void MultiMatching::buildSingleCluster(int i, double zx, double zy, IntegerVector& rowinperm,
                                       NumericVector& clustx, NumericVector& clusty, IntegerVector& clustinfo) {

  for (int j = 0; j < k; j++) {
    IntegerVector permcol = perm(_,j);
    NumericVector ppmatxcol = ppmatx(_,j);
    NumericVector ppmatycol = ppmaty(_,j);
    int navail = sum(perminfo(_,j) == -1);
    if (navail == 0) {
      rowinperm(j) = i;
      clustinfo(j) = 0;  // navail == 0 can only happen if the j-th point in the i-th cluster does not exist
    } else {
      IntegerVector whichrow = which(perminfo(_,j) == -1);
      IntegerVector permcolsub = permcol[whichrow];
      IntegerVector res = closest_dprime2(zx, zy, ppmatxcol[permcolsub], ppmatycol[permcolsub], penp);
        // finds closest point to (zx,zy) in dprime-distance; returns vector of length two:
        // first component is the index of the closest point (first closest point if not unique)
        // second component is the info for this point (1 happy, -1 miserable, 0 at aleph)
      rowinperm(j) = whichrow(res(0));  // translate index back to position within the full permcol
      clustinfo(j) = res(1);
    }

    int ind = permcol(rowinperm(j));
    clustx(j) = ppmatxcol(ind);
    clusty(j) = ppmatycol(ind);
  }

  return;
}


// sample a point from allx/yavail and write it to (sampx, sampy)
int MultiMatching::sampleOneFromAvailable(double& sampx, double& sampy) {

  if (navail == 0) {
    stop("No points available. Cannot sample");
  }

  int ind = sample(navail, 1, false)(0)-1;
  sampx = allxavail[ind];
  sampy = allyavail[ind];

  return(ind);   // info needed so we know which point to remove
}


// called from optimDelete
bool MultiMatching::checkDeletePoint(int i) {

  bool del = false; // the deletion decision

  // Rcout << "Check deleting point " << i << std::endl;

  if (2*khappy(i) <= k) {
    // implies (first inequality in)
    // cost_delete - k_miser C^p = k_hap C^p <= (k_virt + k_miser) C^p <= cost_leave - k_miser C^p
    del = true;
  } else {
    // check if deleting the centerpoint would give a better result
    NumericVector clustx = happyclusterx_kn(_,i);
    NumericVector clusty = happyclustery_kn(_,i);
    NumericVector dvec = cross_dprime2(zetax(i), zetay(i),
                                        clustx[seq(0, khappy(i)-1)], clusty[seq(0, khappy(i)-1)], penp);
    double dist = sum(dvec);

    // the following is correct! Again we have subtracted k_miser * penp from both costs
    dist = dist + (k - khappy(i)) * penp;  // cost if z stays
    double dist2 = khappy(i) * penp;       // cost if z goes
    // Rcout << dist << " now <--> del " << dist2 << std::endl;
    if (dist > dist2) {
      del = true;
    }
  }

  if (del) {
    zetax(i) = NA_REAL;
    zetay(i) = NA_REAL;
    isvirtual(i) = true;
    nvirtual++;

    // update perminfo (needed for updateAllAvailable at the beginning of optimAdd)
    for (int j = 0; j < k; j++) {
      if (perminfo(i,j) == 1)  {
        perminfo(i,j) = -1;
      }
    }
  }

  return del;
}


// called from optimAdd
bool MultiMatching::checkAddPoint(int i, double& propx, double& propy) {

  bool add = false; // the addition decision

  // Rcout << "-----------------------------" << std::endl;
  // Rcout << "Check adding point " << i << std::endl;
  // Rcout << "Sampled proposal: " << propx << " " << propy << std::endl;
  double sampx = propx;  // these are saved for updating allx/yavailable (if UPDATE_ALLXYAVAIL defined)
  double sampy = propy;

  // auxiliary variables for storing information about optimal cluster for (propx, propy)
  IntegerVector rowinperm(k);   // j-th entry = row number in the perm matrix that contains the index
                                // of the j-th cluster points (i.e. point is ppmatx/y(perm(rowinperm(j),j), j))
  NumericVector clustx(k), clusty(k);  // x- and y-coordinates of cluster points
  IntegerVector clustinfo(k);          // their info with respect to (propx, propy)
                                       // 1 happy, -1 miserable, 0 at aleph

  // Find the optimal cluster for (propx, propy), taking from each data pattern the closes -1 point
  // (or the 0 at row i if there are no -1 points)
  buildSingleCluster(i, propx, propy, rowinperm, clustx, clusty, clustinfo);
  // Rcout << rowinperm << std::endl;
  // Rcout << clustx << std::endl;
  // Rcout << clusty << std::endl;
  // Rcout << clustinfo << std::endl;

  // recenter
  LogicalVector newhappy = (clustinfo == 1);
  if (sum(newhappy) > 0) {
    optimClusterCenterEuclid2(clustx[newhappy], clusty[newhappy], propx, propy);
  }
  // Rcout << "Centered proposal: " << propx << " " << propy << std::endl;
  // No update of clustinfo needed! After recentering some of the -1 might have changed
  // to 1 and some of the 1 to -1, but we use only the information about number of nonzeros for now
  // When we really add the point we have to be careful when updating perminfo

  // check if adding the sampled point would give a better result
  NumericVector dvec = cross_dprime2(propx, propy, clustx, clusty, penp);
  // Rcout << "dvec: " << dvec << std::endl;
  double dist = sum(dvec); // cost if centered (propx,propy) is added using the constructed cluster
  double dist2 = sum(clustinfo != 0) * penp; // cost with new cluster if center stays at aleph
  // Rcout << dist2 << " now <--> add " << dist << std::endl << std::endl;
  if (dist < dist2) {
    add = true;
    zetax(i) = propx;
    zetay(i) = propy;
    isvirtual(i) = false;
    nvirtual--;

    for (int j = 0; j < k; j++) {
      // perm
      int i2 = rowinperm(j);
      swap(perm,i,i2, j);

      // ========================================================
      #ifdef UPDATE_ALLXYAVAIL
        int perminfo_pre_i = perminfo(i,j);
        int perminfo_pre_i2 = perminfo(i2,j);
      #endif
      // ========================================================

      // update perminfo(i2,j)
      if (i2 != i) {  // maybe a bit faster (so we only have to compute perminfo once, below)
        if (isvirtual(i2)) {
          perminfo(i2,j) = perminfo(i,j);
          // i is virtual as well, and for all virtual points the targetinfo for a given
          // point is the same as it only depends on whether it is in aleph or mcx
        } else if (perminfo(i,j) == 0) {
          perminfo(i2,j) = 0;
        } else if (dprime2(zetax(i2), zetay(i2), ppmatx(perm(i2,j),j),
                           ppmaty(perm(i2,j),j), penp) == 2*penp) {
          // is the new point that i2 gets miserable?
          // Note that perm is already adapted, so perm(i2,j) is correct
          perminfo(i2,j) = -1;
        } else {
          perminfo(i2,j) = 1;
        }
      }

      // update perminfo(i,j). Note that clustinfo(j) may not be correct due to recentering.
      // Also i is added, so no need to check isvirtual(i)
      if (clustinfo(j) == 0) {
        perminfo(i,j) = 0;
      } else if (dvec(j) == 2*penp) {  // dvec(j) = dprime2(zetax(i), zetay(i),
                                       //                   ppmatx(perm(i,j),j), ppmaty(perm(i,j),j), penp)
        perminfo(i,j) = -1;
      } else {
        perminfo(i,j) = 1;
      }

      // Rcout << "remove available at other: " << perminfo_pre_i << " " << perminfo(i2,j) << std::endl;
      // Rcout << "remove available at new: " << perminfo_pre_i2 << " " << perminfo(i,j) << std::endl;

      // ========================================================
      #ifdef UPDATE_ALLXYAVAIL

      // THE FOLLOWING BLOCK (FIRST PART) SHOULD BE THE FIRST TO OMIT TO SAVE TIME!
      // remove point from allx/yavail if it was miserable in the z_i cluster (before z_i was added)
      // but is now happy among the i2-indices (not really a cluster, just the indices the former points
      // of the z_i cluster were moved to); rather unlikely that (many) removals happen here
      if (i != i2) {  // otherwise the search is repeated below and fails the second time
                      // if we are unlucky and this is where our sampled point comes
                      // from it might even fail the first time (and not the second)
        if (perminfo_pre_i == -1 && perminfo(i2,j) == 1) {
          searchAndDeleteAvailable(i2,j);
            // identify ppmatx/y(perm(i2,j),j) in allx/yavail and delete it
            // WE DONT CATCH THE FAILURE ANY MORE (could do "int fail = "), SEE BELOW...
        }
      }

      // remove point from allx/yavail if it was miserable among the i2-indices
      // but is now happy in the z_i cluster (after z_i is added); there will typically be quite
      // many removals here
      if (perminfo_pre_i2 == -1 && perminfo(i,j) == 1) {
        if (!(ppmatx(perm(i,j),j) == sampx && ppmaty(perm(i,j),j) == sampy)) {
          // the point (propx, propy) must always be included in the new z_i cluster
          // but was already removed from allxavail, allyavail by optimAdd after sampling;
          // So searchAndDeleteAvailable would fail, which per se is not a problem,
          // but we save some computions and are more confident that all works well by
          // checking for failure and just excluding this case
          searchAndDeleteAvailable(i,j);
            // identify ppmatx/y(perm(i,j),j) in allx/yavail and delete it
            // WE DONT CATCH THE FAILURE ANY MORE (could do "int fail = ...; if (fail) ..."):
            // 1. Nothing bad happens if we fail, even if there is no good reason to fail.
            //    The worst consequence would be that we leave a point in allx/yavail of
            //    which we could now already that adding it will not help
            // 2. Tests showed that this function fails only rarely for the following very
            //    good reason: a point was sampled in optimAdd and immediately removed
            //    from allx/yavail, but it turned out that it could not be added to zeta.
            //    We do not want to resample it, but it still has a -1 in perminfo nad
            //    is thus a very legitimate candidate for a new cluster.
            // Any other way of solving this seems to be computationally more intensive
            // than searching this point *once* in vain and/or destroys the flexibility
            // to control the updating of allxavail by defining/undefining UPDATE_ALLXYAVAIL
            // (and keeping the additional blocks more or less clean)
          }
      }

      #endif
      // ========================================================

    }  // j-loop
  }  // if (dist < dist2)

  return add;
}


// find ppmatx/y(perm(i,j),j) in allx/yavail, delete its first occurence and return false
// OR return true if no such occurence is found
bool MultiMatching::searchAndDeleteAvailable(int i, int j) {

  bool fail = true;
  double a = ppmatx(perm(i,j),j);
  double b = ppmaty(perm(i,j),j);
  int ind = whichTwice(allxavail[seq(0,navail-1)], allyavail[seq(0,navail-1)], a, b);

  // delete the entry at position ind gracefully from allxavail, allyavail
  if (ind >= 0) {  //
    navail--;
    swap(allxavail, ind, navail);
    swap(allyavail, ind, navail);
    fail = false;
  }

  return fail;
}
