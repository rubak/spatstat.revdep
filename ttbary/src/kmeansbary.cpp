#include <Rcpp.h>
#include <Rcpp/Benchmark/Timer.h>
#include "domaux.h"
#include "multimatching.h"

using namespace Rcpp;

// [[Rcpp::export]]
List kMeansBary(NumericVector zetax, NumericVector zetay, NumericMatrix ppmatx,
            NumericMatrix ppmaty, double penalty, int add_del = 1e4, int N = 200,
            double eps = 0.001, int verbose=0) { // add_del=1e4 is because Rcpp complains about both
                                                 // INT_MAX and R_PosInf (INT_MAX should be in climits resp. limits.h)
                                                 // R code passes INT_MAX correctly
  NumericVector zetax_out = clone(zetax); // Otherwise mm.zetax points to zetax and the zetax we pass
                                          // in R gets overwritten.
  NumericVector zetay_out = clone(zetay);

  // Construct MultiMatching object, which stores also some intermediate results
  // First call to optimPerm is in the constructor
  MultiMatching mm(zetax_out, zetay_out, ppmatx, ppmaty, penalty, 2.0);

  double sumttdistp_old = mm.getCost();
  double sumttdistp_new = sumttdistp_old + 1.0;  // just to be safe

  double sumttdistp_inter_old = sumttdistp_old;
  double sumttdistp_inter_new = sumttdistp_old + 1.0;


  if (verbose >= 1) {
    Rcout << "Cost after construction: " << sumttdistp_inter_old << std::endl;
    if (verbose == 2) mm.printSome(); else if (verbose >= 3) mm.printAll();
  }

  int ndel = 0, nadd = 0; // count in every iteration number of deleted and added points
  int it = 0;

  do {

    it++;
    if (verbose >= 1) {
      Rcout << std::endl << "Iteration " <<  it << ": " << std::endl;
    }

    // ----------------------------
    mm.optimBary();
    sumttdistp_inter_new = mm.getCost();
    if (verbose >= 1) {
      Rcout << "Cost after optimBary: " << sumttdistp_inter_new << std::endl;
      if (verbose == 2) mm.printSome(); else if (verbose >= 3) mm.printAll();
    }
    if((sumttdistp_inter_new - sumttdistp_inter_old)/sumttdistp_inter_old > 1e-5) {
      Rcout.precision(17);
      Rcout << "old: " << std::fixed << sumttdistp_inter_old << std::endl;
      Rcout << "new: " << std::fixed << sumttdistp_inter_new << std::endl;
      stop("Target function has substantially increased BETWEEN STEPS...");
    }
    sumttdistp_inter_old = sumttdistp_inter_new;
    // ----------------------------

    if (it <= add_del) {
      // ----------------------------
      ndel = mm.optimDelete();
      sumttdistp_inter_new = mm.getCost();
      if (verbose >= 1) {
        Rcout << ndel << " points deleted" << std::endl;
        Rcout << "Cost after optimDelete: " << sumttdistp_inter_new << std::endl;
        if (verbose == 2) mm.printSome(); else if (verbose >= 3) mm.printAll();
      }
      if((sumttdistp_inter_new - sumttdistp_inter_old)/sumttdistp_inter_old > 1e-5) {
        Rcout.precision(17);
        Rcout << "old: " << std::fixed << sumttdistp_inter_old << std::endl;
        Rcout << "new: " << std::fixed << sumttdistp_inter_new << std::endl;
        stop("Target function has substantially increased BETWEEN STEPS...");
      }
      sumttdistp_inter_old = sumttdistp_inter_new;
      // ----------------------------

      // ----------------------------
      nadd = mm.optimAdd();
      sumttdistp_inter_new = mm.getCost();
      if (verbose >= 1) {
        Rcout << nadd << " points added" << std::endl;
        Rcout << "Cost after optimAdd: " << sumttdistp_inter_new << std::endl;
        if (verbose == 2) mm.printSome(); else if (verbose >= 3) mm.printAll();
      }
      if((sumttdistp_inter_new - sumttdistp_inter_old)/sumttdistp_inter_old > 1e-5) {
        Rcout.precision(17);
        Rcout << "old: " << std::fixed << sumttdistp_inter_old << std::endl;
        Rcout << "new: " << std::fixed << sumttdistp_inter_new << std::endl;
        stop("Target function has substantially increased BETWEEN STEPS...");
      }
      sumttdistp_inter_old = sumttdistp_inter_new;
      // ----------------------------
    }

    // ----------------------------
    mm.optimPerm();
    sumttdistp_inter_new = mm.getCost();
    if (verbose >= 1) {
      Rcout << "Cost after optimPerm: " << sumttdistp_inter_new << std::endl;
      if (verbose == 2) mm.printSome(); else if (verbose >= 3) mm.printAll();
    }
    if((sumttdistp_inter_new - sumttdistp_inter_old)/sumttdistp_inter_old > 1e-5) {
      Rcout.precision(17);
      Rcout << "old: " << std::fixed << sumttdistp_inter_old << std::endl;
      Rcout << "new: " << std::fixed << sumttdistp_inter_new << std::endl;
      stop("Target function has substantially increased BETWEEN STEPS...");
    }
    sumttdistp_inter_old = sumttdistp_inter_new;
    // ----------------------------

    sumttdistp_new = mm.getCost();
    /*if ((sumttdistp_new - sumttdistp_old)/sumttdistp_old > 1e-5) { //obsolete right now, only relevant if the _inter_ tests are gone
      Rcout.precision(17);
      Rcout << "old: " << std::fixed << sumttdistp_old << std::endl;
      Rcout << "new: " << std::fixed << sumttdistp_new << std::endl;
      stop("Target function has substantially increased. Something is probably wrong...");
    }*/

    double relimp = (sumttdistp_old - sumttdistp_new)/sumttdistp_old;
    if (verbose >= 1) {
      Rcout << "Relative improvement: " << relimp << std::endl;
    }
    if (it >= N || relimp < eps) {
      if (relimp >= eps) {
        warning("Maximum number of iteration steps reached");
      }
      break;
    }

    sumttdistp_old = sumttdistp_new;

    //Rcout << "x: " << zetax_out << std::endl;
    //Rcout << "y: " << zetay_out << std::endl;

  } while (true);

  if (verbose >= 1) {
    Rcout << std::endl;
  }

  List res = List::create(Named("cost") = sumttdistp_new, _["barycenterx"] = zetax_out, _["barycentery"] = zetay_out, _["iterations"] = it);
  return res;
}



// The following function may be helpful to create starting "barycenters" for the kmeansbary algorithm:
// It samples n elements jointly from two NumericVectors of coordinates
// (Note that there are hints that sampling uniformly on the observation window gives (sometimes) a better
// starting pattern than sampling uniformly from all data points)

// [[Rcpp::export]]
List sampleFromData(int n, NumericVector ppvecx, NumericVector ppvecy) {
  int m = ppvecx.length(); // nrow * ncol

  IntegerVector ind = sample(m,n,true)-1;
  NumericVector zetax = ppvecx[ind];
  NumericVector zetay = ppvecy[ind];

  return(List::create(Named("zetax") = zetax, _["zetay"] = zetay));
}

//same as kmeansbary, but epsvec is part of the input and allows epsilon relaxation in the matchings
// [[Rcpp::export]]
List kMeansBaryEps(NumericVector epsvec, NumericVector zetax, NumericVector zetay, NumericMatrix ppmatx,
            NumericMatrix ppmaty, double penalty, int add_del = 1e4,
            IntegerVector relaxationVars = IntegerVector::create(3,1,2,2), int N = 200,
            double eps = 0.001, int verbose=0) { // add_del=1e4 is because Rcpp complains about both
                                                 // INT_MAX and R_PosInf (INT_MAX should be in climits resp. limits.h)
                                                 // R code passes INT_MAX correctly

  NumericVector zetax_out = clone(zetax); // Otherwise mm.zetax points to zetax and the zetax we pass
                                          // in R gets overwritten.
  NumericVector zetay_out = clone(zetay);

  // Construct MultiMatching object, which stores also some intermediate results
  // First call to optimPerm is in the constructor
  MultiMatching mm(zetax_out, zetay_out, ppmatx, ppmaty, penalty, 2.0);

  double sumttdistp_old = mm.getCost();
  double sumttdistp_new = sumttdistp_old + 1.0;  // just to be safe

  double sumttdistp_inter_old = sumttdistp_old;
  double sumttdistp_inter_new = sumttdistp_old + 1.0;
  
  if (verbose >= 1) {
    Rcout << "Cost after construction: " << sumttdistp_inter_old << std::endl;
    if (verbose == 2) mm.printSome(); else if (verbose >= 3) mm.printAll();
  }

  int ndel = 0, nadd = 0; // count in every iteration number of deleted and added points
  int it = 0;
  int step = std::max(relaxationVars[0], relaxationVars[1]); // make sure that we do not have an empty epsvec in the second block of assignments
  
  int epslen = epsvec.length();
  bool loop = false;
  int extraround = 0;

  do {
    it++;
    if (verbose >= 1) {
      Rcout << std::endl << "Iteration " <<  it-extraround << ": " << std::endl;
    }

    // ----------------------------
    mm.optimBary();
    sumttdistp_inter_new = mm.getCost();
    if (verbose >= 1) {
      Rcout << "Cost after optimBary: " << sumttdistp_inter_new << std::endl;
      if (verbose == 2) mm.printSome(); else if (verbose >= 3) mm.printAll();
    }
    if((sumttdistp_inter_new - sumttdistp_inter_old)/sumttdistp_inter_old > 1e-5) {
      stop("Target function has substantially increased BETWEEN STEPS...");
    }
    sumttdistp_inter_old = sumttdistp_inter_new;
    // ----------------------------

    if (it-extraround <= add_del) {
      // ----------------------------
      ndel = mm.optimDelete();
      sumttdistp_inter_new = mm.getCost();
      if (verbose >= 1) {
        Rcout << ndel << " points deleted" << std::endl;
        Rcout << "Cost after optimDelete: " << sumttdistp_inter_new << std::endl;
        if (verbose == 2) mm.printSome(); else if (verbose >= 3) mm.printAll();
      }
      if((sumttdistp_inter_new - sumttdistp_inter_old)/sumttdistp_inter_old > 1e-5) {
        stop("Target function has substantially increased BETWEEN STEPS...");
      }
      sumttdistp_inter_old = sumttdistp_inter_new;
      // ----------------------------

      // ----------------------------
      nadd = mm.optimAdd();
      sumttdistp_inter_new = mm.getCost();
      if (verbose >= 1) {
        Rcout << nadd << " points added" << std::endl;
        Rcout << "Cost after optimAdd: " << sumttdistp_inter_new << std::endl;
        if (verbose == 2) mm.printSome(); else if (verbose >= 3) mm.printAll();
      }
      if((sumttdistp_inter_new - sumttdistp_inter_old)/sumttdistp_inter_old > 1e-5) {
        stop("Target function has substantially increased BETWEEN STEPS...");
      }
      sumttdistp_inter_old = sumttdistp_inter_new;
      // ----------------------------
    }
    // ----------------------------
    do{ //do the matching more precise until mm.cost does not get worse
      //last steps of matching are ignored in the first iterations
      if  (it <= relaxationVars[0]) {
        int end = std::min(it-1,epslen-1);
        mm.optimPerm(epsvec[seq(0,end)]);
      }
      else{
        if  (step < epslen-1){
          mm.optimPerm(epsvec[seq(relaxationVars[1],step)]);
        }
        else{
          mm.optimPerm(epsvec[seq(relaxationVars[2],epslen-1)]);
        }
      }
      
      sumttdistp_inter_new = mm.getCost();
      if (verbose >= 1) {
        Rcout << "Cost after optimPerm: " << sumttdistp_inter_new << std::endl;
        if (verbose == 2) mm.printSome(); else if (verbose >= 3) mm.printAll();
      }
      if((sumttdistp_inter_new - sumttdistp_inter_old)/sumttdistp_inter_old > 1e-5) {
        if (it <= relaxationVars[0])
        {
          it++; //pretend it is the next iteration -> more precise matching
          extraround++; //we want to subtract these extrarounds from it to get the real number of iterations
          loop = true;
        }
        else{
          if(step < epslen - 1){
            it++;
            extraround++;
            loop = true;
            step += relaxationVars[3];
          }
          else{
            stop("Target function has substantially increased BETWEEN STEPS...");
          }
        }
      }
      else{
        loop = false;
        sumttdistp_inter_old = sumttdistp_inter_new;
      }
    } while (loop);
    // ----------------------------

    sumttdistp_new = mm.getCost();
    /*if ((sumttdistp_new - sumttdistp_old)/sumttdistp_old > 1e-5) { //obsolete right now, only relevant if the _inter_ tests are gone
      if(it > relaxationVars[0] && step < epslen - 1){
        step += relaxationVars[3];
      }
      else{
        stop("Target function has substantially increased. Something is probably wrong...");
      }
    }*/ 

    double relimp = (sumttdistp_old - sumttdistp_new)/sumttdistp_old;
    if (verbose >= 1) {
      Rcout << "Relative improvement: " << relimp << std::endl;
    }
    if (it-extraround >= N || relimp < eps) { //if relimp is negative we wait for new assignments
      if  (it-extraround > relaxationVars[0] && step < epslen - 1){
        step += relaxationVars[3]; //else we get an infinite loop if the mathings do not get worse
      }
      else{
        if (relimp >= eps) {
          warning("Maximum number of iteration steps reached");
        }
        break;
      }
    }

    sumttdistp_old = sumttdistp_new;

  } while (true);

  if (verbose >= 1) {
    Rcout << std::endl;
  }

  it -= extraround;

  List res = List::create(Named("cost") = sumttdistp_new, _["barycenterx"] = zetax_out, _["barycentery"] = zetay_out, _["iterations"] = it);
  return res;
}
