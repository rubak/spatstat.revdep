#include <Rcpp.h>
#include <Rcpp/Benchmark/Timer.h>
#include "domaux.h"
//#include "multimatching.h"
#include "multimatchingnet.h"

using namespace Rcpp;

// [[Rcpp::export]]
List kMeansBaryNet(NumericMatrix dpath, IntegerVector zeta, IntegerMatrix ppmat,
            double penalty, int add_del = 1e4, int N = 200,
            double eps = 0.001) { // add_del=1e4 is because Rcpp complains about both
                                                 // INT_MAX and R_PosInf (INT_MAX should be in climits resp. limits.h)
                                                 // R code passes INT_MAX correctly
  IntegerVector zeta_out = clone(zeta); // Otherwise the zeta we pass in R gets overwritten.

  MultiMatchingNet mmnet(dpath, zeta_out, ppmat, penalty, 1.0);

  double sumttdistp_old = mmnet.getCost();
  double sumttdistp_new = sumttdistp_old + 1.0;  // just to be safe

  double sumttdistp_inter_old = sumttdistp_old;
  double sumttdistp_inter_new = sumttdistp_old + 1.0;

  int ndel = 0, nadd = 0; // count in every iteration number of deleted and added points
  int it = 0;

  do {

    it++;
    //mmnet.printAll();
    // ----------------------------
    mmnet.optimBary();
    sumttdistp_inter_new = mmnet.getCost();
    if(sumttdistp_inter_new - sumttdistp_inter_old > 1e-8) {
      stop("1Target function has substantially increased BETWEEN STEPS...");
    }
    sumttdistp_inter_old = sumttdistp_inter_new;
    // ----------------------------
    //mmnet.printAll();

    /*if (it <= add_del) {
      // ----------------------------
      ndel = mmnet.optimDelete();
      sumttdistp_inter_new = mmnet.getCost();
      if (verbose >= 1) {
        Rcout << ndel << " points deleted" << std::endl;
        Rcout << "Cost after optimDelete: " << sumttdistp_inter_new << std::endl;
        if (verbose == 2) mmnet.printSome(); else if (verbose >= 3) mmnet.printAll();
      }
      if(sumttdistp_inter_new - sumttdistp_inter_old > 1e-8) {
        stop("Target function has substantially increased BETWEEN STEPS...");
      }
      sumttdistp_inter_old = sumttdistp_inter_new;
      // ----------------------------

      // ----------------------------
      nadd = mmnet.optimAdd();
      sumttdistp_inter_new = mmnet.getCost();
      if (verbose >= 1) {
        Rcout << nadd << " points added" << std::endl;
        Rcout << "Cost after optimAdd: " << sumttdistp_inter_new << std::endl;
        if (verbose == 2) mmnet.printSome(); else if (verbose >= 3) mmnet.printAll();
      }
      if(sumttdistp_inter_new - sumttdistp_inter_old > 1e-8) {
        stop("Target function has substantially increased BETWEEN STEPS...");
      }
      sumttdistp_inter_old = sumttdistp_inter_new;
      // ----------------------------
    }*/

    // ----------------------------
    mmnet.optimPerm();
    sumttdistp_inter_new = mmnet.getCost();
    //mmnet.printAll();
    if(sumttdistp_inter_new - sumttdistp_inter_old > 1e-8) {
      stop("2Target function has substantially increased BETWEEN STEPS...");
    }
    sumttdistp_inter_old = sumttdistp_inter_new;
    // ----------------------------

    sumttdistp_new = mmnet.getCost();
    if (sumttdistp_new - sumttdistp_old > 1e-8) {
      stop("3Target function has substantially increased. Something is probably wrong...");
    } //Abbruch wenn Gesamtkosten verschlechtert

    if (it >= N || sumttdistp_old - sumttdistp_new < eps) {
      if (sumttdistp_old - sumttdistp_new >= eps) {
        warning("Maximum number of iteration steps reached");
      }
      break;
    }

    sumttdistp_old = sumttdistp_new; //Gesamtkosten ueberschreiben fuer naechste IT

  } while (true);

  IntegerMatrix permut = mmnet.getPerm();

  //stop("alles ok");
  List res = List::create(Named("cost") = sumttdistp_new, _["barycenter"] = zeta_out, _["perm"] = permut, _["iterations"] = it);
  return res;

  }
