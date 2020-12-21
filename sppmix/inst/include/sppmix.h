//Written by Sakis Micheas, 2015
#ifndef __SPPMIX_H__
#define __SPPMIX_H__

#include <RcppArmadillo.h>
#include <Rcpp.h>
#include <Rmath.h>


using namespace arma;
using namespace Rcpp;

//Data augmentation, with truncation capabilities
//file: DAMCMC2d_sppmix.cpp
//calls the functions: rMultinomial_sppmix,
//invmat2d_sppmix,rWishart_sppmix,ApproxMHRatiomu_sppmix,
//ApproxMHRatiosig_sppmix,ApproxCompMass_sppmix,
//rDirichlet_sppmix,densNormMixatx_sppmix
List DAMCMC2d_sppmix(mat const& data,
                     vec const& xlims,
                     vec const& ylims,
                     int const& m,int const& L,
                     bool const& truncate,
                     vec const& hyperparams,
                     bool const& useKmeans);
List DAMCMC2dRMCP_sppmix(
    mat const& points,
    vec const& xlims,
    vec const& ylims,
    int const& m,int const& L,
    bool const& truncate,
    vec const& d,vec const& mu0,
    mat const& Sigma0,
    int const& df0,
    double const& sig0,
    bool const& useKmeans,
    mat const& startmus);
//Marked IPPP conditional on location
//file: MIPPCondLoc_sppmix.cpp
List MIPPCondLoc_sppmix(mat const& points,
  vec const& marks,vec const& xlims,
  vec const& ylims,int const& L,
  bool const& truncate,vec const& hyperparams,
  vec const& uniquemarks,
  bool const& discrete_mark,double const& r);
List GetProbFieldsCondLoc_sppmix(
  mat const& points,vec const& marks,
  vec const& xlims,
  vec const& ylims,int const& LL,
  vec const& meangamma,vec const& uniquemarks,
  bool const& truncate,double const& r);
List GenMarksProbCondLoc_sppmix(
    mat const& points,int const& L,
    vec const& xlims,vec const& ylims,
    vec const& meangamma,
    vec const& uniquemarks,
    bool const& truncate,double const& r);
List GetProbCondLoc_sppmix(
    mat const& points,vec const& origmarks,
    vec const& xlims,vec const& ylims,
    vec const& meangamma,
    vec const& uniquemarks,
    bool const& truncate,double const& r);
//Birth-Death MCMC
//file: BDMCMC2d_sppmix.cpp
List BDMCMC2d_sppmix(int const& maxnumcomp,
                     mat const& data,
                     vec const& xlims,
                     vec const& ylims,
                     int const& L,
                     bool const& truncate,
                     double const& lamda,
                     double const& lamdab,
                     vec const& hyper);

//Simulation functions
//file: SimulFuncs_sppmix.cpp
//double rUnif_sppmix();
//double rUnifab_sppmix(double const& a,double const& b);
mat rnorm2_sppmix(int const& n, vec const& mu,mat const& sigma);
mat rWishart_sppmix(int const& df, mat const& A);
int rDiscrete_sppmix(int const& start,vec const& probs);
int rBinom_sppmix(int const& n,double const& p);
//double rGamma_sppmix(double const& a,double const& b);
double rExp_sppmix(double const& a);
vec rDirichlet_sppmix(vec const& d);
vec rMultinomial_sppmix(int const& n,vec const& ps);
List rNormMix_sppmix(int const& lamda,List const& mix);
vec rPerm_sppmix(int const& n);
vec rmvnorm_sppmix(vec const& mu,mat const& sigma) ;

//Approximation functions
//file: ApproxFuncs_sppmix.cpp
mat ApproxAvgPostIntensity(List const& genmix,
    vec const& lamdas,int const& LL,int const& burnin,
    vec const& xlims,vec const& ylims,mat const& approxcomp);
double ApproxCompMass_sppmix(vec const& xlims,vec const& ylims,vec const& mu,mat const& sigma);
vec ApproxAllCompMass_sppmix(vec const& xlims,vec const& ylims,List const& mix,bool const& truncated);
double ApproxMHRatiomu_sppmix(vec const& xlims,
    vec const& ylims,vec const& curmu,
    vec const& propmu,mat const& sigma,
    int const& num);
double ApproxMHRatiosig_sppmix(vec const& xlims,
    vec const& ylims,vec const& mu,
    mat const& cursigma,mat const& propsigma,
    int const& num);
mat ApproxBayesianModelAvgIntensity_sppmix(
    List const& genBDmix,vec const& lamdas,
    vec const& numcomp,vec const& distr_numcomp,
    int const& mincomp,int const& maxcomp,
    int const& LL,vec const& xlims,vec const& ylims
  ,mat const& approxcomp);
double ApproxNormCdf2d_sppmix(vec const& uplim,
  vec const& mu,mat const& sigma);

//Operations on posterior realizations
//file: OpsPostGens_sppmix.cpp
vec GetPriorVals_sppmix(mat const& pp,
    List const& allgens,
    int const& priortype,
    vec const& d,vec const& mu0,
    mat const& Sigma0,
    int const& df0,
    double const& sig0);
List GetStats_sppmix(vec const& gens,double const& alpha);
mat GetAllRealiz_ps_sppmix(List const& allgens);
List GetAllRealiz_mus_sppmix(List const& allgens);
List GetAllRealiz_sigmas_sppmix(List const& allgens);
vec GetRealiz_ps_sppmix(List const& allgens,
                        int const& realiz);
mat GetRealiz_mus_sppmix(List const& allgens,
                         int const& realiz);
mat GetRealiz_sigmas_sppmix(List const& allgens,
                            int const& realiz);
List PostGenGetBestPerm_sppmix(List const& allgens);
List GetAllMeans_sppmix(List const& allgens,int const& burnin);
vec GetCompDistr_sppmix(vec const& numcomp,int const& maxnumcomp);
List GetBDCompRealiz_sppmix(List const& genBDmix,vec const& genlamdas,vec const& numcomp,int const& comp);
mat GetAvgLabelsDiscrete2Multinomial_sppmix(mat const& genzs,int const& m);
bool Check4LabelSwitching_sppmix(vec const& chain,int const& lag);
List PostGenGetBestPermIdenConstraint_sppmix(List const& allgens);
mat PermuteZs_sppmix(mat const& allgens_zs,
                     mat const& bestperm);
mat FisherInfoMat_sppmix(mat const& data,
 vec const& map_ps,mat const& map_mus,
 List const& map_sigmas,mat const& map_zs);
List GetDensityValues_sppmix(mat const& data,List const& fit,vec const& xlims, vec const& ylims);
double ComputeBayesFactor_sppmix(mat const& densvals1,mat const& densvals2);

//Helper functions
//file: HelperFuncs_sppmix.cpp
double MatTrace(mat const& M);
int CheckTriangleCLockwise(vec const& a,
                           vec const& b,
                           vec const& c);
double TriangleArea(vec const& a,
                    vec const& b,
                    vec const& c);
bool CheckInPoly(mat const& poly,
                 vec const& x);
double MatrixNorm(mat const& M,double const& p);
double VecLen2(vec const& v);
double VecNorm2(vec const& v);
double SQ_sppmix(double const& x);
double Factorial_sppmix(int x);
mat invmat2d_sppmix(mat const& A);
double densNormMixatx_sppmix(vec const& atx,List const& mix,vec const& approxcomp);
vec densNormMix_atxy_sppmix(mat const& atxy,List const& mix,vec const& approxcomp);
mat dNormMix_sppmix(List const& mix, vec const& x,vec const& y,vec const& approxcomp);
vec Permute_vec_sppmix(vec const& oldvec,vec const& perm);
mat Permute_mat_sppmix(mat const& oldmat,vec const& perm);
mat GetAllPermutations_sppmix(int const& m);
vec GetAPermutation_sppmix(int const& m,int const& which);
List GetGrid_sppmix(int const& len,vec const& xlims,vec const& ylims);
bool EqVec_sppmix(vec const& v1,vec const& v2,double const& tol=0.000001);
double dDirichlet_sppmix(vec const& ps,vec const& ds);
double logGammaFunc_sppmix(double const& x);
double GammaFunc_sppmix(double const& x);
double SumVec_sppmix(vec const& v,int const& start,int const& end);
vec SubstituteVec_sppmix(vec v,vec const& subv,int const& start);
vec SubVec_sppmix(vec const& v,int const& start,int const& end);
double GetMixtureMaxz_sppmix(List const& genmix,int const& len,vec const& xlims,vec const& ylims,
                             vec const& approxcomp);
List MakeMixtureList_sppmix(List const& gens_list,int const& burnin);
List CheckInWindow_sppmix(mat const& points,vec const& xlims,vec const& ylims,bool const& truncate,bool const& show=true);
List GetMax_sppmix(vec const& v);
double dNormal_sppmix(vec const& atx,vec const& mu,
                      mat const& sig);
double Quad_sppmix(vec const& v,mat const& m);
double dNormal1d_sppmix(double const& atx,double const& mu,double const& sigsq);
//Functions that connect to other APIs
//file: OtherAPIs_sppmix.cpp
double ApproxBivNormProb_sppmix(vec const& xlims,
  vec const& ylims,vec const& mu,
  mat const& sigma,int type);
double MultGamma(int const& p,
                 int const& n);
double dInvWishart_sppmix(
    mat const& W,double const& df,
    mat const& alpha);

#endif
