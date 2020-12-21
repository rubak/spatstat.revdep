#include "sppmix.h"
//Written by Sakis Micheas, 2015
//contains simulation functions only, 10 functions

// [[Rcpp::export]]
mat rnorm2_sppmix(int const& n,vec const& mu,
                  mat const& sigma) {
 mat Z = arma::randn(n, 2);
  double sig1=sqrt(sigma(0,0)),
    sig2=sqrt(sigma(1,1));
  double rho=sigma(0,1)/(sig1*sig2);
  //Rcout << " " << sigma << std::endl ;
  //gen bivariate normal, mu, sig
  mat Gens(n, 2);
  colvec Z1 = Z.col(0),Z2 = Z.col(1);
  Gens.col(0) = sig1*Z1+mu(0);
  Gens.col(1) = sig2*(rho*Z1+
    sqrt(1.0-rho*rho)*Z2)+mu(1);
  //  if(det(sigma)<=0){Rcout << "\n"<<"sigma not pd "<<det(sigma)<<"\n"<<std::endl;}
  return Gens;
}

// [[Rcpp::export]]
mat rWishart_sppmix(int const& df, mat const& A){
  mat Gens=rnorm2_sppmix(df, zeros(2),A);
  return Gens.t()*Gens;
}

// [[Rcpp::export]]
int rDiscrete_sppmix(int const& start,vec const& probs)
{
  double u=Rcpp::runif(1)[0], cdf=0;
  int gen1=start,m=probs.size();
  for(int i=0;i<m;i++)
  {
    if (u>cdf && u<=cdf+probs(i))
    {
      gen1=start+i;
      break;
    }
    cdf=cdf+probs(i);
  }
  return gen1;
}

// [[Rcpp::export]]
int rBinom_sppmix(int const& n,
                  double const& p){
  vec probs(n+1);
  probs(0)=pow(1.0-p,1.0*n);
  for(int i=1;i<n;i++)
    probs(i)=((n-i+1.0)/i)*(p/(1-p))*probs(i-1);
  probs(n)=pow(p,1.0*n);
  //  Rcout << probs<< std::endl ;
  return rDiscrete_sppmix(0,probs);
}
/*
// [[Rcpp::export]]
double rGamma_sppmix(double const& a,
                     double const& b){
  double u=Rcpp::runif(1)[0];
  double y=- a * b * log(Rcpp::runif(1)[0]);
  double u0=pow(y,a - 1)*
    exp(-y*(a - 1)/(a * b))
    *exp(a - 1)/pow(a*b,a-1);
  while(u >= u0) {
    y=- a * b * log(Rcpp::runif(1)[0]);
    u0=pow(y,a - 1)*
      exp(-y*(a - 1)/(a * b))
      *exp(a - 1)/pow(a*b,a-1);
    u=Rcpp::runif(1)[0];
  }
  return y;
                     }*/

// [[Rcpp::export]]
double rExp_sppmix(double const& a)
{
  return -log(Rcpp::runif(1)[0])/a;
}


// [[Rcpp::export]]
vec rDirichlet_sppmix(vec const& d){
  int k = d.size();
  vec gens(k);
  if(k==1)
  {
    gens(0)=1;
    return gens;
  }
  double sumd=0;
  for(int i=0;i<k;i++)
  {
    gens(i)=//Rcpp::rgamma(1,d(i),1)[0];//shape,scale
      Rcpp::rgamma(1,d(i),1)[0];//shape,scale
//    Rcpp::rgamma(1,d(i),1/d(i))[0];//shape,scale
    //rExp_sppmix(1.0/d(i));
    sumd+=gens(i);
  }
  return gens/sumd;
}



// [[Rcpp::export]]
vec rMultinomial_sppmix(int const& n,vec const& ps){
  int j,i,k=ps.size();
  vec gen = zeros<vec>(k) ;
  //  if(sum(p)<1)
  //    Rcout << "\n"<<p<<","<<sum(p) << std::endl ;
  gen(0)=//Rcpp::rbinom(1,n,ps(0))[0];//
  //R::rbinom(n,ps(0));//
  rBinom_sppmix(n,ps(0));
  int sum1=0;
  double sump;
  for(j=1;j<k-1;j++)
  {
    sum1=sum1+gen(j-1);
    if (sum1==n)break;
    sump=0;
    for(i=j;i<k;i++)sump=sump+ps(i);
    gen(j)=//Rcpp::rbinom(1,n-sum1,ps(j)/sump)[0];
    //R::rbinom(n-sum1,ps(j)/sump);
    rBinom_sppmix(n-sum1,ps(j)/sump);
  }
  gen(k-1)=n-sum(gen);
  //      Rcout << " passed" << std::endl ;
  return gen;
}


// [[Rcpp::export]]
List rNormMix_sppmix(int const& lamda,
                       List const& mix)
{
  //mix[[k]] is the probability of a comp
  //mix[[k]] is a vec for mu, 2x1
  //mix[[k]] is a mat for sigma is 2x2
  int i,j,m  = mix.size(),comp=0;
  List mth_comp;
  //pick a component
  double u,sump,psj;
  int n=Rcpp::rpois(1,lamda)[0];
  //  Rcout << "mix" << std::endl ;
  //  Rcout << as<vec>(mth_comp["mu"]) << std::endl ;
  mat ret=zeros(n,2);
  vec comps(n);
  for(i=0;i<n;i++)
  {
    u=Rcpp::runif(1)[0];
    sump=0;
    for(j=0;j<m;j++)
    {
      mth_comp = mix[j];
      psj=as<double>(mth_comp["p"]);
      //    Rcout << "psj="<<psj<< std::endl ;
      if (sump<u && u<=sump+psj)
      {
        comp=j;
        break;
      }
      sump=sump+psj;
    }
    comps(i)=comp;
    //   Rcout << "comp="<<comp<< std::endl ;
    mth_comp = mix[comp];
    vec muk = as<vec>(mth_comp["mu"]);
    //  Rcout << muk<< std::endl ;
    mat sigmak = as<mat>(mth_comp["sigma"]);
    //  Rcout << sigmak<< std::endl ;
    ret.row(i)=rnorm2_sppmix(1,muk,sigmak);
  }
  return List::create(
    Named("data") = ret,
    Named("comp") = comps);
}

inline int randWrapper(const int n) { return (int)floor(1.0*unif_rand()*n); }

// [[Rcpp::export]]
vec rPerm_sppmix(int const& n)
{
  vec perm(n);
  for (int i=0; i<n; i++) perm(i)=i+1;
  std::random_shuffle(perm.begin(),perm.end(),randWrapper);
  return perm;
}


// [[Rcpp::export]]
vec rmvnorm_sppmix(vec const& mu,
                  mat const& sigma)
{
  int p=mu.size();
  vec Z=arma::randn(p);
  return mu+trans(chol(sigma))*Z;
}
