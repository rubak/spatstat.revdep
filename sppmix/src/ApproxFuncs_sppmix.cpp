#include "sppmix.h"

//Written by Sakis Micheas, 2015
//approximation functions only
//used by cpp and R functions
//visible in the package only, 7 functions


//[[Rcpp::export]]
mat ApproxAvgPostIntensity(List const& genmix,
  vec const& lamdas,int const& LL,
  int const& burnin,vec const& xlims,
vec const& ylims,mat const& approxcomp)
{
  //compute the average of the posterior surfaces
  //needs to be multiplied by the lamdas
  int countiter=0,i,L =genmix.size();
  mat AvgPostIntensity = zeros(LL,LL);
  //,PostIntensityAvg = zeros(LL,LL);
  vec xy(2);
  vec ticsx=zeros(LL),ticsy=zeros(LL);
  for (i=0;i<LL;i++)
  {
    ticsx(i)=xlims(0)+i*(xlims(1)-xlims(0))/(LL-1);
    ticsy(i)=ylims(0)+i*(ylims(1)-ylims(0))/(LL-1);
  }

  //  Rcout << muk<< std::endl ;
//List tmix=genmix[0];
 // int m=tmix.size();
  //  Rcout <<m<< std::endl ;
  //  List mixcomp(m);//list containing mixture ps,mus,sigs
  //  List mixcomp;//list containing mixture ps,mus,sigs
  double intensityatxy;
  //  countiter=0;
  for(int x1=0;x1<LL;x1++)
    for(int y1=0;y1<LL;y1++)
    {
      Rprintf("\rComputing intensity surfaces: %3.1f%% complete",100.0*countiter/(LL*LL));
      xy(0)=ticsx(x1);
      xy(1)=ticsy(y1);
      //for each realization, compute the intensity
      //surface at the point xy
      for(i=burnin;i<L;i++)
      {
        //        mixcomp=genmix[i];
        //        intensityatxy=densNormMixatx_sppmix(xy,mixcomp);
//        vec approxcomp=ApproxAllCompMass_sppmix(xlims,ylims,genmix[i],truncated);
        intensityatxy=lamdas(i)*densNormMixatx_sppmix(xy,genmix[i],approxcomp.row(i).t());
        AvgPostIntensity(x1,y1)=AvgPostIntensity(x1,y1)+intensityatxy/(L-burnin);
      }
      countiter++;
      Rcpp::checkUserInterrupt();
    }
    Rprintf("\rDone                                                      \n");
  return AvgPostIntensity;
}


// [[Rcpp::export]]
double ApproxCompMass_sppmix(vec const& xlims,
  vec const& ylims,vec const& mu,mat const& sigma)
{
//approximates a single multivariate normal mass
//The grid squares must have small area size
//say, less than 0.01. If for the given LL and
//domain this is not the case, adjust LL
//first, check if a typical grid square has small area
/*  int LL;
  if((xlims(1)-xlims(0))*(ylims(1)-ylims(0))>0.01*LL*LL)
//area is too large, increase LL
  {
    LL=floor(sqrt((xlims(1)-xlims(0))*(ylims(1)-ylims(0))/0.01));
    Rcout <<"need to increase L to "<< LL<<", for a better approximation"<< std::endl ;
  }
  else
    LL=L;
  double quad,approx=0,
    d1=1/sqrt(det(2*datum::pi*sig));
  int i,j;
  vec midtics(2);
  vec ticsx=zeros(LL),ticsy=zeros(LL);
  for (i=0;i<LL;i++)
  {
    ticsx(i)=xlims(0)+i*(xlims(1)-xlims(0))/(LL-1);
    ticsy(i)=ylims(0)+i*(ylims(1)-ylims(0))/(LL-1);
  }

  for(j=0;j<LL-1;j++)
  {
    for (i=0;i<LL-1;i++)
    {
      double area=(ticsx(i+1)-ticsx(i))*(ticsy(j+1)-ticsy(j));
      midtics(0)=(ticsx(i)+ticsx(i+1))/2;
      midtics(1)=(ticsy(j)+ticsy(j+1))/2;
      quad=as_scalar((midtics-mu).t()*siginv*(midtics-mu));
      approx=approx+area*d1* exp(-.5*quad);
    }
  }*/
//this is a Monte-Carlo approximation
/*  vec xs=Rcpp::runif(L,xlims(0),xlims(1));
  vec ys=Rcpp::runif(L,xlims(0),xlims(1));
  vec xy(2);
  double quad,approx=0,
    d1=1/sqrt(det(2*datum::pi*sig));
  int i;
  for (i=0;i<L;i++)
  {
    xy(0)=xs(i);
    xy(1)=ys(i);
    quad=as_scalar((xy-mu).t()*siginv*(xy-mu));
    approx=approx+d1*exp(-.5*quad);
  }
  approx=((xlims(1)-xlims(0))*(ylims(1)-ylims(0)))*approx/L;
*/
//this is a better(?) Monte-Carlo approximation
/*  int L=1000;
  mat xy=rnorm2_sppmix(L,mu,sigma);
  double approx=0;
  for (int i=0;i<L;i++)
  {
    if(xy(i,0)>=xlims(0) && xy(i,0)<=xlims(1) &&
       xy(i,1)>=ylims(0) && xy(i,1)<=ylims(1) )
    approx++;
  }
  return approx/L;
//*/
  return ApproxBivNormProb_sppmix(xlims,
      ylims,mu,sigma,2);
}


// [[Rcpp::export]]
vec ApproxAllCompMass_sppmix(vec const& xlims,
  vec const& ylims,List const& mix,bool const& truncated)
{
  int j,m=mix.length();
  vec approxcomp=ones(m);
  if(truncated)
    for(j=0;j<m;j++)
    {
      List mth_comp = mix[j];
      vec muk = as<vec>(mth_comp["mu"]);
      mat sigmak = as<mat>(mth_comp["sigma"]);
      approxcomp(j)=ApproxBivNormProb_sppmix(xlims,
                 ylims,muk,sigmak,2);
    }
  return approxcomp;
}

// [[Rcpp::export]]
double ApproxMHRatiomu_sppmix(
  vec const& xlims,vec const& ylims,
  vec const& curmu,vec const& propmu,
  mat const& sigma,int const& num)
{
/*  int LL=L;
  if((xlims(1)-xlims(0))*(ylims(1)-ylims(0))>0.01*LL*LL)
    //area is too large, increase LL
  {
    LL=floor(sqrt((xlims(1)-xlims(0))*(ylims(1)-ylims(0))/0.01));
//    Rcout <<"need to increase L to "<< LL<<", for a better approximation"<< std::endl ;
  }
  else
    LL=L;
  double quad,approxFmu=0,approxFpropmu=0,
    d1=1/sqrt(det(2*datum::pi*sig));
  int i,j;
  vec ticsx=zeros(LL),ticsy=zeros(LL);
  for (i=0;i<LL;i++)
  {
    ticsx(i)=xlims(0)+i*(xlims(1)-xlims(0))/(LL-1);
    ticsy(i)=ylims(0)+i*(ylims(1)-ylims(0))/(LL-1);
  }
  vec midtics=zeros(2,1);
  for(j=0;j<LL-1;j++)
  {
    for (i=0;i<LL-1;i++)
    {
      double area=(ticsx(i+1)-ticsx(i))*(ticsy(j+1)-ticsy(j));
      midtics(0)=(ticsx(i)+ticsx(i+1))/2;
      midtics(1)=(ticsy(j)+ticsy(j+1))/2;
      quad=as_scalar((midtics-curmu).t()*siginv*(midtics-curmu));
      approxFmu=approxFmu+area*d1* exp(-.5*quad);
      quad=as_scalar((midtics-propmu).t()*siginv*(midtics-propmu));
      approxFpropmu=approxFpropmu+area*d1* exp(-.5*quad);
    }
  }*/
//  vec uplim(2);
//  uplim(0)=xlims(1);
//  uplim(1)=ylims(1);
//  return pow(ApproxNormCdf2d_sppmix(uplim,
//                                curmu,sigma)/
//          ApproxNormCdf2d_sppmix(uplim,
//                         propmu,sigma),num);
return pow(ApproxCompMass_sppmix(xlims,
 ylims,propmu,sigma)/
  ApproxCompMass_sppmix(xlims,ylims,
                        curmu,sigma),1.0*num);
    //  return approxFmu/approxFpropmu;
/*return pow( ApproxCompMass_sppmix(xlims,ylims,
   curmu,sigma)/ApproxCompMass_sppmix(xlims,
   ylims,propmu,sigma),num);*/
}


// [[Rcpp::export]]
double ApproxMHRatiosig_sppmix(
    vec const& xlims,vec const& ylims,
    vec const& mu,mat const& cursigma,
    mat const& propsigma,int const& num)
{
/*  int LL=L;
  if((xlims(1)-xlims(0))*(ylims(1)-ylims(0))>0.01*LL*LL)
    //area is too large, increase LL
  {
    LL=floor(sqrt((xlims(1)-xlims(0))*(ylims(1)-ylims(0))/0.01));
//    Rcout <<"need to increase L to "<< LL<<", for a better approximation"<< std::endl ;
  }
  else
    LL=L;
  double quad,approxFsig=0,approxFpropsig=0,
    d1=1/sqrt(det(2*datum::pi*sig)),
    d1prop=1/sqrt(det(2*datum::pi*propsigma));
  mat invpropsigma=invmat2d_sppmix(propsigma);
  int i,j;
  vec ticsx=zeros(LL),ticsy=zeros(LL);
  for (i=0;i<LL;i++)
  {
    ticsx(i)=xlims(0)+i*(xlims(1)-xlims(0))/(LL-1);
    ticsy(i)=ylims(0)+i*(ylims(1)-ylims(0))/(LL-1);
  }
  vec midtics=zeros(2,1);
  for(j=0;j<LL-1;j++)
  {
    for (i=0;i<LL-1;i++)
    {
      double area=(ticsx(i+1)-ticsx(i))*(ticsy(j+1)-ticsy(j));
      midtics(0)=(ticsx(i)+ticsx(i+1))/2;
      midtics(1)=(ticsy(j)+ticsy(j+1))/2;
      quad=as_scalar((midtics-mu1).t()*siginv*(midtics-mu1));
      approxFsig=approxFsig+area*d1* exp(-.5*quad);
      quad=as_scalar((midtics-mu1).t()*invpropsigma*(midtics-mu1));
      approxFpropsig=approxFpropsig+area*d1prop* exp(-.5*quad);
    }
  }
  return approxFsig/approxFpropsig;*/
/*  vec uplim(2);
//  uplim(0)=xlims(1);
//  uplim(1)=ylims(1);
  return pow(ApproxNormCdf2d_sppmix(uplim,
                  mu,cursigma)/
               ApproxNormCdf2d_sppmix(uplim,
                    mu,propsigma),num);*/
  return pow(ApproxCompMass_sppmix(
        xlims,ylims,mu,propsigma)/
       ApproxCompMass_sppmix(xlims,ylims,
         mu,cursigma),1.0*num);
/*return pow(ApproxCompMass_sppmix(xlims,ylims,
      mu,cursigma)/ApproxCompMass_sppmix(
    xlims,ylims,mu,propsigma),num);*/
}



//[[Rcpp::export]]
mat ApproxBayesianModelAvgIntensity_sppmix(
    List const& genBDmix,vec const& lamdas,
    vec const& numcomp,vec const& distr_numcomp,
    int const& mincomp,int const& maxcomp,
    int const& LL,vec const& xlims,vec const& ylims
  ,mat const& approxcomp)
{
  //apply burnin before calling this function
  //compute the average of the posterior surfaces
  //needs to be multiplied by the lamdas
  //and weighted by the comp number relative frequency
  //mincomp, maxcomp are 1,2,3,...,maxnumcomp, integers
  int countiter=0,i,j;
  mat AvgPostIntensity = zeros(LL,LL);
  vec ticsx=zeros(LL),ticsy=zeros(LL);
  for (i=0;i<LL;i++)
  {
    ticsx(i)=xlims(0)+i*(xlims(1)-xlims(0))/(LL-1);
    ticsy(i)=ylims(0)+i*(ylims(1)-ylims(0))/(LL-1);
  }
  vec xy(2);
  double weight,intensityatxy;
  Rprintf("\nComputing the Bayesian model average surface...\n");
  for(int kval=mincomp-1;kval<maxcomp;kval++)
  {
    mat PostIntensity = zeros(LL,LL);
    uvec indi=find(numcomp==kval);
    if(sum(indi)==0)//no realizations for this k
      continue;
//    Rprintf("\nTotal operation is %3.1f%% complete\n",100.0*countiter1/(maxcomp-mincomp+1));
    int nn1=indi.size();
    Rprintf("\nNumber of components=%2d, surface count=%5d\n",kval,nn1);
    vec newlamdas=lamdas(indi);
    List newgenmix(nn1);
    mat newapproxcomp(nn1,kval);
    countiter=0;
    for(i=0;i<nn1;i++)
    {
      newgenmix[i]=genBDmix[indi(i)];
      for(j=0;j<kval;j++)
        newapproxcomp(i,j)=
          approxcomp(indi(i),j);
    }
    for(int x1=0;x1<LL;x1++)
    {
      for(int y1=0;y1<LL;y1++)
      {
        Rprintf("\rComponent %2d is %3.1f%% complete",kval,100.0*countiter/(LL*LL));
        xy(0)=ticsx(x1);
        xy(1)=ticsy(y1);
        //for each realization, compute the intensity
        //surface at the point xy
        for(i=0;i<nn1;i++)
        {
          intensityatxy=newlamdas(i)*densNormMixatx_sppmix(xy,newgenmix[i],newapproxcomp.row(i).t());
          PostIntensity(x1,y1)=PostIntensity(x1,y1)+intensityatxy/nn1;
        }
        countiter++;
      }
    }
    weight=distr_numcomp(kval);
    AvgPostIntensity=AvgPostIntensity+weight*PostIntensity;
    Rcpp::checkUserInterrupt();
  }
  Rprintf("\rDone                                                                     \n");
  return AvgPostIntensity;
}

double ApproxNormCdf2d_sppmix(vec const& uplim,
    vec const& mu,mat const& sigma)
{
  vec xlims(2),ylims(2);
  xlims(0)=0;
  xlims(1)=uplim(0);
  ylims(0)=0;
  ylims(1)=uplim(1);
  return ApproxBivNormProb_sppmix(xlims,
        ylims,mu,sigma,0);
}

