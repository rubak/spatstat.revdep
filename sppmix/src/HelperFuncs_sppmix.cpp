#include "sppmix.h"
//Written by Sakis Micheas, 2017
//helper functions only, used by cpp and R functions
//visible in the package only

// [[Rcpp::export]]
double TriangleArea(vec const& a,
                    vec const& b,
                    vec const& c)
{
return 0.5*((b(0)-a(0))*(c(1)-a(1))-
  (b(1)-a(1))*(c(0)-a(0)));
}

// [[Rcpp::export]]
int CheckTriangleCLockwise(vec const& a,
                    vec const& b,
                    vec const& c)
{
  double area=TriangleArea(a,b,c);
  if (area < 0.0) return (-1);//clockwise
  else if (area > 0.0) return (+1);//counterclockwise
  return(0);
}

// [[Rcpp::export]]
bool CheckInPoly(mat const& poly,
                vec const& xy)
{
  //check is x is in poly
  //poly is an nx2 matrix of vertices
  //counterclockwise
  int i,n=poly.n_rows,crossings = 0;
  vec x=poly.col(0),y=poly.col(1);
  for (i=0;i<n-1;i++)
  {
    double slope=(y(i+1)-y(i))/(x(i+1)-x(i));
    bool cond1=(x(i)<=xy(0))&&(xy(0)<x(i+1));
    bool cond2=(x(i+1)<=xy(0))&&(xy(0)<x(i));
    bool above=(xy(1)<slope*(xy(0)-x(i))+y(i));
    if ((cond1 || cond2)  && above ) crossings++;
  }
  bool check=(crossings % 2 != 0 );
  return (check);
}

// [[Rcpp::export]]
double MatrixNorm(mat const& M,double const& p)
{
  double val=0;
  if(p<1)
    stop("Wrong power value");
  int i,j,nrow=M.n_rows,ncol=M.n_cols;
  for(i=0;i<nrow;i++)
    for(j=0;j<ncol;j++)
      val=val+powf(fabs(M(i,j)),p);
  return (powf(val,1.0/p));
}

// [[Rcpp::export]]
double MatTrace(mat const& M)
{
  double val=0;
  int i,nrow=M.n_rows;
  for(i=0;i<nrow;i++)
    val=val+M(i,i);
  return (val);
}

// [[Rcpp::export]]
double Quad_sppmix(vec const& v,mat const& m)
{
  return m(0,0)*v(0)*v(0)+m(1,1)*v(1)*v(1)
  +m(0,1)*v(0)*v(1)+m(1,0)*v(1)*v(0);
}

// [[Rcpp::export]]
double VecNorm2(vec const& v)
{
  return v(0)*v(0)+v(1)*v(1);
}

// [[Rcpp::export]]
double VecLen2(vec const& v)
{
  return sqrt(v(0)*v(0)+v(1)*v(1));
}

// [[Rcpp::export]]
double SQ_sppmix(double const& x)
{
  return x*x;
}

// [[Rcpp::export]]
double Factorial_sppmix(int x)
{
  double num=1;
  while(x>1) num*=x--;
  return num;
}


// [[Rcpp::export]]
mat invmat2d_sppmix(mat const& A){
  mat B=zeros(2,2);
  double det=A(0,0)*A(1,1)-A(0,1)*A(1,0);
  //  if(det<=0){Rcout << "\n"<<"A is not pd " << std::endl;}
  B(0,0)=A(1,1)/det;
  B(0,1)=-A(0,1)/det;
  B(1,0)=-A(1,0)/det;
  B(1,1)=A(0,0)/det;
  return B;
}


// [[Rcpp::export]]
double densNormMixatx_sppmix(vec const& atx,List const& mix
   ,vec const& approxcomp)
{
  //  densNormMix_sppmix(c(0,0),truemix)
  int j,m=mix.size();
  List mth_comp;
  double c1,val=0,psj;
  vec muk,mu1=zeros(2);
  mat invsigk,sigmak;
//    Rcout << m<< std::endl ;
//    Rcout << atx<< std::endl ;
  for(j=0;j<m;j++)
  {
    mth_comp = mix[j];
    psj=as<double>(mth_comp["p"]);
    muk = as<vec>(mth_comp["mu"]);
    //    Rcout << j<< std::endl ;
    //    Rcout << muk<< std::endl ;
    sigmak = as<mat>(mth_comp["sigma"]);
    //    Rcout << sigmak<< std::endl ;
    mu1(0)=atx(0)-muk(0);
    mu1(1)=atx(1)-muk(1);
    //    Rcout << "passed"<< std::endl ;
    invsigk=invmat2d_sppmix(sigmak);
//    Rcout << "detsigk"<< det(sigmak)<< std::endl ;
    c1=1.0/sqrt(det(2*datum::pi*sigmak));
    if(approxcomp(j)>0)
      val=val+(psj/approxcomp(j))*c1*
       exp(-.5*Quad_sppmix(mu1,invsigk));
  }
  return val;
}


// [[Rcpp::export]]
vec densNormMix_atxy_sppmix(mat const& atxy,
  List const& mix,vec const& approxcomp)
{
//atxy is an nx2 matrix
  int i,n=atxy.n_rows;
  vec val=zeros(n);
//  List mix1=mix[0];
//  double psj=as<double>(mix1["p"]);
//  Rcout << psj << std::endl ;
  for(i=0;i<n;i++)
    val(i)=densNormMixatx_sppmix(trans(atxy.row(i)),mix,approxcomp);
  return val;
}


// [[Rcpp::export]]
mat dNormMix_sppmix(List const& mix, vec const& x,
    vec const& y,vec const& approxcomp)
{
  int xnum=x.size(),ynum=y.size();
  //  int m = mix.size();
  //  List mth_comp;
  int i,j;
  mat z=zeros(xnum,ynum);
  vec atxy=zeros(2);
  for(i=0;i<xnum;i++)
    for(j=0;j<ynum;j++)
    {
      atxy(0)=x(i);
      atxy(1)=y(j);
      z(i,j)=densNormMixatx_sppmix(atxy,mix,approxcomp);
    }
    return z;
}


// [[Rcpp::export]]
vec Permute_vec_sppmix(vec const& oldvec,
                       vec const& perm)
{
  int p=perm.size();
  vec newvec(p);
  for (int i=0; i<p; i++)
    newvec(i)=oldvec(perm(i)-1);
  return newvec;
}



// [[Rcpp::export]]
mat Permute_mat_sppmix(mat const& oldmat,
                       vec const& perm)
{
  //permutes rows of a matrix
  int p=perm.size(),pp=oldmat.n_rows,q=oldmat.n_cols;
  if(p!=pp)
  {
    Rcout << "wrong dimensions" << std::endl ;
    return 0;
  }
  mat newmat(p,q);
  for (int i=0; i<p; i++)
    newmat.row(i)=oldmat.row(perm(i)-1);
  return newmat;
}


// [[Rcpp::export]]
mat GetAllPermutations_sppmix(int const& m)
{
  double permnum=Factorial_sppmix(m);
  mat allperms(permnum,m);
  int i,j;
  vec v(m);
  for(j=0;j<m;j++)v(j)=j+1;

  for (i=0;i<permnum;i++)
  {
    std::next_permutation(v.begin(),v.end());
    allperms.row(i)=v.t();
  }
  return allperms;
}


// [[Rcpp::export]]
vec GetAPermutation_sppmix(int const& m,int const& which)
{
  double permnum=Factorial_sppmix(m);
  vec aperm(m);
  int i,j;
  vec v(m);
  for(j=0;j<m;j++)v(j)=j+1;

  for (i=0;i<permnum;i++)
  {
    std::next_permutation(v.begin(),v.end());
    if (i==which)
    {
      aperm=v;
      break;
    }
  }
  return aperm;
}


// [[Rcpp::export]]
List GetGrid_sppmix(int const& len,
            vec const& xlims,vec const& ylims)
{
  //setup grid for truncation
int i,j;
vec ticsx=zeros(len),ticsy=zeros(len);
for (i=0;i<len;i++)
{
  ticsx(i)=xlims(0)+i*(xlims(1)-xlims(0))/(len-1);
  ticsy(i)=ylims(0)+i*(ylims(1)-ylims(0))/(len-1);
}
mat areas=zeros(len,len);
for(j=0;j<len-1;j++)
  for(i=0;i<len-1;i++)
    areas(i,j)=(ticsx(i+1)-ticsx(i))*(ticsy(j+1)-ticsy(j));

return List::create(
  Named("ticsx") = ticsx,
  Named("ticsy") = ticsy,
  Named("areas") = areas);
}


// [[Rcpp::export]]
bool EqVec_sppmix(vec const& v1,vec const& v2,
                  double const& tol)
{
  int len=v1.size();
  for(int i=0;i<len;i++)
  {
    if(std::abs(v1(i)-v2(i))>tol)
      return false;
  }
  return true;
}


// [[Rcpp::export]]
double logGammaFunc_sppmix(double const& x)
{
  // use approximation based on expansion of the Stirling series
  return .5*log(2*datum::pi)+(x-.5)*log(x)-x
  +1/(12*x)+1/(360*x*x*x)+1/(1260*x*x*x*x*x);
}


// [[Rcpp::export]]
double GammaFunc_sppmix(double const& x)
{
  return exp(logGammaFunc_sppmix(x));
}


// [[Rcpp::export]]
double dDirichlet_sppmix(vec const& ps,vec const& ds)
{
  int len=ps.size();
  double val=1;
  for(int i=0;i<len;i++)
    val=val*pow(ps(i),ds(i)-1)/GammaFunc_sppmix(ds(i));
  val=val*GammaFunc_sppmix(sum(ds));
  return val;
}


// [[Rcpp::export]]
double SumVec_sppmix(vec const& v,int const& start
                       ,int const& end)
{
  double val=0;
  for(int i=start;i<=end;i++)
    val+=v(i);
  return val;
}


// [[Rcpp::export]]
vec SubstituteVec_sppmix(vec v,vec const& subv,
                   int const& start)
{
  vec newvec=v;
  int len=subv.size();
  for(int i=0;i<len;i++)
    newvec(i+start)=subv(i);
  return newvec;
}


// [[Rcpp::export]]
vec SubVec_sppmix(vec const& v,int const& start,
                  int const& end)
{
//  vec newvec(end-start+1);
//  for(int i=0;i<end-start+1;i++)
//    newvec(i)=v(start+i);
//return newvec;
  return v.subvec(start,end);
}


// [[Rcpp::export]]
double GetMixtureMaxz_sppmix(List const& genmix,
                             int const& len,
                             vec const& xlims,
                             vec const& ylims,
                             vec const& approxcomp)
{
//assumes no truncation, do grid search
  int i,j;
  //  List tics=GetGrid_sppmix(len,xlims,ylims);
  //  vec ticsx=tics[1], ticsy=tics[2];
  vec ticsx=zeros(len),ticsy=zeros(len);
  for (i=0;i<len;i++)
  {
    ticsx(i)=xlims(0)+i*(xlims(1)-xlims(0))/(len-1);
    ticsy(i)=ylims(0)+i*(ylims(1)-ylims(0))/(len-1);
  }
  vec atxy=zeros(2);
  double zval,zmax=-1;
  for(i=0;i<len;i++)
    for(j=0;j<len;j++)
    {
      atxy(0)=ticsx(i);
      atxy(1)=ticsy(j);
      zval=densNormMixatx_sppmix(atxy,genmix,approxcomp);
      if(zval>zmax)
        zmax=zval;
    }
  return zmax;
/*    List kth_comp;
  double zmax=-100000000000;
  for(j=0;j<m;j++)
  {
    kth_comp = genmix[j];
    double psj=as<double>(kth_comp["p"]);
//    vec muk = as<vec>(kth_comp["mu"]);
    mat sigmak = as<mat>(kth_comp["sigma"]);
    double const1=psj/sqrt(det(2*datum::pi*sigmak));
    if (const1>zmax)
      zmax=const1;
  }
  return zmax;*/
//  return dens.max();
}


// [[Rcpp::export]]
List MakeMixtureList_sppmix(List const& gens_list,
                             int const& burnin)
{
//takes the DAMCMC output and makes a mixture list
//based on its means, sppmix::MakeMixtureList_sppmix(gens)
//find the means first
  int i,j,L=gens_list.size();
  List gen_realiz=gens_list[0];
  int countgens=L-burnin,m=gen_realiz.size();
/*  if(burnin>=L)
  {
    Rcout << "\nm="<<m << std::endl ;
    Rcout << "\nError, burnin>=L" << std::endl ;
    return List::create();
  }*/
  vec sumps=zeros(m);
  mat summus=zeros(m,2);
  mat sumsigmas=zeros(m,4);
  for(i=burnin;i<L;i++)
  {
    sumps=sumps+GetRealiz_ps_sppmix(gens_list,i);
    summus=summus+GetRealiz_mus_sppmix(gens_list,i);
    sumsigmas=sumsigmas+GetRealiz_sigmas_sppmix(gens_list,i);
  }
  //  Rcout << countgens<< L-burnin<<std::endl ;
  List mix(m);
  for(j=0;j<m;j++)
  {
    double psj=sumps(j)/countgens;
    vec muj=summus.row(j).t()/countgens;
    mat sigj(2,2);
    sigj(0,0)=sumsigmas(j,0)/countgens;
    sigj(0,1)=sumsigmas(j,1)/countgens;
    sigj(1,0)=sumsigmas(j,2)/countgens;
    sigj(1,1)=sumsigmas(j,3)/countgens;
    mix[j]=List::create(
      Named("p") = psj,
      Named("mu") = muj,
      Named("sigma") = sigj);
  }
  return (mix);
}



// [[Rcpp::export]]
List CheckInWindow_sppmix(mat const& points,
    vec const& xlims,vec const& ylims,
    bool const& truncate,
    bool const& show)
{
  int n=points.n_rows;
  //check how many points are in the window
  int countinW=0;
  vec indicesinW=zeros(n);
  for(int iin=0;iin<n;iin++)
  {
    if (points(iin,0)>=xlims(0)
          && points(iin,0)<=xlims(1)
          && points(iin,1)>=ylims(0)
          && points(iin,1)<=ylims(1))
    {
      countinW=countinW+1;
      indicesinW(iin)=1;
    }
  }
  //check if we're going to truncate
  if(!truncate)
    countinW=n;
  mat data=zeros(countinW,2);
  if(truncate)
  {
    if(show)
      Rcout <<n-countinW<<" points are outside W=["
          <<xlims(0)<<","<<xlims(1)<<"]x["
          <<ylims(0)<<","<<ylims(1)<<"]"
          << std::endl ;
    data=points.rows(find(indicesinW==1));
  }
  else
  {
    data=points;
  }
  return List::create(
    Named("count_inW") = countinW,
    Named("data_inW") = data);
}


// [[Rcpp::export]]
List GetMax_sppmix(vec const& v)
{
  int len=v.size();
  double max_value=v.max();
  int pos=-1;
  for(int i=0;i<len;i++)
  {
    if (v(i)==max_value)
    {
      pos=i;
      break;
    }
  }
  return List::create(
    Named("max") = max_value,
    Named("pos") = pos);
}


// [[Rcpp::export]]
double dNormal1d_sppmix(double const& atx,
                        double const& mu,double const& sigsq)
{
  return exp(-.5*(atx-mu)*(atx-mu)/sigsq)/sqrt(2*3.141593*sigsq);
}


// [[Rcpp::export]]
double dNormal_sppmix(vec const& atx,vec const& mu,
                      mat const& sig)
{
  double rho=sig(0,1)/sqrt(sig(0,0)*sig(1,1));
/*  Rcout <<"\nrho="<<rho<< std::endl ;
  double val1=dNormal1d_sppmix(atx(0),mu(0),sig(0,0));
  double val2=dNormal1d_sppmix(atx(1),
    mu(1)+rho*sqrt(sig(1,1)/sig(0,0))*
    (atx(0)-mu(0)),sig(1,1)*(1-rho*rho));
  Rcout <<"\nval1="<<val1<< std::endl ;
  Rcout <<"\nval2="<<val2<< std::endl ;
  Rcout <<"\nval1*val2="<<val1*val2<< std::endl ;
  return val1*val2;*/
  return dNormal1d_sppmix(atx(0),mu(0),sig(0,0))*
    dNormal1d_sppmix(atx(1),
    mu(1)+rho*sqrt(sig(1,1)/sig(0,0))*
      (atx(0)-mu(0)),sig(1,1)*(1-rho*rho));
/*  vec v(2);
  v(0)=atx(0)-mu(0);
  v(1)=atx(1)-mu(1);
  double detsig=sig(0,0)*sig(1,1)-sig(0,1)*sig(1,0);
  mat invsig(2,2);
  invsig(0,0)=sig(1,1)/detsig;
  invsig(0,1)=-sig(0,1)/detsig;
  invsig(1,0)=-sig(1,0)/detsig;
  invsig(1,1)=sig(0,0)/detsig;
  */
 // double logdens=-log(2*datum::pi)-
//    .5*log(detsig)-.5*Quad_sppmix(v,invsig);
//  return exp(logdens);
 // if(detsig<0.001)
//    Rcout <<"\nnear singular sigma, detsig="<<detsig<< std::endl ;
//  return exp(-.5*Quad_sppmix(v,invsig))/(2*3.141593*sqrt(detsig));
//  mat q=v.t()*invmat2d_sppmix(sig)*v;
//  return exp(-.5*q(0,0))/sqrt(det(2*3.141593*sig));
//Rcout <<"\nrho="<<rho<< std::endl ;
//Rcout <<"\ndetsig="<<detsig<< std::endl ;
//double val=exp(-.5*((v(0)*v(0)/sig(0,0)
//  +v(1)*v(1)/sig(1,1))-2*rho*v(0)*v(1)/
//  sqrt(sig(0,0)*sig(1,1)))/(1-rho*rho))
//  /(2*3.141593*sqrt(sig(0,0)*sig(1,1)*(1-rho*rho)));
//return val;
}

// [[Rcpp::export]]
double MultGamma(int const& p,
                 int const& n)
{
  int i;
  double sum1=0;
  for(i=0;i<p;i++)
    sum1=sum1+logGammaFunc_sppmix((n-1)/2-i/2);
  return (exp((p*(p-1)/4)*log(datum::pi)+sum1));
}

// [[Rcpp::export]]
double dInvWishart_sppmix(
    mat const& W,double const& df,
    mat const& alpha)
{
  int i,p=W.n_rows;
  double dens,const1=
    powf(arma::det(alpha),0.5*df)*
    powf(2.0,-df*0.5*p)/MultGamma(p,0.5*df);
  dens=const1*powf(det(W),-0.5*(df+p+1))*exp(-0.5*MatTrace(alpha*invmat2d_sppmix(W)));
  return(dens);
}
