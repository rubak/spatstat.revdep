#include "sppmix.h"
//Written by Sakis Micheas, 2017

//[[Rcpp::export]]
List DAMCMC2d_sppmix(mat const& points,
                     vec const& xlims,
                     vec const& ylims,
                     int const& m,int const& L,
                     bool const& truncate,
                     vec const& hyperparams,
                     bool const& useKmeans)
{
  List checkin=CheckInWindow_sppmix(points,xlims,ylims,truncate);
  mat data=as<mat>(checkin["data_inW"]);
  int n=data.n_rows;
  Rcout << "Dataset has " << n <<" points" << std::endl ;
  cube genmus=zeros(m,2,L);
  //all elements in the field are 2x2 matrices
  field<mat> gensigmas(L,m),geninvsigmas(L,m);
  mat genz=zeros(L,n),genps=zeros(L,m);
  mat consts=zeros(L,m),approxmass=ones(L,m);
  vec meanps=zeros(m);
  mat meanmus=zeros(m,2);
  cube meansigmas=zeros(2,2,m);
//  bool singularities=false;
  mat sumsig,newdata,propmu=zeros(1,2),
    genmutemp=zeros(1,2),
    mutemp1=zeros(1,2),propsigma,
    sumxmu,ps1,ps2,beta,cov1;
  vec qij(m),approx=zeros(m),
    ds,qs=zeros(m),mu1=zeros(2),
    newmu=zeros(2);
  int i,j,r,dat;
  ivec which=randi(m, distr_param(0,n-1));
  vec mins=zeros(2),maxs=zeros(2);
  if(xlims.size()==1)//passed 0, use data
  {
    mins(0)=min(data.col(0));
    mins(1)=min(data.col(1));
    maxs(0)=max(data.col(0));
    maxs(1)=max(data.col(1));
  }else
  {
    mins(0)=xlims(0);
    mins(1)=ylims(0);
    maxs(0)=xlims(1);
    maxs(1)=ylims(1);
  }
  //  Rcout << "mins="<<mins<< std::endl ;
  double Rx=max(data.col(0))-min(data.col(0)),
    Ry=max(data.col(1))-min(data.col(1));
  vec ksi=zeros(2);
  ksi(0)=sum(data.col(0))/n;
  ksi(1)=sum(data.col(1))/n;
//    Rcout << ksi<<"\n"<< std::endl ;
  mat kappa(2,2),kappainv(2,2),Idenmat(2,2);
  Idenmat(0,1)=0;
  Idenmat(1,0)=0;
  Idenmat(0,0)=1;
  Idenmat(1,1)=1;
  kappa(0,1)=0;
  kappa(1,0)=0;
  kappa(0,0)=100/(Rx*Rx);
  kappa(1,1)=100/(Ry*Ry);
  kappainv(0,1)=0;
  kappainv(1,0)=0;
  //  kappainv.eye(2,2);
  kappainv(0,0)=Rx*Rx/100;
  kappainv(1,1)=Ry*Ry/100;
  //hypers 3:a=3, 4:g=.3, 5:gam=1
  double //a=3,g=1,gam=1;
    a=hyperparams(0),g=hyperparams(1),
      gam=hyperparams(2);
  mat hmat(2,2),hmatinv(2,2);
  hmat(0,1)=0;
  hmat(1,0)=0;
  //  hmat.eye(2,2);
  hmat(0,0)=100*g/(a*Rx*Rx);
  hmat(1,1)=100*g/(a*Ry*Ry);
  //  hmatinv.eye(2,2);
  hmatinv(0,1)=0;
  hmatinv(1,0)=0;
  hmatinv(0,0)=a*Rx*Rx/(100*g);
  hmatinv(1,1)=a*Ry*Ry/(100*g);
  //starting values
//  return List::create();
  sumxmu=zeros(2,2);
  newmu(0)=sum(data.col(0))/n;
  newmu(1)=sum(data.col(1))/n;
  for(r=0;r<n;r++)
    sumxmu=sumxmu+trans(data.row(r)-newmu.t())*(data.row(r)-newmu.t());
  //      Rcout << sumxmu<<"\n"<<std::endl;
  ps2=invmat2d_sppmix(2*Idenmat+sumxmu/(n-1));

  if(useKmeans)
  {
    mat kmeansmat;
//The means parameter is the
//output matrix for storing the
//resulting centroids of the sets,
//with each centroid stored as a
//column vector
    bool status = arma::kmeans(
      kmeansmat, data.t(), m,
      arma::random_subset, 10, false);
    if(status == false)
    {
      stop("clustering failed");
    }
    for (i=0;i<m;i++)
    {
      genmus(i,0,0)=kmeansmat(0,i);
      genmus(i,1,0)=kmeansmat(1,i);
    }
  }
  for (i=0;i<m;i++)
  {
    //    genmus.slice(0);
    if(!useKmeans)
    {
      genmus(i,0,0)=data(which(i),0);
      genmus(i,1,0)=data(which(i),1);
    }
    geninvsigmas(0,i)=rWishart_sppmix(2*a,ps2);//kappa;//invmat2d_sppmix(gensigmas(0,i));
    gensigmas(0,i)=invmat2d_sppmix(geninvsigmas(0,i));//kappainv;
    if(det(gensigmas(0,i))<0.000001)// || det(geninvsigmas(0,i))<.001)
    {
      Rcout <<"\nNear singular covariance for starting value. Det(Sigma)="<<det(gensigmas(0,i))<< std::endl ;
      //      return List::create();
    }//*/
//    gensigmas(0,i)=kappainv;
//    geninvsigmas(0,i)=kappa;//invmat2d_sppmix(gensigmas(0,i));
    genps(0,i)=1.0/m;
    mu1(0)=genmus(i,0,0);
    mu1(1)=genmus(i,1,0);
    if(truncate)
      approxmass(0,i)=
        ApproxBivNormProb_sppmix(xlims,
                               ylims,mu1,gensigmas(0,i),2);
    else
      approxmass(0,i)=1;
    consts(0,i)=1.0/(approxmass(0,i)*sqrt(det(2*datum::pi*gensigmas(0,i))));
  }
  mat prevz,zmultinomial(n,m);
//     Rcout <<"passed"<< std::endl ;
  if(m==1)
    for (dat=0;dat<n;dat++)
    {
      zmultinomial(dat,0)=1;
      genz(0,dat)=0;
    }
  else
  for (dat=0;dat<n;dat++)
  {
    zmultinomial.row(dat)=reshape(rMultinomial_sppmix(1,rDirichlet_sppmix(ones(m))),1,m);
    uvec q1=find(zmultinomial.row(dat)==1);
    genz(0,dat)=q1[0];
//    Rcout <<"q1="<<q1<<"\n"<< std::endl ;
  }
//  return List::create();
  //  Rcout << sum(zmultinomial.col(0))<< std::endl ;
  //  Rcout << sum(zmultinomial.col(1))<< std::endl ;
  //Rcout << zmultinomial<< std::endl ;
  //    return List::create();

  double MHjump=0,ratio=1,sumd;//,approxcompj;//,quad;
  //setup grid for truncation
  int sum1;
/*
  vec ticsx=zeros(LL),ticsy=zeros(LL);
  for (i=0;i<LL;i++)
  {
    ticsx(i)=mins(0)+i*(maxs(0)-mins(0))/(LL-1);
    ticsy(i)=mins(1)+i*(maxs(1)-mins(1))/(LL-1);
  }
  mat areas=zeros(LL,LL);
  for(j=0;j<LL-1;j++)
    for(i=0;i<LL-1;i++)
      areas(i,j)=(ticsx(i+1)-ticsx(i))*(ticsy(j+1)-ticsy(j));
*/
//  int l=2;
//  int countiter=0;
  //start MCMC
  Rcout << "Preliminaries done. Starting MCMC..." << std::endl ;
  for(i=0;i<L-1;i++)
  {
    //sample ksi and kappa
/*    if(1)
    {
      mu1(0)=sum(genmus.slice(i).col(0))/m;
      mu1(1)=sum(genmus.slice(i).col(1))/m;
      cov1=invmat2d_sppmix(m*kappa);
//      Rcout << "passed1"<< std::endl ;
      ksi=trans(rnorm2_sppmix(1,mu1,cov1));
//      Rcout << ksi<<'\n'<< std::endl ;
      sumxmu=zeros(2,2);
      for(r=0;r<m;r++)
        sumxmu=sumxmu+(genmus.slice(i).row(r)-ksi.t()).t()*(genmus.slice(i).row(r)-ksi.t());
      ps1=invmat2d_sppmix(l*Idenmat+sumxmu);
      kappa=rWishart_sppmix(l+m,ps1);
    }*/
    prevz=zmultinomial;
    Rprintf("\rWorking: %3.1f%% complete",100.0*i/(L-2));
    //Rprintf("\rWorking: %3.1f%% complete (iteration %d/%d)",100.0*i/(L-2),i,L);
    //  sample B matrix
    sumsig=zeros(2,2);
    for (j=0;j<m;j++)
      sumsig=sumsig+geninvsigmas(i,j);
    ps1=invmat2d_sppmix(2*hmat+2*sumsig);
    //   if(det(ps1)<=0){Rcout << "\n"<<"ps1 not pd "<<det(ps1)<<"\n"<<std::endl;}
    beta=rWishart_sppmix(2*g+2*m*a,ps1);
    //samples mus and sigmas
    ds=zeros(m,1);
    for(j=0;j<m;j++)
    {
      //samples mus
      if(m>1)
        sum1=sum(zmultinomial.col(j));
      else
        sum1=n;
      //     indi=find(zmultinomial.col(j)==1);
//      Rcout << sum1<<'\n'<< std::endl ;
      if(sum1>0)
        {
          newdata=data.rows(find(zmultinomial.col(j)==1));
          //        Rcout << newdata<< std::endl ;
          //    if(newdata.n_rows!=sum1){Rcout << "\n"<<"not reading data properly "<<newdata.n_rows<< " "<<sum1 << std::endl ;return List::create();}
          newmu(0)=sum(newdata.col(0))/sum1;
          newmu(1)=sum(newdata.col(1))/sum1;
          //        Rcout << newmu << std::endl ;
        }
        else
        {
          //       break;
          // Rcout <<"\n"<< i << " "<< j << " "<<sum1 << std::endl ;
          //         return List::create();
          newmu=zeros(2,1);
          Rcout << "\n"<<"sum1==0, Component with less than 2 points, exiting" << std::endl ;
          return List::create();
        }
      cov1=invmat2d_sppmix(sum1*geninvsigmas(i,j)+kappa);
      //    if(det(cov1)<=0){Rcout << "\n"<<"Cov1 not pd, "<<cov1 << std::endl ;}
      //      Rcout << newmu << sum1<< std::endl ;
      mu1=cov1*(sum1*geninvsigmas(i,j)*newmu+kappa*ksi);
      //      Rcout << mu1<<cov1<< std::endl ;
      //proposed mu
      genmutemp=rnorm2_sppmix(1,mu1,cov1);
      //      Rcout << "genmu" <<genmutemp<< std::endl ;
      if(truncate)
      {
        mu1(0)=genmus(j,0,i);
        mu1(1)=genmus(j,1,i);
        //check if proposed mu is in domain
//        if(mu1(0)<xlims(0) ||mu1(0)>xlims(1) ||
//           mu1(1)<ylims(0) ||mu1(1)>ylims(1))
//          ratio=0;
//        else
          ratio=ApproxMHRatiomu_sppmix(xlims,
            ylims,mu1,trans(genmutemp),
            gensigmas(i,j),sum1);
//        Rcout << ratio << '\n'<< std::endl ;

      //pow(ApproxMHRatiomu_sppmix(LL,xlims,ylims,
      //         mu1,trans(genmutemp),gensigmas(i,j),
        //         geninvsigmas(i,j)),sum1);
      }
      else
        ratio=1;
      if(Rcpp::runif(1)[0]<ratio)
      {
        genmus(j,0,i+1)=genmutemp(0);
        genmus(j,1,i+1)=genmutemp(1);
      }
      else
      {
        genmutemp(0)=genmus(j,0,i);
        genmutemp(1)=genmus(j,1,i);
        genmus(j,0,i+1)=genmus(j,0,i);
        genmus(j,1,i+1)=genmus(j,1,i);
      }
      //samples sigmas
      sumxmu=zeros(2,2);
      if (sum1>0)
        for(r=0;r<sum1;r++)
          sumxmu=sumxmu+(newdata.row(r)-genmutemp).t()*(newdata.row(r)-genmutemp);
      //      Rcout << sumxmu<<"\n"<<std::endl;
      ps2=invmat2d_sppmix(2*beta+sumxmu);
      cov1=rWishart_sppmix(2*a+sum1,ps2);
      //      if(det(ps2)<=0){Rcout << "\n"<<        "ps2 not pd "<<det(ps2)<<"\n"<<        det(beta)<< " "<<det(sumxmu)<<std::endl;}
      propsigma=invmat2d_sppmix(cov1);
      if(truncate)
      {
        mu1(0)=genmus(j,0,i+1);
        mu1(1)=genmus(j,1,i+1);
        ratio=ApproxMHRatiosig_sppmix(xlims,
            ylims,mu1,gensigmas(i,j),
            propsigma,sum1);
      //pow(ApproxMHRatiosig_sppmix(LL,xlims,ylims,
        //      mu1,propsigma,gensigmas(i,j),
          //    geninvsigmas(i,j)),sum1);
      }
      else
        ratio=1;
//      if (det(propsigma)<0.001)
//        ratio=0;
      if(Rcpp::runif(1)[0]<ratio)
      {
        geninvsigmas(i+1,j)=cov1;//invmat2d_sppmix(propsigma);
        gensigmas(i+1,j)=propsigma;
      }
      else
      {
        geninvsigmas(i+1,j)=geninvsigmas(i,j);
        gensigmas(i+1,j)=gensigmas(i,j);
      }
      ds(j)=gam+sum1;
      if (truncate)
      {
        mu1(0)=genmus(j,0,i+1);
        mu1(1)=genmus(j,1,i+1);
        approx[j]=ApproxCompMass_sppmix(
          xlims,ylims,mu1,gensigmas(i+1,j));
      }
      else
        approx[j]=1;
//      approx[j]=approxcompj;
      approxmass(i+1,j)=approx[j];
      consts(i+1,j)=1.0/(approx[j]*
        sqrt(det(2*datum::pi*gensigmas(i+1,j))));
    }
    //sample component probs
    genps.row(i+1)=rDirichlet_sppmix(ds).t();
    //    Rcout << genps.row(i+1)<< std::endl ;
    //    Rcout << ds<< std::endl ;
    //   Rcout << "1 "<<sum(zmultinomial.col(0))<< std::endl ;
    //   Rcout << "2 "<<sum(zmultinomial.col(1))<< std::endl ;
    //    Rcout << consts.row(i+1)<< std::endl ;

//    if(m>1)
      //sample indicators zij
    for(dat=0;dat<n;dat++)
    {
      sumd=0;
      for(j=0;j<m;j++)
      {
        //        genmutemp(0)=genmus(j,0,i+1);
        //       genmutemp(1)=genmus(j,1,i+1);
        //        Rcout << genmutemp<< std::endl ;
        //        Rcout << data.row(dat)<< std::endl ;
        //        mutemp1=data.row(dat)-genmutemp;
        //        quad=as_scalar(mutemp1*geninvsigmas(i+1,j)*trans(mutemp1));
        //        Rcout << quad<< std::endl ;
        //        Rcout << mutemp1<< std::endl ;
/*        mu1(0)=genmus(j,0,i+1);
        mu1(1)=genmus(j,1,i+1);
        qij(j)=genps(i+1,j)*dNormal_sppmix(data.row(dat).t(),
                mu1,gensigmas(i+1,j))/approx[j];*/
/*        if (det(gensigmas(i+1,j))<0.000001)
        {
          Rcout <<i<< std::endl ;
          singularities=true;
        }//*/
        mu1(0)=data(dat,0)-genmus(j,0,i+1);
        mu1(1)=data(dat,1)-genmus(j,1,i+1);
        qij(j)=genps(i+1,j)*consts(i+1,j)*
          exp(-.5*Quad_sppmix(mu1,geninvsigmas(i+1,j)));
 /*       if (qij(j)>1.0f)
        {
          Rcout <<"\ni="<<i<< std::endl ;
          Rcout <<"\ndet(gensigmas(i+1,j))="<<det(gensigmas(i+1,j))<< std::endl ;
//          Rcout <<"\nconsts(i+1,j)="<<consts(i+1,j)<< std::endl ;
          Rcout <<"\ngenps(i+1,j)="<<genps(i+1,j)<< std::endl ;
          Rcout <<"\nquad="<<as_scalar(mutemp1*geninvsigmas(i+1,j)*trans(mutemp1))<< std::endl ;
          Rcout <<"\ndat="<<dat<< std::endl ;
          Rcout <<"\nj="<<j<< std::endl ;
          Rcout <<"\nm="<<m<< std::endl ;
          Rcout <<"\n Error, density gives probability >1. Rerun with truncated set to TRUE, value="<<qij(j)<< std::endl ;
          return List::create();
        }//*/
        sumd=sumd+qij(j);
        //        Rcout << qij(dat,j)<< std::endl ;
      }
  /*    if (sumd>1.0f)
      {
        Rcout <<"\n Error, model gives probability >1. Rerun with truncated set to TRUE, value="<<sumd<< std::endl ;
        return List::create();
      }//*/
      qs=qij/sumd;
      //if(sum(qs)<1)
      //Rcout << sum(qs)<< std::endl ;
      //      if (sum(qij)==1)
      if (m>1)
        zmultinomial.row(dat)=reshape(rMultinomial_sppmix(1,qs),1,m);
      else
        zmultinomial.row(dat)=ones<vec>(m);
//      zmultinomial.row(dat)=reshape(rMultinomial_sppmix(1,qs),1,m);
      //      zmultinomial.row(dat)=reshape(rMultinomial_sppmix(1,qij/sum(qij)),1,m);
//            Rcout << sum(zmultinomial.row(dat))<< std::endl ;
    }
    ratio=1;
    if(m>1)
    for(j=0;j<m;j++)
    {
      //component with less that 2 points is not accepted
      if (sum(zmultinomial.col(j))<2)
      {
          //        Rcout <<i<<zmultinomial<< std::endl ;
          //       Rcout <<"\nComponent "<<j+1<<" has less than 2 pts"<< std::endl ;
          //       Rcout << "1 "<<sum(zmultinomial.col(0))<< std::endl ;
          //       Rcout << "2 "<<sum(zmultinomial.col(1))<< std::endl ;
          //       return List::create();
        ratio=0;
        break;
      }
    }
/*    if(0)//singularities)
    {
      ratio=0;
      singularities=false;
    }*/
    if(Rcpp::runif(1)[0]<ratio)
    {
      MHjump=MHjump+1;
    }
    else
    {
      genps.row(i+1)=genps.row(i);
      for(j=0;j<m;j++)
      {
        geninvsigmas(i+1,j)=geninvsigmas(i,j);
        gensigmas(i+1,j)=gensigmas(i,j);
        genmus(j,0,i+1)=genmus(j,0,i);
        genmus(j,1,i+1)=genmus(j,1,i);
        approxmass(i+1,j)=approxmass(i,j);
      }
      zmultinomial=prevz;
    }
    for(dat=0;dat<n;dat++)
    {
      if(m>1)
      {
        uvec q1=find(zmultinomial.row(dat)==1);
        genz(i+1,dat)=q1[0];
      }
      else
        genz(i+1,dat)=0;
    }
    Rcpp::checkUserInterrupt();
  }
  Rprintf("\rDone. Metropolis-Hastings acceptance: %3.1f%%\n",100.0*MHjump/L);

  //create a list, with each element corresponding
  //to a single realization, which itself is a list
  //with each element containing the mixture ps, mus, sigs,
  //as a list of m elements
  List allgens(L);
  //list containing mixture ps,mus,sigs
  for(i=0;i<L;i++)
  {
    List mix(m);
    for(j=0;j<m;j++)
    {
      mu1(0)=genmus(j,0,i);
      mu1(1)=genmus(j,1,i);
      mix[j]=List::create(
        Named("p") = genps(i,j),
        Named("mu") = mu1,
        Named("sigma") = gensigmas(i,j));
    }
    allgens[i]=mix;
  }

  //sample lambdas
  double alamda=1,blamda=10000;
  vec lamdas=rgamma(L,n+alamda,1/(1+1/blamda));
//  if(singularities)
//    Rcout <<"\nWarning: there were near singular covariance matrices during computation"<< std::endl ;

  return List::create(
    Named("allgens_List") = allgens,
    Named("genps") = genps,
    Named("genmus") = genmus,
    Named("gensigmas") = gensigmas,
    Named("genzs") = genz,
    Named("genlamdas") = lamdas,
    Named("ApproxCompMass")=approxmass);
}

//[[Rcpp::export]]
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
    mat const& startmus)
{
  List checkin=CheckInWindow_sppmix(points,xlims,ylims,truncate);
  mat data=as<mat>(checkin["data_inW"]);
  int n=data.n_rows;
  Rcout << "Dataset has " << n <<" points" << std::endl ;
  cube genmus=zeros(m,2,L);
  //all elements in the field are 2x2 matrices
  field<mat> gensigmas(L,m),geninvsigmas(L,m);
  mat genz=zeros(L,n),genps=zeros(L,m);
  mat consts=zeros(L,m),approxmass=ones(L,m);
  vec meanps=zeros(m);
  mat meanmus=zeros(m,2);
  cube meansigmas=zeros(2,2,m);
  //  bool singularities=false;
  mat sumsig,newdata,propmu=zeros(1,2),
    genmutemp=zeros(1,2),
    mutemp1=zeros(1,2),propsigma,
    sumxmu,ps1,ps2,cov1;
  vec qij(m),approx=zeros(m),
    ds,qs=zeros(m),mu1=zeros(2),
    newmu=zeros(2);
  int i,j,r,dat;
  ivec which=randi(m, distr_param(0,n-1));
  vec mins=zeros(2),maxs=zeros(2);
  if(xlims.size()==1)//passed 0, use data
  {
    mins(0)=min(data.col(0));
    mins(1)=min(data.col(1));
    maxs(0)=max(data.col(0));
    maxs(1)=max(data.col(1));
  }else
  {
    mins(0)=xlims(0);
    mins(1)=ylims(0);
    maxs(0)=xlims(1);
    maxs(1)=ylims(1);
  }
//  Rcout << "mins="<<mins<< std::endl ;
  double Rx=max(data.col(0))-min(data.col(0)),
    Ry=max(data.col(1))-min(data.col(1));
  vec ksi=zeros(2);
  ksi(0)=mu0(0);//sum(data.col(0))/n;
  ksi(1)=mu0(1);//sum(data.col(1))/n;
  //    Rcout << ksi<<"\n"<< std::endl ;
//  Rcout << "\nPass1"<<std::endl;
  mat //kappa(2,2),kappainv(2,2),
    Idenmat(2,2);
  Idenmat(0,1)=0;
  Idenmat(1,0)=0;
  Idenmat(0,0)=1;
  Idenmat(1,1)=1;
  mat beta=sig0*sig0*Idenmat;
  mat kappainv=Sigma0;
  mat kappa=invmat2d_sppmix(Sigma0);
//  kappa(0,1)=0;  kappa(1,0)=0;  kappa(0,0)=100/(Rx*Rx);  kappa(1,1)=100/(Ry*Ry);
//  kappainv(0,1)=0;  kappainv(1,0)=0;  kappainv(0,0)=Rx*Rx/100;  kappainv(1,1)=Ry*Ry/100;
  //hypers 3:a=3, 4:g=.3, 5:gam=1
//  double //a=3,g=1,gam=1;
//    a=hyperparams(0),g=hyperparams(1),gam=hyperparams(2);
//  mat hmat(2,2),hmatinv(2,2);
//  hmat(0,1)=0;  hmat(1,0)=0;  hmat(0,0)=100*g/(a*Rx*Rx);hmat(1,1)=100*g/(a*Ry*Ry);
  //  hmatinv.eye(2,2);
  //hmatinv(0,1)=0;  hmatinv(1,0)=0;  hmatinv(0,0)=a*Rx*Rx/(100*g);  hmatinv(1,1)=a*Ry*Ry/(100*g);
  //starting values
  //  return List::create();
//  Rcout << "\nPass2"<<std::endl;
  sumxmu=zeros(2,2);
  newmu(0)=sum(data.col(0))/n;
  newmu(1)=sum(data.col(1))/n;
  for(r=0;r<n;r++)
    sumxmu=sumxmu+trans(data.row(r)-newmu.t())*(data.row(r)-newmu.t());
//  Rcout << sumxmu<<"\n"<<std::endl;
//Rcout << "\nPass1"<<std::endl;
  ps2=invmat2d_sppmix(2*Idenmat+sumxmu/(n-1));
//Rcout << "\nPass2"<<std::endl;
  mat ksis(m,2);
  if(useKmeans)
  {
    mat kmeansmat;
    //The means parameter is the
    //output matrix for storing the
    //resulting centroids of the sets,
    //with each centroid stored as a
    //column vector
    bool status = arma::kmeans(
      kmeansmat, data.t(), m,
      arma::random_subset, 10, false);
    if(status == false)
    {
      stop("clustering failed");
    }
    for (i=0;i<m;i++)
    {
      genmus(i,0,0)=kmeansmat(0,i);
      genmus(i,1,0)=kmeansmat(1,i);
    }//genmus=zeros(m,2,L)
  }
  for (i=0;i<m;i++)
  {
    //    genmus.slice(0);
    if(!useKmeans)
    {
      genmus(i,0,0)=startmus(i,0);// data(which(i),0);
      genmus(i,1,0)=startmus(i,1);//data(which(i),1);
      ksis.row(i)=startmus.row(i);
    }
    geninvsigmas(0,i)=rWishart_sppmix(df0,ps2);//kappa;//invmat2d_sppmix(gensigmas(0,i));
    gensigmas(0,i)=invmat2d_sppmix(geninvsigmas(0,i));//kappainv;
    if(det(gensigmas(0,i))<0.000001)// || det(geninvsigmas(0,i))<.001)
    {
      Rcout <<"\nNear singular covariance for starting value. Det(Sigma)="<<det(gensigmas(0,i))<< std::endl ;
      //      return List::create();
    }//*/
    //    gensigmas(0,i)=kappainv;
    //    geninvsigmas(0,i)=kappa;//invmat2d_sppmix(gensigmas(0,i));
    mu1(0)=genmus(i,0,0);
    mu1(1)=genmus(i,1,0);
    if(truncate)
      approxmass(0,i)=
        ApproxBivNormProb_sppmix(xlims,
                                 ylims,mu1,gensigmas(0,i),2);
    else
      approxmass(0,i)=1;
    consts(0,i)=1.0/(approxmass(0,i)*sqrt(det(2*datum::pi*gensigmas(0,i))));
  }
  genps.row(0)=rDirichlet_sppmix(d).t();
  mat prevz,zmultinomial(n,m);
  //     Rcout <<"passed"<< std::endl ;
  if(m==1)
    for (dat=0;dat<n;dat++)
    {
      zmultinomial(dat,0)=1;
      genz(0,dat)=0;
    }
    else
      for (dat=0;dat<n;dat++)
      {
 //       Rcout << "\nPass1"<<std::endl;
        zmultinomial.row(dat)=reshape(rMultinomial_sppmix(1,rDirichlet_sppmix(d)),1,m);
//        Rcout << "\nPass12"<<std::endl;
        uvec q1=find(zmultinomial.row(dat)==1);
        genz(0,dat)=q1[0];
        //    Rcout <<"q1="<<q1<<"\n"<< std::endl ;
      }
      //  return List::create();
      //  Rcout << sum(zmultinomial.col(0))<< std::endl ;
      //  Rcout << sum(zmultinomial.col(1))<< std::endl ;
      //Rcout << zmultinomial<< std::endl ;
      //    return List::create();

      double MHjump=0,ratio=1,sumd;//,approxcompj;//,quad;
    //setup grid for truncation
    int sum1;
    /*
    vec ticsx=zeros(LL),ticsy=zeros(LL);
    for (i=0;i<LL;i++)
    {
    ticsx(i)=mins(0)+i*(maxs(0)-mins(0))/(LL-1);
    ticsy(i)=mins(1)+i*(maxs(1)-mins(1))/(LL-1);
    }
    mat areas=zeros(LL,LL);
    for(j=0;j<LL-1;j++)
    for(i=0;i<LL-1;i++)
    areas(i,j)=(ticsx(i+1)-ticsx(i))*(ticsy(j+1)-ticsy(j));
    */
    //  int l=2;
    //  int countiter=0;
    //start MCMC
    Rcout << "Preliminaries done. Starting MCMC..." << std::endl ;
    for(i=0;i<L-1;i++)
    {
      prevz=zmultinomial;
      Rprintf("\rWorking: %3.1f%% complete",100.0*i/(L-2));
      //samples mus and sigmas
      ds=zeros(m,1);
      for(j=0;j<m;j++)
      {
        //samples mus
        if(m>1)
          sum1=sum(zmultinomial.col(j));
        else
          sum1=n;
        //     indi=find(zmultinomial.col(j)==1);
        //      Rcout << sum1<<'\n'<< std::endl ;
        if(sum1>0)
        {
          newdata=data.rows(find(zmultinomial.col(j)==1));
          //        Rcout << newdata<< std::endl ;
          //    if(newdata.n_rows!=sum1){Rcout << "\n"<<"not reading data properly "<<newdata.n_rows<< " "<<sum1 << std::endl ;return List::create();}
          newmu(0)=sum(newdata.col(0))/sum1;
          newmu(1)=sum(newdata.col(1))/sum1;
          //        Rcout << newmu << std::endl ;
        }
        else
        {
          //       break;
          // Rcout <<"\n"<< i << " "<< j << " "<<sum1 << std::endl ;
          //         return List::create();
          newmu=zeros(2,1);
          Rcout << "\n"<<"sum1==0, Component with less than 2 points, exiting" << std::endl ;
          return List::create();
        }
        cov1=invmat2d_sppmix(sum1*geninvsigmas(i,j)+kappa);
        //    if(det(cov1)<=0){Rcout << "\n"<<"Cov1 not pd, "<<cov1 << std::endl ;}
        //      Rcout << newmu << sum1<< std::endl ;
        mu1=cov1*(sum1*geninvsigmas(i,j)*newmu+kappa*ksis.row(j).t());
        //      Rcout << mu1<<cov1<< std::endl ;
        //proposed mu
        genmutemp=rnorm2_sppmix(1,mu1,cov1);
        //      Rcout << "genmu" <<genmutemp<< std::endl ;
        if(truncate)
        {
          mu1(0)=genmus(j,0,i);
          mu1(1)=genmus(j,1,i);
          //check if proposed mu is in domain
          //        if(mu1(0)<xlims(0) ||mu1(0)>xlims(1) ||
          //           mu1(1)<ylims(0) ||mu1(1)>ylims(1))
          //          ratio=0;
          //        else
          ratio=ApproxMHRatiomu_sppmix(xlims,
                                       ylims,mu1,trans(genmutemp),
                                       gensigmas(i,j),sum1);
          //        Rcout << ratio << '\n'<< std::endl ;

          //pow(ApproxMHRatiomu_sppmix(LL,xlims,ylims,
          //         mu1,trans(genmutemp),gensigmas(i,j),
          //         geninvsigmas(i,j)),sum1);
        }
        else
          ratio=1;
        if(Rcpp::runif(1)[0]<ratio)
        {
          genmus(j,0,i+1)=genmutemp(0);
          genmus(j,1,i+1)=genmutemp(1);
        }
        else
        {
          genmutemp(0)=genmus(j,0,i);
          genmutemp(1)=genmus(j,1,i);
          genmus(j,0,i+1)=genmus(j,0,i);
          genmus(j,1,i+1)=genmus(j,1,i);
        }
        //samples sigmas
      //  Rcout << "\nPass2"<<std::endl;
        sumxmu=zeros(2,2);
        if (sum1>0)
          for(r=0;r<sum1;r++)
            sumxmu=sumxmu+(newdata.row(r)-genmutemp).t()*(newdata.row(r)-genmutemp);
        //      Rcout << sumxmu<<"\n"<<std::endl;
        ps2=invmat2d_sppmix(beta+sumxmu);
        cov1=rWishart_sppmix(df0+sum1,ps2);
        //      if(det(ps2)<=0){Rcout << "\n"<<        "ps2 not pd "<<det(ps2)<<"\n"<<        det(beta)<< " "<<det(sumxmu)<<std::endl;}
        propsigma=invmat2d_sppmix(cov1);
        if(truncate)
        {
          mu1(0)=genmus(j,0,i+1);
          mu1(1)=genmus(j,1,i+1);
          ratio=ApproxMHRatiosig_sppmix(xlims,
                                        ylims,mu1,gensigmas(i,j),
                                        propsigma,sum1);
          //pow(ApproxMHRatiosig_sppmix(LL,xlims,ylims,
          //      mu1,propsigma,gensigmas(i,j),
          //    geninvsigmas(i,j)),sum1);
        }
        else
          ratio=1;
        //      if (det(propsigma)<0.001)
        //        ratio=0;
        if(Rcpp::runif(1)[0]<ratio)
        {
          geninvsigmas(i+1,j)=cov1;//invmat2d_sppmix(propsigma);
          gensigmas(i+1,j)=propsigma;
        }
        else
        {
          geninvsigmas(i+1,j)=geninvsigmas(i,j);
          gensigmas(i+1,j)=gensigmas(i,j);
        }
        ds(j)=d(j)+sum1;
        if (truncate)
        {
          mu1(0)=genmus(j,0,i+1);
          mu1(1)=genmus(j,1,i+1);
          approx[j]=ApproxCompMass_sppmix(
            xlims,ylims,mu1,gensigmas(i+1,j));
        }
        else
          approx[j]=1;
        //      approx[j]=approxcompj;
        approxmass(i+1,j)=approx[j];
        consts(i+1,j)=1.0/(approx[j]*
          sqrt(det(2*datum::pi*gensigmas(i+1,j))));
      }
      //sample component probs
      genps.row(i+1)=rDirichlet_sppmix(ds).t();
      //    Rcout << genps.row(i+1)<< std::endl ;
      //    Rcout << ds<< std::endl ;
      //   Rcout << "1 "<<sum(zmultinomial.col(0))<< std::endl ;
      //   Rcout << "2 "<<sum(zmultinomial.col(1))<< std::endl ;
      //    Rcout << consts.row(i+1)<< std::endl ;

      //    if(m>1)
      //sample indicators zij
      for(dat=0;dat<n;dat++)
      {
        sumd=0;
        for(j=0;j<m;j++)
        {
          //        genmutemp(0)=genmus(j,0,i+1);
          //       genmutemp(1)=genmus(j,1,i+1);
          //        Rcout << genmutemp<< std::endl ;
          //        Rcout << data.row(dat)<< std::endl ;
          //        mutemp1=data.row(dat)-genmutemp;
          //        quad=as_scalar(mutemp1*geninvsigmas(i+1,j)*trans(mutemp1));
          //        Rcout << quad<< std::endl ;
          //        Rcout << mutemp1<< std::endl ;
          /*        mu1(0)=genmus(j,0,i+1);
          mu1(1)=genmus(j,1,i+1);
          qij(j)=genps(i+1,j)*dNormal_sppmix(data.row(dat).t(),
              mu1,gensigmas(i+1,j))/approx[j];*/
          /*        if (det(gensigmas(i+1,j))<0.000001)
          {
          Rcout <<i<< std::endl ;
          singularities=true;
          }//*/
          mu1(0)=data(dat,0)-genmus(j,0,i+1);
          mu1(1)=data(dat,1)-genmus(j,1,i+1);
          qij(j)=genps(i+1,j)*consts(i+1,j)*
            exp(-.5*Quad_sppmix(mu1,geninvsigmas(i+1,j)));
          /*       if (qij(j)>1.0f)
          {
          Rcout <<"\ni="<<i<< std::endl ;
          Rcout <<"\ndet(gensigmas(i+1,j))="<<det(gensigmas(i+1,j))<< std::endl ;
          //          Rcout <<"\nconsts(i+1,j)="<<consts(i+1,j)<< std::endl ;
          Rcout <<"\ngenps(i+1,j)="<<genps(i+1,j)<< std::endl ;
          Rcout <<"\nquad="<<as_scalar(mutemp1*geninvsigmas(i+1,j)*trans(mutemp1))<< std::endl ;
          Rcout <<"\ndat="<<dat<< std::endl ;
          Rcout <<"\nj="<<j<< std::endl ;
          Rcout <<"\nm="<<m<< std::endl ;
          Rcout <<"\n Error, density gives probability >1. Rerun with truncated set to TRUE, value="<<qij(j)<< std::endl ;
          return List::create();
          }//*/
          sumd=sumd+qij(j);
          //        Rcout << qij(dat,j)<< std::endl ;
        }
        /*    if (sumd>1.0f)
        {
        Rcout <<"\n Error, model gives probability >1. Rerun with truncated set to TRUE, value="<<sumd<< std::endl ;
        return List::create();
        }//*/
        qs=qij/sumd;
        //if(sum(qs)<1)
        //Rcout << sum(qs)<< std::endl ;
        //      if (sum(qij)==1)
        if (m>1)
          zmultinomial.row(dat)=reshape(rMultinomial_sppmix(1,qs),1,m);
        else
          zmultinomial.row(dat)=ones<vec>(m);
        //      zmultinomial.row(dat)=reshape(rMultinomial_sppmix(1,qs),1,m);
        //      zmultinomial.row(dat)=reshape(rMultinomial_sppmix(1,qij/sum(qij)),1,m);
        //            Rcout << sum(zmultinomial.row(dat))<< std::endl ;
}
      ratio=1;
      if(m>1)
        for(j=0;j<m;j++)
        {
          //component with less that 2 points is not accepted
          if (sum(zmultinomial.col(j))<2)
          {
            //        Rcout <<i<<zmultinomial<< std::endl ;
            //       Rcout <<"\nComponent "<<j+1<<" has less than 2 pts"<< std::endl ;
            //       Rcout << "1 "<<sum(zmultinomial.col(0))<< std::endl ;
            //       Rcout << "2 "<<sum(zmultinomial.col(1))<< std::endl ;
            //       return List::create();
            ratio=0;
            break;
          }
        }
        /*    if(0)//singularities)
        {
        ratio=0;
        singularities=false;
        }*/
        if(Rcpp::runif(1)[0]<ratio)
        {
          MHjump=MHjump+1;
        }
        else
        {
          genps.row(i+1)=genps.row(i);
          for(j=0;j<m;j++)
          {
            geninvsigmas(i+1,j)=geninvsigmas(i,j);
            gensigmas(i+1,j)=gensigmas(i,j);
            genmus(j,0,i+1)=genmus(j,0,i);
            genmus(j,1,i+1)=genmus(j,1,i);
            approxmass(i+1,j)=approxmass(i,j);
          }
          zmultinomial=prevz;
        }
        for(dat=0;dat<n;dat++)
        {
          if(m>1)
          {
            uvec q1=find(zmultinomial.row(dat)==1);
            genz(i+1,dat)=q1[0];
          }
          else
            genz(i+1,dat)=0;
        }
        Rcpp::checkUserInterrupt();
        }
    Rprintf("\rDone. Metropolis-Hastings acceptance: %3.1f%%\n",100.0*MHjump/L);

    //create a list, with each element corresponding
    //to a single realization, which itself is a list
    //with each element containing the mixture ps, mus, sigs,
    //as a list of m elements
    List allgens(L);
    //list containing mixture ps,mus,sigs
    for(i=0;i<L;i++)
    {
      List mix(m);
      for(j=0;j<m;j++)
      {
        mu1(0)=genmus(j,0,i);
        mu1(1)=genmus(j,1,i);
        mix[j]=List::create(
          Named("p") = genps(i,j),
          Named("mu") = mu1,
          Named("sigma") = gensigmas(i,j));
      }
      allgens[i]=mix;
    }

    //sample lambdas
    double alamda=1,blamda=10000;
    vec lamdas=rgamma(L,n+alamda,1/(1+1/blamda));
    //  if(singularities)
    //    Rcout <<"\nWarning: there were near singular covariance matrices during computation"<< std::endl ;

    return List::create(
      Named("allgens_List") = allgens,
      Named("genps") = genps,
      Named("genmus") = genmus,
      Named("gensigmas") = gensigmas,
      Named("genzs") = genz,
      Named("genlamdas") = lamdas,
      Named("ApproxCompMass")=approxmass);
    }
