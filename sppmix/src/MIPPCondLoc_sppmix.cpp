#include "sppmix.h"
//Written by Sakis Micheas, 2016
//data (locations), marks, x-y limits,
//m=num of comps to fit
//L=num iter
//just get realizations no plotting here

mat GetNeighborMarks(
    mat const& data,
    vec const& marks,
    double const& r,
    vec const& uniquemarks,
    int const& marknum)
{
  int i,j,k,n=data.n_rows;
  mat neighbors=zeros(marknum,n);
  for(i=0;i<n;i++)
  {
    //count neighbors
    double summark=0;
    int ncount=0;
    for(j=0;j<n;j++)
    {
//      if(i==j)continue;
      double dist1=VecLen2(
        trans(data.row(j)-data.row(i)));
      if(dist1<=r)
      {
        ncount++;
        summark+=marks(j);
      }
    }
    for(k=0;k<marknum;k++)
      neighbors(k,i)=
        SQ_sppmix(uniquemarks(k)-
        summark/ncount);
  }
  return(neighbors);
}

//[[Rcpp::export]]
List MIPPCondLoc_sppmix(mat const& points,
                        vec const& marks,
                     vec const& xlims,
                     vec const& ylims,
                     int const& L,
                     bool const& truncate,
                     vec const& hyperparams,
                     vec const& uniquemarks,
                     bool const& discrete_mark,
                     double const& r)
{
  List checkin=CheckInWindow_sppmix(points,xlims,ylims,truncate);
  mat data=as<mat>(checkin["data_inW"]);
  int i,j,k,m,n=data.n_rows;
  Rcout << "\nDataset has " << n <<" points" << std::endl ;

  //find neighborhoods
  int marknum=uniquemarks.size();
  mat gen_gammas=zeros(L,marknum),
    iden=arma::eye(marknum,marknum),
    neighbors=zeros(marknum,n);
  neighbors=GetNeighborMarks(data,marks,
                   r,uniquemarks,marknum);
  //hyperparams[0] is the std
  double MHjump=0,
    sigma0=hyperparams[0];
//  vec mu0=zeros(2);
//  mat sigma0=hyperparams[0]*arma::eye(marknum,marknum);
  vec v1=hyperparams.subvec(1,marknum);
//  Rcout<<"start value="<<v1<<"\n";
  gen_gammas.row(0)=v1.t();
  //norm(marknum,mu0,sigma0);
  // start the main loop
  for (i=0;i<L-1;i++)
  {
    Rcout<<"\rWorking: "<<100.0*(i+1)/L<<"% complete          ";
    //    Rprintf("\rWorking: %3.1f%% complete",100.0*i/(L-2));
    gen_gammas.row(i+1)=gen_gammas.row(i);
    //individual updates
    for(m=0;m<marknum;m++)
    {
      double numer,denumer,prod1=1,
        ratio=1,totcur,totprop,
        prodcur=1,prodprop=1;
      gen_gammas(i+1,m)=gen_gammas(i,m)+sigma0*arma::randn(1)[0];
      for(k=0;k<n;k++)
      {
        totcur=totprop=0;
        for(j=0;j<marknum;j++)
        {
          totcur+=exp(-gen_gammas(i,j)*neighbors(j,k));
          totprop+=exp(-gen_gammas(i+1,j)*neighbors(j,k));
        }
        for(j=0;j<marknum;j++)
          if(marks(k)==uniquemarks(j))
          {
            prodcur*=exp(-gen_gammas(i,j)*neighbors(j,k))/totcur;
            prodprop*=exp(-gen_gammas(i+1,j)*neighbors(j,k))/totprop;
            break;
          }
      }
//      numer=exp(-(.5/sigma0)*VecNorm2(gen_gammas.row(i+1).t()));
//      denumer=exp(-(.5/sigma0)*VecNorm2(gen_gammas.row(i).t()));
      numer=exp(-(.5/1)*SQ_sppmix(gen_gammas(i+1,m)));
      denumer=exp(-(.5/1)*SQ_sppmix(gen_gammas(i,m)));
      prod1=numer/denumer;
      ratio=prod1*prodprop/prodcur;
      if (Rcpp::runif(1)[0]<ratio)
        MHjump=MHjump+1;
      else
        gen_gammas(i+1,m)=gen_gammas(i,m);
    }
//block update
/*
    double numer,denumer,prod1=1,
      ratio=1,totcur,totprop,
  //    prodcur=1,prodprop=1,prods=1,
      logsum=0;
    v1=gen_gammas.row(i).t()+sigma0*arma::randn(marknum);
    gen_gammas.row(i+1)=v1.t();
    for(k=0;k<n;k++)
    {
//      Rcout<<"neighbors k="<<k
//           <<"\n"<<neighbors.col(k);
      totcur=totprop=0;
      for(j=0;j<marknum;j++)
      {
        totcur+=exp(-gen_gammas(i,j)*neighbors(j,k));
        totprop+=exp(-gen_gammas(i+1,j)*neighbors(j,k));
      }
      for(j=0;j<marknum;j++)
        if(marks(k)==uniquemarks(j))
        {
//          prodcur*=exp(-gen_gammas(i,j)*neighbors(j,k))/totcur;
  //        prodprop*=exp(-gen_gammas(i+1,j)*neighbors(j,k))/totprop;
//prods*=totcur*exp((-gen_gammas(i+1,j)+gen_gammas(i,j))*neighbors(j,k))/totprop;//prodprop/prodcur
          logsum+=totcur-totprop+(-gen_gammas(i+1,j)+gen_gammas(i,j))*neighbors(j,k);
          break;
        }
    }
//    Rcout<<"prods="<<prods<<"\n";
    numer=exp(-(.5/(sigma0))*VecNorm2(gen_gammas.row(i+1).t()));
    denumer=exp(-(.5/(sigma0))*VecNorm2(gen_gammas.row(i).t()));
    prod1=numer/denumer;
//    Rcout//<<"prod1="<<prod1<<"\n"
//         <<"prodprop/prodcur="<<prodprop/prodcur<<"\n";
    ratio=prod1*exp(logsum);//prods;//prodprop/prodcur;
    if (Rcpp::runif(1)[0]<ratio)
      MHjump=MHjump+1;
    else
      gen_gammas.row(i+1)=gen_gammas.row(i);
*/
    Rcpp::checkUserInterrupt();
  }
  Rprintf("\rDone. Metropolis-Hastings acceptance: %3.1f%%                 \n",100.0*MHjump/(marknum*L));

//  for (i=0;i<L;i++)
//    gen_gammas.row(i)=arma::sort(gen_gammas.row(i).t()).t();
  return List::create(
    Named("gen_gammas") = gen_gammas,
 //   Named("prob_field") = ps_fields,
 //   Named("MHjump") = MHjump/L,
    Named("discrete_mark") = discrete_mark);
}

//[[Rcpp::export]]
List GetProbFieldsCondLoc_sppmix(
    mat const& points,vec const& marks,
    vec const& xlims,vec const& ylims,
    int const& LL,vec const& meangamma,
    vec const& uniquemarks,
    bool const& truncate,double const& r)
{
  List checkin=CheckInWindow_sppmix(points,xlims,ylims,truncate,false);
  mat data=as<mat>(checkin["data_inW"]);
  int i,j,k,dat,n=data.n_rows;
  vec ticsx=zeros(LL),ticsy=zeros(LL);
  for (i=0;i<LL;i++)
  {
    ticsx(i)=xlims(0)+i*(xlims(1)-xlims(0))/(LL-1);
    ticsy(i)=ylims(0)+i*(ylims(1)-ylims(0))/(LL-1);
  }
  //find neighborhoods
  int marknum=uniquemarks.size();
  cube ps_fields=zeros(LL,LL,marknum);
  for(i=0;i<LL;i++)
    for(j=0;j<LL;j++)
    {
      vec pt(2);
      pt(0)=ticsx(i);
      pt(1)=ticsy(j);
      double summark=0;
      int ncount=0;
      for(dat=0;dat<n;dat++)
      {
        double dist1=VecLen2(
          trans(data.row(dat)-pt.t()));
        if(dist1<=r)
        {
          ncount++;
          summark+=marks(dat);
        }
      }
      double sum1=0;
      if(ncount>0)
        for(k=0;k<marknum;k++)
          sum1+=exp(-meangamma(k)*
            SQ_sppmix(uniquemarks(k)-
            summark/ncount));
      else
        sum1=marknum;
      for(k=0;k<marknum;k++)
      {
        if(ncount>0)
          ps_fields(i,j,k)=exp(-meangamma(k)*
            SQ_sppmix(uniquemarks(k)-
            summark/ncount))/sum1;
        else
          ps_fields(i,j,k)=0;
//        if(ps_fields(i,j,k)>1)ps_fields(i,j,k)=1;
//        else if(ps_fields(i,j,k)<0)ps_fields(i,j,k)=0;
      }
    }
  return List::create(
    Named("x") = ticsx,
    Named("y") = ticsy,
    Named("ps_fields") = ps_fields);
}

//[[Rcpp::export]]
List GenMarksProbCondLoc_sppmix(
    mat const& points,
    int const& L,
    vec const& xlims,vec const& ylims,
    vec const& meangamma,
    vec const& uniquemarks,
    bool const& truncate,double const& r)
{
  List checkin=CheckInWindow_sppmix(points,xlims,ylims,truncate,false);
  mat data=as<mat>(checkin["data_inW"]);
  int i,k,j,l,n=data.n_rows;
  vec marks(n);
  //find neighborhoods
  int marknum=uniquemarks.size();
  double sum1=0;
  mat probs=zeros(n,marknum);
  //get starting values, discrete uniform
  marks(0)=rDiscrete_sppmix(1,ones(marknum)/marknum);
  Rprintf("start");
  for(i=1;i<n;i++)
  {
    //count neighbors
    int ncount=0;
    double summark=0;
    for(j=0;j<i;j++)
    {
      double dist1=VecLen2(
        trans(data.row(j)-data.row(i)));
      if(dist1<=r)
      {
        ncount++;
        summark+=marks(j);
      }
    }
    if(ncount==0)
    {
      marks(i)=rDiscrete_sppmix(1,ones(marknum)/marknum);
      continue;
    }
    sum1=0;
    for(k=0;k<marknum;k++)
      sum1+=exp(-meangamma(k)*
       SQ_sppmix(uniquemarks(k)-
        summark/ncount));
    for(k=0;k<marknum;k++)
      probs(i,k)=exp(-meangamma(k)*
       SQ_sppmix(uniquemarks(k)-
        summark/ncount))/sum1;
    marks(i)=rDiscrete_sppmix(1,probs.row(i).t());
  }

  for(l=0;l<L;l++)
  {
    Rprintf("\rWorking: %3.1f%% complete",100.0*l/(L-1));
    for(i=0;i<n;i++)
    {
      //count neighbors
//      mat neighbors=zeros(marknum);
      double summark=0;
      int ncount=0;
      for(j=0;j<n;j++)
      {
        double dist1=VecLen2(
          trans(data.row(j)-data.row(i)));
        if(dist1<=r)
        {
          ncount++;
          summark+=marks(j);
        }
      }
      if(ncount==0)
      {
        marks(i)=rDiscrete_sppmix(1,ones(marknum)/marknum);
        continue;
      }
      sum1=0;
      for(k=0;k<marknum;k++)
        sum1+=exp(-meangamma(k)*
         SQ_sppmix(uniquemarks(k)-
          summark/ncount));
      for(k=0;k<marknum;k++)
        probs(i,k)=exp(-meangamma(k)*
         SQ_sppmix(uniquemarks(k)-
          summark/ncount))/sum1;
      marks(i)=rDiscrete_sppmix(1,probs.row(i).t());
    }
    Rcpp::checkUserInterrupt();
  }
  Rprintf("\rDone.                                          ");
  return  List::create(
    Named("probs") = probs,
    Named("marks") = marks);
}

//[[Rcpp::export]]
List GetProbCondLoc_sppmix(
    mat const& points,vec const& origmarks,
    vec const& xlims,vec const& ylims,
    vec const& meangamma,
    vec const& uniquemarks,
    bool const& truncate,double const& r)
{
  List checkin=CheckInWindow_sppmix(points,xlims,ylims,truncate,false);
  mat data=as<mat>(checkin["data_inW"]);
  int i,k,j,n=data.n_rows,LL=100;
  vec marks=origmarks;
  List probfields=GetProbFieldsCondLoc_sppmix(
    data,marks,xlims,ylims,LL,
    meangamma,uniquemarks,truncate,r);
  vec ticsx=probfields[0],
      ticsy=probfields[1];
  cube ps_fields=probfields[2];
  int marknum=uniquemarks.size();
  mat probs=zeros(n,marknum);
  for(i=0;i<n;i++)
  {
    mat neighbors=zeros(marknum);
    int indx=0,indy=0;
    for(j=1;j<LL;j++)
    {
      if(ticsx(j-1)<=data(i,0) &&
         data(i,0)<ticsx(j))
      {
        indx=j-1;
        break;
      }
    }
    for(j=1;j<LL;j++)
    {
      if(ticsy(j-1)<=data(i,1) &&
         data(i,1)<ticsy(j))
      {
        indy=j-1;
        break;
      }
    }
    for(k=0;k<marknum;k++)
      probs(i,k)=ps_fields(indx,indy,k);
  }
  return  List::create(
    Named("probs") = probs);
}
