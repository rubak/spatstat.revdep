#include "sppmix.h"
//Written by Sakis Micheas, 2017
//Operations on posterior realizations
//used by cpp and R functions
//visible in the package only

//[[Rcpp::export]]
vec GetPriorVals_sppmix(mat const& pp,
    List const& allgens,
    int const& priortype,
    vec const& d,vec const& mu0,
    mat const& Sigma0,
    int const& df0,
    double const& sig0)
{
//for each realization, compute and return the
//prior density at those values
  int i,j,k,L=allgens.size();
  List mix1=allgens[0];
  int m=mix1.size(),n=pp.n_rows;
  mat allps=GetAllRealiz_ps_sppmix(allgens);
  List allmus=GetAllRealiz_mus_sppmix(allgens);
  List allsigmas=GetAllRealiz_sigmas_sppmix(allgens);
//  pattern <- cbind(pp$x,pp$y)
  vec priorvals(L),
    ppx=pp.col(0),ppy=pp.col(1);
  double Rx=arma::max(ppx)-arma::min(ppx),
    Ry=arma::max(ppy)-arma::min(ppy);
  vec ksi(2);
  ksi(0)=arma::mean(ppx);
  ksi(1)=arma::mean(ppy);
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
  kappainv(0,0)=Rx*Rx/100;
  kappainv(1,1)=Ry*Ry/100;
//  mat sumxmu=zeros(2,2);
//  for(k=0;k<n;k++)
//    sumxmu=sumxmu+(pp.row(k).t()-ksi)*(pp.row(k)-ksi.t());
  mat ps2=//invmat2d_sppmix(
    sig0*sig0*Idenmat;//);//+sumxmu/(n-1));
  for(i=0;i<L;i++)
  {
  //  Rprintf("\rComputing prior values... %3.1f%% complete",100.0*(i+1)/L);
    vec cur_ps=allps.row(i).t();//GetRealiz_ps_sppmix(allgens,i);
    mat cur_mus=allmus[i];//GetRealiz_mus_sppmix(allgens,i);
    mat cur_sigmas=allsigmas[i];//GetRealiz_sigmas_sppmix(allgens,i);
    double dnorm=dDirichlet_sppmix(cur_ps,d),dwish=1;
    for(j=0;j<m;j++)
    {
      vec mu1=cur_mus.row(j).t();
      vec sig1=cur_sigmas.row(j).t();
      mat W(2,2);
      W(0,0)=sig1(0);
      W(0,1)=sig1(1);
      W(1,0)=sig1(2);
      W(1,1)=sig1(3);
      dnorm*=dNormal_sppmix(
        mu1,ksi,kappainv);
      dwish*=dInvWishart_sppmix(
          W,df0,ps2);
    }
    priorvals(i)=dnorm*dwish;
  }
//  Rprintf("\rDone                                               ");
  return(priorvals);
}

// [[Rcpp::export]]
List GetStats_sppmix(vec const& gens,
                     double const& alpha)
{
  //apply burnin before calling this function
  double mu=mean(gens);
  int L=gens.size();
  vec CS(2);
  vec sortedgens=sort(gens);
  //  Rcout << sortedgens<< std::endl ;
  CS(1)=sortedgens((int)floor((1-alpha/2.0)*L));
  CS(0)=sortedgens((int)floor(alpha*L/2.0));
  return List::create(
    Named("Min") = gens.min(),
    Named("Mean") = mu,
    Named("Max") = gens.max(),
    Named("CredibleSet") = CS,
    Named("CredibleSetConfidence") = 100*(1-alpha));
}

// [[Rcpp::export]]
mat GetAllRealiz_ps_sppmix(List const& allgens)
{
  int i,j,L=allgens.size();
  List mixcomp,gen_realiz=allgens[0];
  int m=gen_realiz.size();
  mat allps(L,m);
  for(i=0;i<L;i++)
  {
    gen_realiz=allgens[i];
    for(j=0;j<m;j++)
    {
      mixcomp=gen_realiz[j];
      allps(i,j)=as<double>(mixcomp["p"]) ;
    }
  }
  return allps;
}

// [[Rcpp::export]]
List GetAllRealiz_mus_sppmix(List const& allgens)
{
  //list of L mx2 matrices
  int i,j,L=allgens.size();
  List mixcomp,gen_realiz=allgens[0];
  int m=gen_realiz.size();
  List allmus(L);
  for(i=0;i<L;i++)
  {
    gen_realiz=allgens[i];
    mat musmat(m,2);
    for(j=0;j<m;j++)
    {
      mixcomp=gen_realiz[j];
      musmat.row(j)=trans(as<vec>(mixcomp["mu"]));
    }
    allmus[i]=musmat;
  }
  return allmus;
}

// [[Rcpp::export]]
List GetAllRealiz_sigmas_sppmix(List const& allgens)
{
  //list of L mx4 matrices
  int i,j,L=allgens.size();
  List mixcomp,gen_realiz=allgens[0];
  int m=gen_realiz.size();
  List allsigs(L);
  for(i=0;i<L;i++)
  {
    gen_realiz=allgens[i];
    mat sigmat(m,4);
    for(j=0;j<m;j++)
    {
      mixcomp=gen_realiz[j];
      mat sigmak = as<mat>(mixcomp["sigma"]);
      sigmat(j,0)=sigmak(0,0);
      sigmat(j,1)=sigmak(0,1);
      sigmat(j,2)=sigmak(1,0);
      sigmat(j,3)=sigmak(1,1);
    }
    allsigs[i]=sigmat;
  }
  return allsigs;
}

// [[Rcpp::export]]
vec GetRealiz_ps_sppmix(List const& allgens,
                        int const& realiz)
{
  int j,L=allgens.size();
  if(realiz>=L)
  {
    Rcout << "index out of bounds" << std::endl ;
    return 0;
  }
  List mixcomp,gen_realiz=allgens[realiz];
  int m=gen_realiz.size();
  vec ps(m);
  for(j=0;j<m;j++)
  {
    mixcomp=gen_realiz[j];
    ps(j)=as<double>(mixcomp["p"]) ;
  }
  return ps;
}


// [[Rcpp::export]]
mat GetRealiz_mus_sppmix(List const& allgens,
                         int const& realiz)
{
  int j,L=allgens.size();
  if(realiz>=L)
  {
    Rcout << "index out of bounds" << std::endl ;
    return 0;
  }
  List mixcomp,gen_realiz=allgens[realiz];
  int m=gen_realiz.size();
  mat mus(m,2);
  for(j=0;j<m;j++)
  {
    mixcomp=gen_realiz[j];
    mus.row(j)= trans(as<vec>(mixcomp["mu"]));
  }
  return mus;
}


// [[Rcpp::export]]
mat GetRealiz_sigmas_sppmix(List const& allgens,
                            int const& realiz)
{
  int j,L=allgens.size();
  if(realiz>=L)
  {
    Rcout << "index out of bounds" << std::endl ;
    return 0;
  }
  List mixcomp,gen_realiz=allgens[realiz];
  int m=gen_realiz.size();
  //vectorize the matrix and return it
  mat sigmas(m,4);
  for(j=0;j<m;j++)
  {
    mixcomp=gen_realiz[j];
    mat sigmak = as<mat>(mixcomp["sigma"]);
    sigmas(j,0)=sigmak(0,0);
    sigmas(j,1)=sigmak(0,1);
    sigmas(j,2)=sigmak(1,0);
    sigmas(j,3)=sigmak(1,1);
    //  Rcout << muk<< std::endl ;
  }
  return sigmas;
}


// [[Rcpp::export]]
List PostGenGetBestPerm_sppmix(List const& allgens)
{
//sppmix::PostGenGetBestPerm_sppmix(gens$allgens)
  int i,j,k,L=allgens.size();
  double mind,loss1;
  List permuted_gens(L),mix1=allgens[0];
  int m=mix1.size();
  double permnum=Factorial_sppmix(m);
  mat //stateperm=zeros(L,m),
    allperms=GetAllPermutations_sppmix(m),
    current_perm=zeros(L,m),//will contain the best perm
    previous_perm=ones(L,m);
//  vec oldindex(m),found(m);
  int minindex=0,done=0,count=0;
//  vec action(m);//estimating the ds hyperparam
//repeat until we find the L permutations
//that yield the smallest MC risk wrt choice of
//hyperparams ds for ps
//loss function is:-log(Dirichlet(ps,ds))
  vec cur_sig(4),mu1(2);//,mu11(2);
//  printf("\nNote: Slow operation for moderate-large number of components\n");
//  Rcout << "\nm="<<m<< std::endl ;
  mat allps=GetAllRealiz_ps_sppmix(allgens);
  List allmus=GetAllRealiz_mus_sppmix(allgens);
  List allsigmas=GetAllRealiz_sigmas_sppmix(allgens);
  //    Rcout << "\nm1="<<m<< std::endl ;
  while (done==0)
  {
    //squared error loss, minimized at mean(thetas)
    vec sum_p=zeros(m);
    mat sum_mu=zeros(m,2);
    mat sum_sigma=zeros(m,4);
    current_perm=previous_perm;
    for(i=0;i<L;i++)
    {
//      dDirichlet_sppmix(perm_ps,action);
      vec v1=current_perm.row(i).t(); //rPerm_sppmix(m);
//      current_perm.row(i)=v1.t();
      //permute all generated p's,mu's and sigma's
      //and then get the average
      vec cur_ps=allps.row(i).t();//GetRealiz_ps_sppmix(allgens,i);
      //      Rcout << cur_ps << std::endl ;
      mat cur_mus=allmus[i];//GetRealiz_mus_sppmix(allgens,i);
      mat cur_sigmas=allsigmas[i];//GetRealiz_sigmas_sppmix(allgens,i);
      vec perm_ps=Permute_vec_sppmix(cur_ps,v1);
      mat perm_mus=Permute_mat_sppmix(cur_mus,v1);
      mat perm_sigmas=Permute_mat_sppmix(cur_sigmas,v1);
      sum_p=sum_p+perm_ps;
      sum_mu=sum_mu+perm_mus;
      sum_sigma=sum_sigma+perm_sigmas;
    }
    sum_p=sum_p/L;
    sum_mu=sum_mu/L;
    sum_sigma=sum_sigma/L;
//    Rcout << sum_p << std::endl ;
//    Rcout << sum_mu << std::endl ;
//    Rcout << sum_sigma << std::endl ;
    done=1;
//    Rcout << "passed1" << std::endl ;
    for(i=0;i<L;i++)
    {
//      vec v1=current_perm.row(i).t();
      vec cur_ps=allps.row(i).t();//GetRealiz_ps_sppmix(allgens,i);
      //      Rcout << cur_ps << std::endl ;
      mat cur_mus=allmus[i];//GetRealiz_mus_sppmix(allgens,i);
      mat cur_sigmas=allsigmas[i];//GetRealiz_sigmas_sppmix(allgens,i);
      //go through all permutations and find the one
      //that minimizes the loss for the given action
      //there are m! permutations...
      mind=1000000000000;
      for(k=0;k<permnum;k++)
      {
        vec v2=allperms.row(k).t();// GetAPermutation_sppmix(m,k);
        vec perm_ps=Permute_vec_sppmix(cur_ps,v2);
        mat perm_mus=Permute_mat_sppmix(cur_mus,v2);
        mat perm_sigmas=Permute_mat_sppmix(cur_sigmas,v2);
        loss1=0;
        for(j=0;j<m;j++)
        {
          mu1=sum_mu.row(j).t()-perm_mus.row(j).t();
          cur_sig=sum_sigma.row(j).t()-perm_sigmas.row(j).t();
//          mu11=cur_mus.row(j).t();
//          if(perm_ps(j)<.000001)
//            loss1=loss1+10000000000;
//          else
          loss1=loss1
//            +sqrt(mu11(0)*mu11(0)+mu11(1)*mu11(1))
            +VecNorm2(mu1)
            +dot(cur_sig,cur_sig)
            +SQ_sppmix(sum_p(j)-perm_ps(j));
        }
        if(loss1<mind)
        {
          mind=loss1;
          minindex=k;
        }
      }
      current_perm.row(i)=allperms.row(minindex);//GetAPermutation_sppmix(m,minindex).t();
      Rcpp::checkUserInterrupt();
    }

    for(i=0;i<L;i++)
    {
      for(j=0;j<m;j++)
      {
        if (current_perm(i,j)!=previous_perm(i,j))
        {
          done=0;
          break;
        }
      }
      if(done==0)break;
    }
//    oldindex=stateperm.row(minindex).t();
    previous_perm=current_perm;
    if(done)
    {
//      found=oldindex;
      break;
    }
    count++;
//    printf("\rApplying Permutations, iteration %d",count);
    if (count>100000)
    {
      Rcout << "Didnt find it in 100000 iterations"<< std::endl ;
      break;
    }
  }
//  Rcout << "Best permutation:"<<current_perm//found.t()
//    << std::endl ;
//apply this permuation to all of the original
//realizations and return the permuted
// ps,mus and sigmas
  cube permgenmus=zeros(m,2,L);
  field<mat> permgensigmas(L,m);
  mat permgenps=zeros(L,m);

  for(i=0;i<L;i++)
  {
    vec cur_ps=allps.row(i).t();//GetRealiz_ps_sppmix(allgens,i);
    mat cur_mus=allmus[i];//GetRealiz_mus_sppmix(allgens,i);
    mat cur_sigmas=allsigmas[i];//GetRealiz_sigmas_sppmix(allgens,i);
    vec mu1(2),perm_ps=Permute_vec_sppmix(cur_ps,current_perm.row(i).t());
    mat perm_mus=Permute_mat_sppmix(cur_mus,current_perm.row(i).t());
    mat sig=zeros(2,2),perm_sigmas=Permute_mat_sppmix(cur_sigmas,current_perm.row(i).t());
//    vec mu1(2),perm_ps=Permute_vec_sppmix(cur_ps,found);
//    mat perm_mus=Permute_mat_sppmix(cur_mus,found);
//    mat sig=zeros(2,2),perm_sigmas=Permute_mat_sppmix(cur_sigmas,found);
    List mix2(m);
    for(j=0;j<m;j++)
    {
      sig(0,0)=perm_sigmas(j,0);
      sig(0,1)=perm_sigmas(j,1);
      sig(1,0)=perm_sigmas(j,2);
      sig(1,1)=perm_sigmas(j,3);
//     Rcout << perm_mus.row(j)<< std::endl ;
      mu1(0)=perm_mus(j,0);
      mu1(1)=perm_mus(j,1);
 //     Rcout << "passed2"<< std::endl ;
      mix2[j]=List::create(
        Named("p") = perm_ps(j),
        Named("mu") = mu1,
        Named("sigma") = sig);
      permgensigmas(i,j)=sig;
    }
    permuted_gens[i]=mix2;
    permgenps.row(i)=perm_ps.t();
    permgenmus.slice(i)=perm_mus;
  }
//apply the permutation to the realizations
//in the R code, just return the best one
//  printf("\rDone                                                      \n");
  return List::create(
//    Named("best_perm") = found,
    Named("best_perm") = current_perm,
    Named("permuted_gens") = permuted_gens,
    Named("permuted_ps") = permgenps,
    Named("permuted_mus") = permgenmus,
    Named("permuted_sigmas") = permgensigmas);
//  return found;
}


// [[Rcpp::export]]
List GetAllMeans_sppmix(List const& allgens,
                        int const& burnin)
{
//allgens is the list of lists of realizations
//  z=GetAllMeans_sppmix(gens$allgens,burnin)
  int i,j,L=allgens.size();
  List gen_realiz=allgens[0];
  int countgens=L-burnin,m=gen_realiz.size();
  vec sumps=zeros(m);
  mat summus=zeros(m,2);
  mat sumsigmas=zeros(m,4);

  for(i=burnin;i<L;i++)
  {
    sumps=sumps+GetRealiz_ps_sppmix(allgens,i);
    summus=summus+GetRealiz_mus_sppmix(allgens,i);
    sumsigmas=sumsigmas+GetRealiz_sigmas_sppmix(allgens,i);
//    countgens++;
  }
//  Rcout << countgens<< L-burnin<<std::endl ;
  mat sig(2,2);
  cube meansigmas(2,2,m);
  for(j=0;j<m;j++)
  {
    sig(0,0)=sumsigmas(j,0)/countgens;
    sig(0,1)=sumsigmas(j,1)/countgens;
    sig(1,0)=sumsigmas(j,2)/countgens;
    sig(1,1)=sumsigmas(j,3)/countgens;
    meansigmas.slice(j)=sig;
  }
  return List::create(
    Named("meanps") = sumps/countgens,
    Named("meanmus") = summus/countgens,
    Named("meansigmas") = meansigmas);
}


// [[Rcpp::export]]
vec GetCompDistr_sppmix(vec const& numcomp,
    int const& maxnumcomp)
{
  //returns the distribution of the # of components
  //apply burnin before calling this function
  int L=numcomp.size();
  vec distr_numcomp(maxnumcomp);
//  vec newcomps=SubVec_sppmix(numcomp,burnin,L);
  //   numcomp.subvec(burnin,L);
  for(int j=0;j<maxnumcomp;j++)
  {
    uvec q1 = find(numcomp==j);
    //    Rcout << q1 << std::endl ;
    distr_numcomp(j)=1.0*q1.size()/L;
  }
  return distr_numcomp;
}


// [[Rcpp::export]]
List GetBDCompRealiz_sppmix(List const& genBDmix,
    vec const& genlamdas,vec const& numcomp,
    int const& comp)
{
  //returns the gens for a specific # of components
  //apply burnin before calling this function
//  int L=numcomp.size();
  uvec indi=find(numcomp==comp);
  if(sum(indi)==0)//no realizations for this k
    return List::create();
  int nn1=indi.size();
  vec newlamdas=genlamdas(indi);
//  vec newnumcomp=numcomp(indi);
  List newgenmix(nn1);
  for(int i=0;i<nn1;i++)
    newgenmix[i]=genBDmix[indi(i)];

  return List::create(
    Named("newgenBD") = newgenmix,
    Named("newlamdas") = newlamdas);
    //,Named("newnumcomp") = newnumcomp);
}


// [[Rcpp::export]]
mat GetAvgLabelsDiscrete2Multinomial_sppmix(mat
      const& genzs,int const& m)
{
  //returns the membership matrix from
  //the realizations matrix Lxn (the output
  //from DAMCMC)
  //apply burnin before calling this function
  int i,dat,iters=genzs.n_rows,n=genzs.n_cols;
  mat zmultinomial=zeros(n,m);
//  vec ptavglabel=zeros(n);
  for(i=0;i<iters;i++)
  {
    //build the label matrix for this realization
    mat cur_z=zeros(n,m);
    //put 1 at location genzs(i,dat)
    for(dat=0;dat<n;dat++)
      cur_z(dat,genzs(i,dat))=1;
    zmultinomial+=cur_z;
//    ptavglabel+=genzs.row(i).t();
  }
  return zmultinomial/iters;
}


// [[Rcpp::export]]
bool Check4LabelSwitching_sppmix(vec const& chain,int const& lag)
{
  //checks for sharps changes in the generated
  //chains, works on a specific chain
  //similar problem to change-point detection
  //works using the past lag-realizations
  //apply burnin before calling this function
  int i,iters=chain.size();
//  int lag=floor(0.05*iters);//use 5% of iterations
  bool LabelSwitchingPresent=false;
  for(i=lag;i<iters-lag;i++)
  {
    vec lagv=chain.subvec(i-lag,i);
    double meanlag=mean(lagv),
      stdlag=sqrt(var(lagv));
    vec flagv=chain.subvec(i,i+lag);
    double fmeanlag=mean(flagv),
      fstdlag=sqrt(var(flagv));
    if((fmeanlag<meanlag-3*stdlag ||
       fmeanlag>meanlag+3*stdlag) &&
       (meanlag<fmeanlag-3*fstdlag ||
       meanlag>fmeanlag+3*fstdlag))
    {
      LabelSwitchingPresent=true;
//      Rcout << i << '\n'<< chain(i) << '\n'<<lag << '\n'<<meanlag << '\n'<<stdlag<< '\n'<<std::endl ;
      break;
    }
  }
  return LabelSwitchingPresent;
}


// [[Rcpp::export]]
List PostGenGetBestPermIdenConstraint_sppmix(List const& allgens)
{
  //allgens is all of the List output from DAMCMC
  //apply burnin before calling this function
  //it uses the identifiability constraint
  //to permute all realizations
  List allgens_List=allgens[0];
  int i,j,k,L=allgens_List.size();
//    Rcout << "L="<<L<<std::endl ;
  List permuted_gens(L),mix1=allgens_List[0];
  mat allgens_zs=allgens[4];
  int m=mix1.size(),n=allgens_zs.n_cols;
  cube permgenmus=zeros(m,2,L);
  field<mat> permgensigmas(L,m);
  mat permgenps=zeros(L,m);
  mat permgenzs=zeros(L,n),sig1(2,2);
//  Rcout << "m="<<m<<std::endl ;
//  Rcout << "L="<<L<<std::endl ;
  mat best_perm=zeros(L,m);//will contain the best perm
  vec cur_sig(4),mu1(2);
  mat allps=GetAllRealiz_ps_sppmix(allgens_List);
  List allmus=GetAllRealiz_mus_sppmix(allgens_List);
  List allsigmas=GetAllRealiz_sigmas_sppmix(allgens_List);
  for(i=0;i<L;i++)
  {
    //get the vector that orders the ps
    vec cur_ps=allps.row(i).t();//GetRealiz_ps_sppmix(allgens,i);
    mat cur_mus=allmus[i];//GetRealiz_mus_sppmix(allgens,i);
    mat cur_sigmas=allsigmas[i];//GetRealiz_sigmas_sppmix(allgens,i);
    //vec cur_ps=GetRealiz_ps_sppmix(allgens_List,i);
    //mat cur_mus=GetRealiz_mus_sppmix(allgens_List,i);
    //mat cur_sigmas=GetRealiz_sigmas_sppmix(allgens_List,i);
//    Rcout << "passed"<<std::endl ;
//    Rcout << "\ncur_ps="<< cur_ps.t()<<std::endl ;
    vec dists(m),permind(m);
    for(j=0;j<m;j++)
    {
      //find distance from origin and get ordering
      mu1=cur_mus.row(j).t();
      cur_sig=trans(cur_sigmas.row(j));
      sig1(0,0)=cur_sig(0);
      sig1(0,1)=cur_sig(1);
      sig1(1,0)=cur_sig(2);
      sig1(1,1)=cur_sig(3);
      dists(j)=//cur_ps(j)*exp(-.5*
     //   Quad_sppmix(mu1,invmat2d_sppmix(sig1));
        //)/sqrt(det(2*datum::pi*sig1));
      //cur_ps(j)*cur_ps(j)+
        sqrt(mu1(0)*mu1(0)+
          mu1(1)*mu1(1));
//        +cur_sig(0)*cur_sig(3)
 //       -cur_sig(1)*cur_sig(2);
/*        +cur_sig(0)*cur_sig(0)
        +cur_sig(1)*cur_sig(1)
        +cur_sig(2)*cur_sig(2)
        +cur_sig(3)*cur_sig(3);*/
    }
    uvec indices=sort_index(dists);// = sort_index(cur_ps);
    //    Rcout << "\nindices="<< indices.t()<<std::endl ;
    for(j=0;j<m;j++)
    {
      permind(j)=indices(j)+1;
    }
    best_perm.row(i)=permind.t();
 //   Rcout << "passed"<<std::endl ;
    vec perm_ps=Permute_vec_sppmix(cur_ps,permind);
//    Rcout << "\nperm_ps="<< perm_ps.t()<<std::endl ;
    //    Rcout << "passed1"<<std::endl ;
    mat perm_mus=Permute_mat_sppmix(cur_mus,permind);
 //   Rcout << "passed"<<std::endl ;
    mat sig=zeros(2,2),perm_sigmas=Permute_mat_sppmix(cur_sigmas,permind);
    List mix2(m);
    for(j=0;j<m;j++)
    {
      sig(0,0)=perm_sigmas(j,0);
      sig(0,1)=perm_sigmas(j,1);
      sig(1,0)=perm_sigmas(j,2);
      sig(1,1)=perm_sigmas(j,3);
      //     Rcout << perm_mus.row(j)<< std::endl ;
      mu1(0)=perm_mus(j,0);
      mu1(1)=perm_mus(j,1);
      //     Rcout << "passed2"<< std::endl ;
      mix2[j]=List::create(
        Named("p") = perm_ps(j),
        Named("mu") = mu1,
        Named("sigma") = sig);
      permgensigmas(i,j)=sig;
    }
    permuted_gens[i]=mix2;
    permgenps.row(i)=perm_ps.t();
    permgenmus.slice(i)=perm_mus;
    for(k=0;k<n;k++)
    {
 //genz=zeros(L,n)
 //zmultinomial(n,m)
 //build the label matrix for this realization
      vec permz(m),cur_z=zeros(m);
      //put 1 at location genzs(i,dat)
      cur_z(allgens_zs(i,k))=1;
      permz=Permute_vec_sppmix(cur_z,permind);
      for(j=0;j<m;j++)
        if(permz(j)==1)
        {
          permgenzs(i,k)=j;
          break;
        }
    }
  }
//  printf("\rDone applying identifiability constraint                      \n");
  return List::create(
    Named("allgens_List") = permuted_gens,
    Named("genps") = permgenps,
    Named("genmus") = permgenmus,
    Named("gensigmas") = permgensigmas,
    Named("genzs") = permgenzs,
    Named("genlamdas") = allgens[5],
    Named("best_perm") = best_perm);
}


// [[Rcpp::export]]
mat PermuteZs_sppmix(mat const& allgens_zs,
                     mat const& bestperm)
{
  //apply burnin before calling this function
//genz=zeros(L,n)
  int i,j,k,L=bestperm.n_rows,
     m=bestperm.n_cols;
  int n=allgens_zs.n_cols;
  mat permgenzs=zeros(L,n);
//    Rcout << "m="<<m<<std::endl ;
//    Rcout << "L="<<L<<std::endl ;
//    Rcout << "n="<<n<<std::endl ;
  vec cur_sig(4),mu1(2);
  for(i=0;i<L;i++)
  {
    for(k=0;k<n;k++)
    {
      vec permz(m),cur_z=zeros(m,1);
      //put 1 at location genzs(i,dat)
      cur_z(allgens_zs(i,k))=1;
      permz=Permute_vec_sppmix(cur_z,
          bestperm.row(i).t());
      for(j=0;j<m;j++)
        if(permz(j)==1)
        {
          permgenzs(i,k)=j;
          break;
        }
    }
  }
  return permgenzs;
}


// [[Rcpp::export]]
mat FisherInfoMat_sppmix(mat const& data,
  vec const& map_ps,mat const& map_mus,
  List const& map_sigmas,mat const& map_zs)
{
  //mus is a mx2 matrix
  //sigmas is a List with m elements
  //each being a 2x2 matrix
  //zs is an nxm matrix
  //apply burnin before calling this function
  int i,j,n=map_zs.n_rows,m=map_zs.n_cols;
  mat Fisher=zeros(6*m-1,6*m-1);
  cube sigma_inv(2,2,m);
  cube sigmas(2,2,m);
  //  Rcout << "m="<<m<<std::endl ;
  mat sig(2,2),invsig(2,2);
  List mixcomp;
  for(j=0;j<m;j++)
  {
    mixcomp=map_sigmas[j];
    sig = as<mat>(mixcomp["sigma"]);
    sigmas.slice(j)=sig;
//    Rcout << "map_sigmas="<<sig<<std::endl ;
    sigma_inv.slice(j)=invmat2d_sppmix(sig);
  }
//  Rcout << "n="<<n<<std::endl ;
  vec mu1(2),mu2(2);
  for(i=0;i<n;i++)
  {
    vec scorevec(6*m-1);
    int curdim=0;
    for(j=0;j<m-1;j++)
    {
      scorevec(curdim)=map_zs(i,j)/map_ps(j)-map_zs(i,m-1)/map_ps(m-1);
      curdim++;
    }
    for(j=0;j<m;j++)
    {
      mu1(0)=data(i,0)-map_mus(j,0);
      mu1(1)=data(i,1)-map_mus(j,1);
      invsig=sigma_inv.slice(j);
      mu2=-map_zs(i,j)*invsig*mu1;
      scorevec(curdim)=mu2(0);
      curdim++;
      scorevec(curdim)=mu2(1);
      curdim++;
    }
    for(j=0;j<m;j++)
    {
      mu1(0)=data(i,0)-map_mus(j,0);
      mu1(1)=data(i,1)-map_mus(j,1);
      sig=sigmas.slice(j);
      invsig=sigma_inv.slice(j);
      double qcol1=mu1(0)*sig(0,0)
        +mu1(1)*sig(1,0);
      double qcol2=mu1(0)*sig(0,1)
        +mu1(1)*sig(1,1);
      scorevec(curdim)=0.5*map_zs(i,j)*
        (-invsig(0,0)+qcol1*qcol1);
      curdim++;
      scorevec(curdim)=map_zs(i,j)*
        (-invsig(0,1)+qcol1*qcol2);
      curdim++;
      scorevec(curdim)=0.5*map_zs(i,j)*
        (-invsig(1,1)+qcol2*qcol2);
      curdim++;
    }
    Fisher+=scorevec*scorevec.t();
  }
  return Fisher;
}


// [[Rcpp::export]]
List GetDensityValues_sppmix(
  mat const& data,List const& fit,
  vec const& xlims, vec const& ylims)
{
  //allgens is all of the List output
  //(fit object) from DAMCMCExtras
  //apply burnin before calling this function
  //GetDensityValues_sppmix(as.matrix(genPPP$x,genPPP$y,genPPP$n,2),drop_realization(fit))
//  Rcout <<"start"<<std::endl ;
  int i,j,dat;
  List allgens_List=fit[0];
  int L=allgens_List.size();
  List mix1=allgens_List[0];
  mat allgens_zs=fit[4];//Lxn
  mat ApproxCompMass=fit[6];//Lxm
  int m=mix1.size(),n=allgens_zs.n_cols;
/*  double Rx=max(data.col(0))-min(data.col(0)),
    Ry=max(data.col(1))-min(data.col(1));
  vec ksi=zeros(2);
  ksi(0)=sum(data.col(0))/n;
  ksi(1)=sum(data.col(1))/n;
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
  mat sumxmu=zeros(2,2);
//  Rcout <<"sumxmu"<<std::endl ;
  for(int r=0;r<n;r++)
    sumxmu=sumxmu+trans(data.row(r)-ksi.t())*(data.row(r)-ksi.t());
  mat ps2=invmat2d_sppmix(2*Idenmat+sumxmu/(n-1));
*/
//  Rcout <<"data"<<data<<std::endl ;

  vec xy(2),mu1(2),genlambdas=fit[5],
     Density(L),logDensity(L);
  cube CompDensityAtXi(n,m,L);
//,CompDensityAtXiwrtPrior(n,m,L);
  mat DensityAtXi(n,L),
    //DensityAtXiwrtPrior(n,L),
    sig(2,2),invsig(2,2);
  double meanlamda=mean(genlambdas),
    constij,sumdev=0,lognfac=
    log(sqrt(2.0*datum::pi*n))+n*log(1.0*n)-n;
  //use stirling formula for
  //n!~=sqrt(2*datum::pi*n)*n^n*exp(-n)
  //for priors
//  vec ApproxCompwrtPrior(m),
//    prior_lambdas=rgamma(L,1,1),
//    prior_ps(m);
//  mat prior_mu(m,2),prior_sigma(m,4);
//  Rcout <<"start"<<std::endl ;
  mat allps=GetAllRealiz_ps_sppmix(allgens_List);
  List allmus=GetAllRealiz_mus_sppmix(allgens_List);
  List allsigmas=GetAllRealiz_sigmas_sppmix(allgens_List);
//  Rprintf("\nComputing density values...\n");
  for(i=0;i<L;i++)
  {
//    Rprintf("\rWorking: %3.1f%% complete",100.0*i/(L-1));
    //sample from prior
//    prior_ps=rDirichlet_sppmix(ones(m));
//    Rcout <<prior_ps<<std::endl ;
//    mat prior_mus=rnorm2_sppmix(m,ksi,kappainv);
 //   Rcout <<prior_mus<<std::endl ;
/*    cube gensigmas=zeros(2,2,m),
      geninvsigmas=zeros(2,2,m);
    for (j=0;j<m;j++)
    {
      geninvsigmas.slice(j)=rWishart_sppmix(2*3,ps2);//kappa;//invmat2d_sppmix(gensigmas(0,i));
      gensigmas.slice(j)=invmat2d_sppmix(geninvsigmas.slice(j));//kappainv;
      ApproxCompwrtPrior(j)=
        ApproxCompMass_sppmix(xlims,ylims,
          prior_mus.row(j).t(),gensigmas.slice(j));
    }*/
    //posterior samples
    vec cur_ps=allps.row(i).t();//GetRealiz_ps_sppmix(allgens,i);
    mat cur_mus=allmus[i];//GetRealiz_mus_sppmix(allgens,i);
    mat cur_sigmas=allsigmas[i];//GetRealiz_sigmas_sppmix(allgens,i);
//    vec cur_ps=GetRealiz_ps_sppmix(allgens_List,i);
//    mat cur_mus=GetRealiz_mus_sppmix(allgens_List,i);
//    mat cur_sigmas=GetRealiz_sigmas_sppmix(allgens_List,i);
//    mix1=allgens_List[i];
    logDensity(i)=0;
    for(dat=0;dat<n;dat++)
    {
      DensityAtXi(dat,i)=0;
      xy=trans(data.row(dat));
//      Rcout <<"\nxy\n"<< xy<<std::endl ;
      for(j=0;j<m;j++)
      {
        mu1(0)=cur_mus(j,0)-xy(0);
        mu1(1)=cur_mus(j,1)-xy(1);
//        Rcout <<"\nmu1\n"<< mu1<<std::endl ;
        sig(0,0)=cur_sigmas(j,0);
        sig(0,1)=cur_sigmas(j,1);
        sig(1,0)=cur_sigmas(j,2);
        sig(1,1)=cur_sigmas(j,3);
        invsig=invmat2d_sppmix(sig);
        if(ApproxCompMass(i,j)>0.00001)
        {
          constij=1.0/(ApproxCompMass(i,j)*
            sqrt(det(2*datum::pi*sig)));
          CompDensityAtXi(dat,j,i)=
            cur_ps(j)*constij*exp(-.5*
            Quad_sppmix(mu1,invsig));
        }
        else
        {
          CompDensityAtXi(dat,j,i)=0;
        }
//        Rcout <<"\nCompDensityAtXi(dat,j,i)\n"<< CompDensityAtXi(dat,j,i)<<std::endl ;
        DensityAtXi(dat,i)+=CompDensityAtXi(dat,j,i);
/*        //sample from the priors
        invsig=geninvsigmas.slice(j);
        sig=gensigmas.slice(j);
        constij=1.0/(ApproxCompwrtPrior(j)*
          sqrt(det(2*datum::pi*sig)));
        CompDensityAtXiwrtPrior(dat,j,i)=
          prior_ps(j)*constij*exp(-.5*
          Quad_sppmix(prior_mus.row(j).t(),
                      invsig));
        DensityAtXiwrtPrior(dat,i)+=CompDensityAtXiwrtPrior(dat,j,i);
*/
      }
      if(DensityAtXi(dat,i)>0.0001)
        logDensity(i)+=log(DensityAtXi(dat,i));
//      Rcout <<"\nDensityAtXi="<< DensityAtXi(dat,i)<<std::endl ;
    }
    Density(i)=exp(logDensity(i));
//    Rcout <<"\nDensityi="<< Density(i)<<std::endl ;
//    if(Density(i)>0.00001)
//      sumdev+=-2*log(Density(i));
    sumdev+=-2*logDensity(i);
//    Rcout <<"\nsumdev"<<sumdev<<"\n"<< std::endl ;
    Rcpp::checkUserInterrupt();
  }
//  double deviance_MCapprox=sumdev/L;
//  Rcout <<"\ndeviance_MCapprox"<<deviance_MCapprox<<"\n"<< std::endl ;
  //  Rcout <<"\nComputing the Entropy\n"<< std::endl ;
//  Rprintf("\rDone.                    ");

  double likMAP=0//,likAvg=0//,entropy=0,
  //  entropyavg=0
  ,entropyMAP=0;
//  Rprintf("\rComputing the entropy value...\n");
  for(dat=0;dat<n;dat++)
  {
//    Rprintf("\rWorking: %3.1f%% complete",100.0*dat/(n-1));
    double maxdens=-1;//max works better
    double avgdens=0;//average surface
//    int ind_maxdens=-1;
    for(i=0;i<L;i++)
    {
      if(DensityAtXi(dat,i)>maxdens)
      {
        maxdens=DensityAtXi(dat,i);
 //       ind_maxdens=i;
      }
      avgdens+=DensityAtXi(dat,i);
    }
    likMAP+=log(maxdens);
    //likAvg+=log(avgdens/L);
    //entropy
    double entropyMAPsum=0;//,entropysum=0,entropysum1=0;
    for(j=0;j<m;j++)
    {
      double zij,//minzij=2,
        MAPzij=-1;//,avgzij=0;
      for(i=0;i<L;i++)
      {
        zij=CompDensityAtXi(dat,j,i)/DensityAtXi(dat,i);
//        avgzij+=zij;
//        if(zij<minzij)        minzij=zij;
        if(zij>MAPzij)
          MAPzij=zij;
      }
      if(MAPzij>0)
        entropyMAPsum+=MAPzij*log(MAPzij);
  //    if(minzij>0)        entropysum+=minzij*log(minzij);
//      if(avgzij>0)    entropysum1+=avgzij*log(avgzij/L)/L;
    }
    entropyMAP+=-entropyMAPsum;
//    entropy+=-entropysum;
//    entropyavg+=-entropysum1;
  }
//  Rprintf("\rDone.                    ");

  double loglikelihoodMAP=n*log(meanlamda)-meanlamda-lognfac+likMAP;
//  double loglikelihoodAvg=n*log(meanlamda)-meanlamda-lognfac+likAvg;
  //  Rcout <<"\nEntropy aprox="<<entropy<< std::endl ;
//  Rcout <<"Neg LogLik approx="<<-loglikelihood<< std::endl ;

  //marginal calculations, choose any theta, say the MAP
  //  m(x)=f(x\theta)*pi(theta)/pi(theta\x)
  //log(f(x\theta)) is approximated by loglikelihood above
  //so that f(x\theta)=exp(loglikelihood)
  //pi(theta) is a constant, it does not depend on x
  //pi(theta\x) is approximated by
  //sum(over zij of pi(theta\x,zij))/L
  //use MAP estimators for theta
  double marginal=0;//,marginalwrtprior=0;
/*  for(i=0;i<L;i++)
  {
    double sumd=1;//,sumd1=0;
    for(dat=0;dat<n;dat++)
    {
      sumd+=log(DensityAtXi(dat,i));
//      sumd1+=log(DensityAtXiwrtPrior(dat,i));
    }
    marginal+=exp(sumd);
//    marginalwrtprior+=exp(sumd1);
  }
  marginal=marginal/L;*/
  marginal=mean(Density);
//  marginalwrtprior=marginalwrtprior/L;
  //  Rcout <<"\nMarginal Monte Carlo Approx="<<marginal<< std::endl ;
  //  marginal=exp(loglikelihood);
//  Rcout <<"Marginal aprox="<<marginal<< std::endl ;

  return List::create(
    Named("Marginal") = marginal,
    Named("LogLikelihood") = //loglikelihoodAvg,//
      loglikelihoodMAP,
    Named("CompDensityAtXi") = CompDensityAtXi,
    Named("DensityAtXi") = DensityAtXi,
    Named("EntropyMAP") = entropyMAP,
//    Named("DevianceMC") = deviance_MCapprox,
    Named("Density") = Density,
    Named("logDensity") =logDensity);
}


// [[Rcpp::export]]
double ComputeBayesFactor_sppmix(
  mat const& densvals1,mat const& densvals2)
{
  int n=densvals1.n_rows,
      L=densvals1.n_cols,i,j,r;
  double val=0,sum2=0;
  for(r=0;r<L;r++)
  {
    double sum1=0;
    for(j=0;j<L;j++)
    {
      double prod=1;
      for(i=0;i<n;i++)
      {
        if(densvals1(i,j)==0)
        {
          prod=0;
          break;
        }
        if(densvals1(i,r)==0)
        {
          prod=-1;
          break;
        }
        prod*=densvals1(i,j)/densvals2(i,r);
      }
      if(prod!=-1)
        sum1+=prod;
      else
      {
        sum1=0;
        break;
      }
    }
    if(sum1!=0)
     sum2+=1/sum1;
  }
  val=1/sum2;
  return val;
}
