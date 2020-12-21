#' Estimate a mixture model parameters using MCMC
#'
#' @description
#' These functions fit a Poisson point process with a mixture
#' intensity, in a Bayesian framework, using
#' the Data Augmentation MCMC (DAMCMC) method for a fixed number of components
#' or the Birth-Death MCMC (BDMCMC) method for a random,
#' finite number of components. Estimation of the parameters
#' of the models is accomplished via calculation of
#' the posterior means, based on posterior realizations obtained
#' by these computational methods.
#'
#' For DAMCMC examples see
#'
#' \url{http://faculty.missouri.edu/~micheasa/sppmix/sppmix_all_examples.html
#' #est_mix_damcmc}
#'
#' For BDMCMC examples see
#'
#' \url{http://faculty.missouri.edu/~micheasa/sppmix/sppmix_all_examples.html
#' #est_mix_bdmcmc}
#'
#' @param pp Point pattern object of class \code{\link[spatstat.geom]{ppp}}.
#' @param m Either the number of components to fit in DAMCMC or
#' the maximum number of components requested for a BDMCMC fit.
#' @param truncate Logical variable indicating whether or not we
#' normalize the densities of the mixture components
#' to have all their mass within the window defined
#' in the point pattern \code{pp}.
#' @param L Number of iterations for the DAMCMC (default is 10000) or BDMCMC (default is 30000). If the value passed is less than the default, then it is set to the default value.
#' @param hyper_da Hyperparameters for DAMCMC, default is (3, 1, 1).
#' @param useKmeans Logical variable. If TRUE use a kmeans clustering method to obtain the starting values for the component means, otherwise, randomly sample a point from the pattern and use it as a component mean.
#' @return An object of type \code{damcmc_res}, containing MCMC realizations.
#' \item{allgens_List}{A list of size \code{L} containing posterior realizations, with elements that are lists of size \code{m}, with each element being the posterior realization of the parameters of a component, i.e., \code{p}, \code{mu} and \code{sigma}.}
#' \item{genps}{An \code{Lxm} matrix containing the L posterior realizations of the m component probabilities of the normal mixture.}
#' \item{genmus}{An \code{mx2xL} array containing the L posterior realizations of m component mean vectors of the normal mixture.}
#' \item{gensigmas}{An \code{Lxm} list, with elements that are \code{2x2} matrices, the posterior realizations of the covariance matrices of the components of the normal mixture.}
#' \item{genzs}{An \code{Lxn} matrix containing the membership indicators of the \code{n} points of the point pattern \code{pp}, over all iterations of the MCMCM.}
#' \item{genlambdas}{An \code{Lx1} vector containing the posterior realizations of the of the lambda parameter (the average number of points over the window).}
#' \item{genlambdas}{An \code{Lx1} vector containing the posterior realizations of the of the lambda parameter (the average number of points over the window).}
#' \item{ApproxCompMass}{An \code{Lxm} matrix containing the approximate mass of the m mixture components for each iteration.}
#' \item{data}{The original point pattern.}
#' \item{L}{Same as input.}
#' \item{m}{Same as input.}
#' @author Jiaxun Chen, Sakis Micheas, Yuchen Wang
#' @details
#' The DAMCMC follows the sampling scheme of Diebolt and Robert (1994).
#' @references Diebolt, J., and Robert, C. P. (1994). Estimation of
#' Finite Mixture Distributions through Bayesian Sampling.
#' Journal of the Royal Statistical Society B, 56, 2, 363-375.
#' @rdname est_mix
#' @seealso \code{\link{PlotUSAStates}},
#' \code{\link{plotmix_2d}},
#' \code{\link{GetPMEst}},
#' \code{\link{plot_chains}},
#' \code{\link{check_labels}},
#' \code{\link{GetBDTable}},
#' \code{\link{GetBDCompfit}},
#' \code{\link{plot.intensity_surface}},
#' \code{\link{plot.bdmcmc_res}},
#' \code{\link{to_int_surf}},
#' \code{\link[spatstat.geom]{owin}},
#' \code{\link{plot_CompDist}},
#' \code{\link{drop_realization}},
#' \code{\link{plot_chains}},
#' \code{\link{plot_ind}},
#' \code{\link{FixLS_da}},
#' \code{\link{rnormmix}}
#' @examples
#' \donttest{
#' fit <- est_mix_damcmc(spatstat::redwood, m = 3)
#' fit
#' plot(fit)
#' #We work with the California Earthquake data. We fit an IPPP with intensity surface modeled
#' #by a mixture with 5 normal components.
#' CAfit=est_mix_damcmc(CAQuakes2014.RichterOver3.0, m=5, L = 20000)
#' #Now retrieve the surface of Maximum a Posteriori (MAP) estimates of the mixture parameter.
#' #Note that the resulting surface is not affected by label switching.
#' MAPsurf=GetMAPEst(CAfit)
#' #Plot the states and the earthquake locations along with the fitted MAP IPPP intensity
#' #surface.
#' ret=PlotUSAStates(states=c('California','Nevada','Arizona'), showcentroids=FALSE,
#'  shownames=TRUE,main= "Earthquakes in CA, 2014",pp=CAQuakes2014.RichterOver3.0, surf=MAPsurf,
#'  boundarycolor="white",namescolor="white")
#' plot(CAfit)
#' #check labels
#' check_labels(CAfit)
#' # Fix label switching, start with approx=TRUE
#' post_fixed = FixLS_da(CAfit, plot_result = TRUE)
#' plot_chains(post_fixed)
#' plot_chains(post_fixed,separate = FALSE)
#' #this one works generally better but it is slow for large m
#' post_fixed = FixLS_da(CAfit,approx=FALSE, plot_result = TRUE)
#' plot_chains(post_fixed)
#' plot_chains(post_fixed,separate = FALSE)}
#'
#' @export
est_mix_damcmc <- function(pp, m, truncate = FALSE,
                           L = 10000, hyper_da = c(3, 1, 1),
                           useKmeans=FALSE)
{
  #if(L<10000)  {    warning("Number of realizations is less than 10k. Setting it to 10k.");L=10000  }
  # Start the clock!
  ptm1 <- proc.time()
  fit <- DAMCMC2d_sppmix(points = cbind(pp$x, pp$y),
                         xlims = Window(pp)$xrange, ylims = Window(pp)$yrange,
                         m = m, truncate = truncate,
                         L = L, hyperparams = hyper_da,
                         useKmeans=useKmeans)
  # Stop the clock
  ptm<-proc.time() - ptm1
  cat(paste("\nComputation time in seconds:",ptm[[1]],"\n"))
  fit$data <- pp
  class(fit$data) <- c("sppmix", "ppp")
  fit$L <- L
  fit$m <- m
  class(fit) <- "damcmc_res"
  return(fit)
}


#' @export
#' @method print damcmc_res
print.damcmc_res <- function(x,...) {
  fit=x
  cat("Normal Mixture fit via Data Augmentation MCMC \n",
      "DAMCMC iterations:", fit$L, "\n",
      "Number of components:", fit$m, "\n")
}

#' @param lambda1 Parameter for the truncated Poisson prior; \code{lambda1} = 1 by default.
#' @param lambda2 Birth rate; \code{lambda2} = 10 by default.
#' @param hyper Hyperparameters for the hierarchical prior. See 'Details' for more information.
#' @rdname est_mix
#' @return Additional return values from the BDMCMC fit (this is a \code{bdmcmc_res} object).
#' \item{numcomp}{An \code{Lx1} vector containing the number of components the BDMCMC chooses to fit at each iteration.}
#' \item{maxnumcomp}{Same as input m.}
#' \item{Badgen}{An \code{Lx1} vector with 0-1 values, where 1 indicates that the realization should be dropped, or 0 if the realization should be kept. The BDMCMC can produce "bad" births and degenerate realizations (zero component mass) that end up in the chain for a specific number of components. Although these realizations do not affect averages of surfaces (no label switching problem and the corresponding component probability is zero), they have a detrimental effect when working for mixture deconvolution within the chain of a specific number of components, i.e., computing the posterior averages of the mixture parameters corresponding to the mixture with MAP number of components. If this parameter is 1 for some realizations, then dropping these degenerate realizations allows us to use the label switching algorithms efficiently and achieve mixture deconvolution. Obviously, the number of posterior realizations can drop significantly in number when we first apply burnin and then drop the bad realizations, so it is good practice to run the BDMCMC for at least 20000 iterations.}
#' @details
#' The BDMCMC uses the sampling scheme of Stephens (2000).
#' The definitions of the hyperparameters can be found
#' in equations (21)-(24). The number of components entertained
#' is from 1 to \code{m}, and the moves for the BDMCMC chain for the
#' number of components, are either one up (add a component) or one down
#' (remove a component). Make sure you plot the BDMCMC fit to obtain
#' additional information on the fit.
#' @references Stephens, M. (2000). Bayesian analysis of mixture models with an unknown
#' number of components-an alternative to reversible jump methods.
#' The Annals of Statistics, 28, 1, 40-74.
#' @examples
#' \donttest{
#' fitBD <- est_mix_bdmcmc(spatstat::redwood, m = 5)
#' fitBD
#' plotsBDredwood=plot(fitBD)
#' #Earthquakes example
#' CAfitBD=est_mix_bdmcmc(pp = CAQuakes2014.RichterOver3.0, m = 5)
#' BDtab=GetBDTable(CAfitBD)#retrieve frequency table and MAP estimate for number of components
#' BDtab
#' MAPm=BDtab$MAPcomp
#' plotsCAfitBD=plot(CAfitBD)
#' #get the surface of posterior means with MAP components and plot it
#' plotmix_2d(GetPMEst(CAfitBD,MAPm),CAQuakes2014.RichterOver3.0)
#' #retrieve all BDMCMC realizations corresponding to a mixture with MAP components
#' BDfitMAPcomp=GetBDCompfit(CAfitBD,MAPm)
#' BDfitMAPcomp
#' plot(BDfitMAPcomp$BDsurf,main=paste("Mixture intensity surface with",MAPm, "components"))
#' #Example of Dropping bad realizations and working with the MAP surface
#' open_new_plot=FALSE
#' truncate=FALSE
#' truemix4=rnormmix(m = 4, sig0 = .1, df = 5,xlim= c(-2,2), ylim = c(-2,2))
#' plot(truemix4,xlim= c(-2,2), ylim = c(-2,2),whichplots=0, open_new_window=
#'  open_new_plot)+add_title("True mixture of normals density")
#' trueintsurfmix4=to_int_surf(truemix4,lambda = 150,win =spatstat.geom::owin( c(-2,2),c(-2,2)))
#' #not truncating so let us use a larger window
#' bigwin=spatstat.geom::owin(c(-4,4),c(-4,4))
#' ppmix4 <- rsppmix(intsurf = trueintsurfmix4,truncate = truncate,win=bigwin)# draw points
#' print(plotmix_2d(trueintsurfmix4,ppmix4, open_new_window=open_new_plot,
#'  win=spatstat.geom::owin(c(-4,4),c(-4,4)))+add_title(
#'  "True Poisson intensity surface along with the point pattern, W=[-4,4]x[-4,4]",
#'  lambda =trueintsurfmix4$lambda,m=trueintsurfmix4$m,n=ppmix4$n))
#' BDMCMCfit=est_mix_bdmcmc(pp = ppmix4, m = 5,L=30000,truncate = truncate)
#' #check the original distribution of the number of components
#' plot_CompDist(BDMCMCfit, open_new_window=open_new_plot)
#' #get the realizations corresponding to the MAP number of components
#' BDtab=GetBDTable(BDMCMCfit,FALSE)#retrieve frequency table and MAP estimate for the
#' #number of components
#' MAPm=BDtab$MAPcomp
#' BDMCMCfitMAPcomp=GetBDCompfit(BDMCMCfit,MAPm)
#' BDMCMCfitMAPcompgens=BDMCMCfitMAPcomp$BDgens
#' #look at the range of the means with the degenerate realizations included
#' print(range(BDMCMCfitMAPcompgens$genmus[,1,]))
#' print(range(BDMCMCfitMAPcompgens$genmus[,2,]))
#' #use the original output of the BDMCMC and apply 10% burnin (default)
#' BDMCMCfit=drop_realization(BDMCMCfit)
#' #now we drop the bad realizations
#' BDMCMCfit=drop_realization(BDMCMCfit, (BDMCMCfit$Badgen==1))
#' #we see how many realizations are left
#' plot_CompDist(BDMCMCfit, open_new_window=open_new_plot)
#' #get the realizations for the MAP only
#' BDMCMCfitMAPcomp=GetBDCompfit(BDMCMCfit,MAPm)
#' BDMCMCfitMAPcompgens=BDMCMCfitMAPcomp$BDgens
#' #check again the range of values for the x-y coords of the component means; they
#' #should be within the window
#' print(range(BDMCMCfitMAPcompgens$genmus[,1,]))
#' print(range(BDMCMCfitMAPcompgens$genmus[,2,]))
#' #check the MAP surface
#' plotmix_2d(BDMCMCfitMAPcomp$BDsurf,ppmix4, open_new_window=open_new_plot,
#'  win=bigwin) +add_title("MAP Poisson intensity surface along with the point pattern",
#'  lambda =BDMCMCfitMAPcomp$BDsurf$lambda, m=BDMCMCfitMAPcomp$BDsurf$m, n=ppmix4$n,
#'  L=BDMCMCfitMAPcomp$BDgens$L)
#' plot_chains(BDMCMCfitMAPcompgens, open_new_window=open_new_plot, separate = FALSE)
#' #burnin has been applied, set to zero and check for label switching
#' labelswitch=check_labels(BDMCMCfitMAPcompgens,burnin=0)
#' #use the identifiability constraint approach first
#' post_fixedBDMCMCfitIC = FixLS_da(BDMCMCfitMAPcompgens,burnin=0)
#' plot_chains(post_fixedBDMCMCfitIC, open_new_window=open_new_plot, separate = FALSE)
#' print(plot_ind(post_fixedBDMCMCfitIC, burnin=0, open_new_window=
#'  open_new_plot)+add_title("Posterior means of the membership indicators (IC permuted labels)",
#'  m = post_fixedBDMCMCfitIC$m, n = post_fixedBDMCMCfitIC$data$n))
#' permSurfaceofPostMeansIC=GetPMEst(post_fixedBDMCMCfitIC, burnin=0)
#' print(plotmix_2d(permSurfaceofPostMeansIC,ppmix4, open_new_window=open_new_plot,
#'  win=bigwin)+add_title("Poisson surface of posterior means (IC)",
#'  lambda=permSurfaceofPostMeansIC$lambda, m=permSurfaceofPostMeansIC$m, n=ppmix4$n,
#'  L=post_fixedBDMCMCfitIC$L))
#' #use the decision theoretic approach via SEL to find the best permutation; this one should
#' #work much better
#' post_fixedBDMCMCfitSEL = FixLS_da(BDMCMCfitMAPcompgens,approx=FALSE, burnin=0)
#' plot_chains(post_fixedBDMCMCfitSEL, open_new_window=open_new_plot, separate = FALSE)
#' print(plot_ind(post_fixedBDMCMCfitSEL, burnin=0, open_new_window=open_new_plot)+add_title(
#'  "Posterior means of the membership indicators (best permutation)",
#'  m=post_fixedBDMCMCfitSEL$m,n = post_fixedBDMCMCfitSEL$data$n))
#' permSurfaceofPostMeansSEL=GetPMEst(post_fixedBDMCMCfitSEL, burnin=0)
#' print(plotmix_2d(permSurfaceofPostMeansSEL,ppmix4, open_new_window=open_new_plot,
#'  win=bigwin)+ add_title("Poisson surface of posterior means (best permutation)",
#'  lambda=permSurfaceofPostMeansSEL$lambda, m=permSurfaceofPostMeansSEL$m,
#'  n=ppmix4$n,L=post_fixedBDMCMCfitSEL$L))}
#'
#' @export
est_mix_bdmcmc <- function(pp, m, truncate = FALSE,
                           lambda1 = 1, lambda2 = 10, hyper = c(m/2,1/36,3,2,1,1),
                           L = 30000)
{
#  if(L<10){warning("Number of realizations is very small. Setting it to 1000.");L=1000  }
  # Start the clock!
  ptm1 <- proc.time()
  fit <- BDMCMC2d_sppmix(m, points = cbind(pp$x, pp$y),
                         xlims = Window(pp)$xrange, ylims = Window(pp)$yrange,
                         truncate = truncate,
                         lamda = lambda1, lamdab = lambda2, hyper = hyper,
                         L = L)
  ptm<-proc.time() - ptm1
  cat(paste("\nComputation time in seconds:",ptm[[1]],"\n"))
  fit$data <- pp
  class(fit$data) <- c("sppmix", "ppp")
  fit$L <- L
  class(fit) <- "bdmcmc_res"
  return(fit)
}


#' @export
#' @method print bdmcmc_res
print.bdmcmc_res <- function(x,...) {
  fit=x
  cat("Normal Mixture fit via Birth-Death MCMC\n",
      "BDMCMC iterations:", fit$L, "\n",
      "Number of components:")
  print(table(fit$numcomp))
}

#' Checking convergence: diagnostics
#'
#' @description
#' This function reports the Gelman-Rubin
#' convergence diagnostic \code{R} (also known as
#' the potential scale reduction), by producing
#' \code{k} DAMCMC fits and computing the
#' within-chain and between-chain variances. Values approximately
#' equal to 1 indicate convergence, otherwise we need to run the chain for
#' a longer number of iterations to get convergence.
#'
#' For DAMCMC examples see
#'
#' \url{http://faculty.missouri.edu/~micheasa/sppmix/sppmix_all_examples.html
#' #Get_Rdiag}
#'
#' @param pp Point pattern object of class \code{\link[spatstat.geom]{ppp}}.
#' @param m The number of components to fit.
#' @param truncate Logical variable indicating whether or not we
#' normalize the densities of the mixture components
#' to have all their mass within the window defined
#' in the point pattern \code{pp}.
#' @param L Number of iterations to use for
#' each chain created; default is 20000. Note
#' that half of them will be dropped so use a large number.
#' @param numofchains Number of chains to create; default is 2.
#' @param permute Request to generate chains that are
#' unpermuted (\code{permute=0}), identifiability
#' constraint (IC) permuted (\code{permute=1}), or minimum squared
#' error loss (SEL) permuted (\code{permute=2}).
#'
#' @author Sakis Micheas
#' @seealso \code{\link{est_mix_damcmc}},
#' \code{\link{rmixsurf}},
#' \code{\link{rsppmix}}
#' @examples
#' \donttest{
#' truemix_surf <- rmixsurf(m = 3, lambda=100, xlim = c(-3,3), ylim = c(-3,3))
#' plot(truemix_surf)
#' genPPP=rsppmix(intsurf = truemix_surf, truncate = FALSE)
#' Get_Rdiag(pp = genPPP, m = 3)}
#'
#' @export
Get_Rdiag<-function(pp, m, truncate = FALSE,
    L = 20000,numofchains=2,permute=0)
{
  L1=L
  truncate1=truncate
  AllChains_lambdas=matrix(0,numofchains,L)
  AllChains_ps=array(0,c(numofchains,L,m))
  AllChains_musx=array(0,c(numofchains,L,m))
  AllChains_musy=array(0,c(numofchains,L,m))
  AllChains_sig11=array(0,c(numofchains,L,m))
  AllChains_sig12=array(0,c(numofchains,L,m))
  AllChains_sig22=array(0,c(numofchains,L,m))
  for(k in 1:numofchains)
  {
    cat("\nProcessing chain",k,"\n")
    fit <- est_mix_damcmc(pp = pp, m = m,
          truncate = truncate1, L = L1)
    if(permute==1)
      fit =FixLS_da(fit, burnin=0,approx=TRUE,run_silent=TRUE)
    if(permute==2)
      fit =FixLS_da(fit, burnin=0, approx=FALSE,run_silent=TRUE)
    AllChains_lambdas[k,]=t(fit$genlamdas)
    for(i in 1:L)
      for(j in 1:m)
      {
        AllChains_ps[k,i,j]=fit$genps[i,j]
        AllChains_musx[k,i,j]=fit$genmus[j,1,i]
        AllChains_musy[k,i,j]=fit$genmus[j,2,i]
        AllChains_sig11[k,i,j]=fit$gensigmas[i,j][[1]][1,1]
        AllChains_sig12[k,i,j]=fit$gensigmas[i,j][[1]][1,2]
        AllChains_sig22[k,i,j]=fit$gensigmas[i,j][[1]][2,2]
      }
  }
  cat("\nDone. Calculating R-diagnostics...")
  Rps=rep(0,m)
  Rmusx=rep(0,m)
  Rmusy=rep(0,m)
  Rsig11=rep(0,m)
  Rsig12=rep(0,m)
  Rsig22=rep(0,m)
  for(j in 1:m)
  {
    Rps[j]=Get_Rdiag_single(AllChains_ps[,,j])
    Rmusx[j]=Get_Rdiag_single(AllChains_musx[,,j])
    Rmusy[j]=Get_Rdiag_single(AllChains_musy[,,j])
    Rsig11[j]=Get_Rdiag_single(AllChains_sig11[,,j])
    Rsig12[j]=Get_Rdiag_single(AllChains_sig12[,,j])
    Rsig22[j]=Get_Rdiag_single(AllChains_sig22[,,j])
  }
  Rlambda=Get_Rdiag_single(AllChains_lambdas)
  cat(" Done\n")
  return(list(Rlambda=Rlambda,
              Rps=Rps,
              Rmusx=Rmusx,
              Rmusy=Rmusy,
              Rsig11=Rsig11,
              Rsig12=Rsig12,
              Rsig22=Rsig22))
}

Get_Rdiag_single<-function(chains)
{
  #chains is numofchainsxL
  numofchains=nrow(chains)
  L=ncol(chains)
  chains_means=rep(0,numofchains)
  chains_vars=rep(0,numofchains)
  for(k in 1:numofchains)
  {
    chains_means[k]=mean(chains[k,])
    chains_vars[k]=(L-1)*var(chains[k,])
  }
  grandmean=mean(chains_means)
  B=L*sum((chains_means-grandmean)^2)/(numofchains-1)
  W=mean(chains_vars)
  Varhat_plus=(L-1)*W/L+B/L
  R=sqrt(Varhat_plus/W)
  return(R)
}

#' @export
est_mix_RMCPdamcmc <- function(
  pp, m, truncate = FALSE,L = 10000,
  d,mu0,Sigma0,df0,sig0,
  useKmeans=TRUE,startmus)
{
#  if(L<10000){warning("Number of realizations is less than 10k. Setting it to 10k.")L=10000}
  # Start the clock!
  ptm1 <- proc.time()
  if(missing(Sigma0))
    Sigma0=stats::cov(cbind(pp$x,pp$y))
  if(missing(mu0))
    mu0=c(mean(pp$x),mean(pp$y))
  if(missing(sig0))
    sig0=1
  if(missing(df0))
    df0=10
  if(missing(d))
    d=rep(1,m)
  if(length(d)!=m)
    stop("d is not the same length as m")
  if(!useKmeans)
  {
    if(missing(startmus))
    {
      startmus=matrix(0,m,2)
      for(i in 1:m)
      {
        ind=sample(x=1:pp$n,size=1,replace=TRUE)
        startmus[i,]=c(pp$x[ind],pp$y[ind])
      }
    }
  }
  fit <- DAMCMC2dRMCP_sppmix(
    points = cbind(pp$x, pp$y),
    xlims = Window(pp)$xrange,
    ylims = Window(pp)$yrange,
    m = m, L = L, truncate = truncate,
    d=d,mu0=mu0,Sigma0=Sigma0,
    df0=df0,sig0=sig0,
    useKmeans=useKmeans,
    startmus=startmus)
  # Stop the clock
  ptm<-proc.time() - ptm1
  cat(paste("\nComputation time in seconds:",ptm[[1]],"\n"))
  fit$data <- pp
  class(fit$data) <- c("sppmix", "ppp")
  fit$L <- L
  fit$m <- m
  class(fit) <- "damcmc_res"
  return(fit)
}
