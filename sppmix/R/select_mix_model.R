#' Mixture Model Selection
#'
#' @description
#' This function suggests the best number of components by computing
#' model selection criteria, including
#' AIC (Akaike Information Criterion),
#' BIC (Bayesian Information Criterion),
#' ICLC (Integrated Classification Likelihood Criterion).
#'
#' Since the only parameter of interest is the number
#' of components of the mixture, we consider several fixed
#' numbers of components defined in the vector Ms,
#' and we entertain mixture models with their other parameters
#' approximated via the MAP estimators of DAMCMC runs.
#'
#' For examples see
#'
#' \url{http://faculty.missouri.edu/~micheasa/sppmix/sppmix_all_examples.html
#' #selectMix}
#'
#' @param pp Point pattern object of class \code{ppp}.
#' @param Ms A vector of integers, representing different numbers of components
#' to assess for the mixture model for the intensity function.
#' @param L Number of iterations for the DAMCMC we run for each number of components in \code{Ms}; default is 30000.
#' @param burnin Number of initial realizations to discard. By default, it is 1/10 of the total number of iterations.
#' @param truncate Logical variable indicating whether or not we
#' normalize the densities of the mixture components
#' to have all their mass within the window defined
#' in the point pattern \code{pp}.
#' @param runallperms Set to 0 to use an approximation to the
#' Likelihood and Entropy within the MCMC (not affected by label switching).
#' Set to 1 to use an identifiability constraint to permute the labels and use the posterior means of
#' the parameters to compute the criteria. Set to 2 to use the
#' decision theoretic approach (minimize Squared Error Loss) in order to permute the labels.
#' The latter setting can take a long time to run for m>7.
#' @return A list containing the following components:
#' \item{AIC}{the values of the AIC criterion}
#' \item{BIC}{the values of the BIC criterion}
#' \item{ICLC}{the values of the ICLC criterion}
#' \item{Marginal}{the values of the marginal density}
#' \item{LogLikelihood}{the values of the LogLikelihood}
#'
#' @details For each integer in the vector \code{Ms}, we
#' fit a mixture with that many components
#' using DAMCMC. Then the criteria are computed and presented at the end of the
#' calculations.
#'
#' Note that the AIC and BIC do not account for constraints in the parameter
#' space of the mixture model parameters. The ICLC uses the estimated
#' entropy of the distribution of the membership indicators and therefore
#' should be trusted more in identifying the true number of components, instead of AIC and BIC.
#'
#' In addition, we run Stephens' BDMCMC and present the posterior
#' distribution for the number of components.
#'
#' All these methods should serve us in making
#' an informed choice about the true number of components
#' and then proceed to fit the DAMCMC with the
#' chosen number of components and take care of
#' label switching (if present), in order to
#' achieve mixture deconvolution. If we simply
#' want the surface, then the Bayesian model
#' average from the BDMCMC fit is the best solution.
#'
#' @references
#' Stephens, M. (2000). Bayesian analysis of mixture models with an unknown number of components: an alternative to reversible jump methods. The Annals of Statistics, 28, 1, 40-74.
#'
#' McLachlan, G., and Peel, D. (2000). Finite Mixture Models. Wiley-Interscience.
#'
#' Jasra, A., Holmes, C.C. and Stephens, D. A. (2005). Markov Chain Monte Carlo Methods and the Label Switching Problem in Bayesian Mixture. Statistical Science, 20, 50-67.
#' @author Jiaxun Chen, Sakis Micheas
#' @seealso \code{\link{normmix}},
#' \code{\link[spatstat]{square}},
#' \code{\link{est_mix_damcmc}},
#' \code{\link{est_mix_bdmcmc}},
#' \code{\link{GetBMA}},
#' \code{\link{FixLS_da}},
#' \code{\link{rsppmix}}
#' @examples
#' \donttest{
#' # create the true mixture intensity surface
#' truesurf <- normmix(ps=c(.2, .6,.2), mus=list(c(0.3, 0.3), c(0.7, 0.7), c(0.5, 0.5)),
#'  sigmas = list(.01*diag(2), .01*diag(2), .01*diag(2)), lambda=100, win=spatstat::square(1))
#' plot(truesurf)
#' # generate the point pattern, truncate=TRUE by default
#' pp <- rsppmix(truesurf,truncate=FALSE)
#' plot(pp,mus=truesurf$mus)
#' # compute model selection criteria via an approximation that is not affected by label
#' # switching and will typically work well for large L
#' ModelSel=selectMix(pp,1:5,truncate=FALSE)
#' # show info
#' ModelSel
#' #generate the intensity surface randomly
#' truesurf <- rmixsurf(5,100,xlim = c(-3,3), ylim = c(-3,3), rand_m = TRUE)
#' truesurf
#' pp <- rsppmix(truesurf,truncate=FALSE)
#' ModelSel0=selectMix(pp,1:5,runallperms = 0, truncate=FALSE)
#' ModelSel1=selectMix(pp,1:5,runallperms = 1, truncate=FALSE)
#' ModelSel2=selectMix(pp,1:5,runallperms = 2, truncate=FALSE)}
#'
#' @export
selectMix <- function(pp, Ms, L = 30000, burnin = .1*L,
    truncate = FALSE,runallperms=0)
{
  if (any(Ms < 1)) {
    stop("Operation exited. Must have at least 1 component.")
  }
  n <- pp$n
  pattern <- cbind(pp$x,pp$y)
  lognfac <- log(sqrt(2*pi*n)) + n*log(n) - n
  mmax <- length(Ms)
  AIC <- rep(0, mmax)
  BIC <- rep(0, mmax)
  ICLC <- rep(0, mmax)
  entropy <- rep(0, mmax)
  marginal <- rep(0, mmax)
  loglikelihood <- rep(0, mmax)
  # Start the clock!
  ptmtot <- proc.time()
  densvals=vector("list",mmax)
  for (m in 1:mmax)
  {
    cat(paste("\n================ # Components=",Ms[m],"===============\n"))
    post_real1<-est_mix_damcmc(pp, m = Ms[m], L = L,truncate = truncate)
    post_real=drop_realization(post_real1,drop=burnin)
    meanlamda=mean(post_real$genlamdas)
    if(runallperms==0)
    {
      cat("\nCalculating the entropy and the marginal...")
      vals=GetDensityValues_sppmix(pattern,
                                   post_real,as.vector(pp$window$xrange),
                                   as.vector(pp$window$yrange))
      densvals[[m]]=vals$DensityAtXi
      marginal[m]=vals$Marginal
      entropy[m]=vals$EntropyMAP
      loglikelihood[m]=vals$LogLikelihood
      deviance_atmean=-2*loglikelihood[m]
    }
    else if(runallperms<=2)
    {
      cat("\nFixing labels, calculating the entropy and the marginal...")
      fixed=post_real
      if(Ms[m]>1)
        fixed<-FixLS_da(post_real,
                      burnin=burnin,
                      approx = (runallperms==1),
                      plot_result = FALSE,
                      run_silent=FALSE)

      vals=GetDensityValues_sppmix(pattern,
                                   fixed,as.vector(pp$window$xrange),
                                   as.vector(pp$window$yrange))
      densvals[[m]]=vals$DensityAtXi
      marginal[m]=vals$Marginal
      mix_of_postmeans<-MakeMixtureList(fixed$allgens_List,burnin=0)
      fixednormmix=MakeNormMixFromMixtureList(mix_of_postmeans)
      approxcomp=rep(1,Ms[m])
      if(truncate)
        approxcomp=approx_normmix(fixednormmix,xlim= as.vector(pp$window$xrange),
                                ylim = as.vector(pp$window$yrange))
      den <-densNormMix_atxy_sppmix(pattern, mix_of_postmeans,approxcomp)
      loglikelihood[m] <- -lognfac + n*log(meanlamda) - meanlamda +
        sum(log(den))
      deviance_atmean=-2*loglikelihood[m]

      entropy[m] <-0
      if(Ms[m]>1)
      {
        zs <- GetAvgLabelsDiscrete2Multinomial_sppmix(
        fixed$genzs, Ms[m])
        zsn0 <- zs[zs!=0]
        entropy[m] <- sum(-zsn0*log(zsn0))
      }
      cat(" Done\n")
    }
    else
    {
      stop("Bad choice for parameter runallperms")
    }
    r <-6 * Ms[m]
    #standard deviance=-2*loglikelihood at the means
    AIC[m] <- 2*r+deviance_atmean
    BIC[m] <- r*log(n)+deviance_atmean
    ICLC[m] <- r*log(n)+deviance_atmean+2*entropy[m]
    cat(paste("\nAIC=",AIC[m]))
    cat(paste("\nBIC=",BIC[m]))
    cat(paste("\nICLC=",ICLC[m]))
  }

  cat("\n===============================\n")
  RVAL <- list(AIC = AIC,
               BIC = BIC,
               ICLC = ICLC,
               Marginal=marginal,
               LogLikelihood=loglikelihood)
  #run BDMCMC
  BDMCMCfit=est_mix_bdmcmc(pp, m = max(Ms), L = L,truncate = truncate)
  BDMCMCfit=drop_realization(BDMCMCfit)
  BDMCMCfit=drop_realization(BDMCMCfit, (BDMCMCfit$Badgen==1))
  plot_CompDist(BDMCMCfit,FALSE)
  BDtab=GetBDTable(BDMCMCfit,FALSE)#retrieve frequency table and MAP estimate for number of components
  MAPm=BDtab$MAPcomp
  PostMeanm=BDtab$MeanComp
  # Stop the clock
  ptmtotal<-proc.time() - ptmtot
  cat(paste("\nTotal computation time in seconds:",ptmtotal[[1]]))

  cat(paste("\nBDMCMC suggests (MAP estimator)",MAPm,"components"))
  cat(paste("\nBDMCMC suggests (posterior mean estimator)",PostMeanm,"components"))
  cat(paste("\nAIC minimum for",Ms[which.min(AIC)],"components"))
  cat(paste("\nBIC minimum for",Ms[which.min(BIC)],"components"))
  cat(paste("\nICLC minimum for",Ms[which.min(ICLC)],"components"))
  return(RVAL)
}
