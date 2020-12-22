#' Retrieve the Surface of Posterior Means
#'
#' @description
#' The function first calculates the posterior means
#' of the parameters of the components of the
#' mixture intensity, based on a DAMCMC or
#' BDMCMC fit. Then the surface of posterior
#' means is calculated using the posterior
#' means of the parameters. For a BDMCMC fit,
#' the number of components should be
#' specified, and all realizations with
#' that number of components are gathered
#' to calculate the posterior intensity surface.
#'
#' For examples see
#'
#' \url{http://faculty.missouri.edu/~micheasa/sppmix/sppmix_all_examples.html
#' #GetPMEst}
#'
#' @param fit Object of class \code{damcmc_res} or \code{bdmcmc_res}.
#' @param burnin Number of initial realizations to discard. By default, it is 1/10 of the total number of iterations.
#' @param num_comp Number of components requested. The posterior will be calculated only
#' based on the posterior realizations that have this many mixture components. If missing the realizations
#' corresponding to the MAP number of components are returned. This parameter
#' is ignored if \code{fit} is of class \code{damcmc_res}.
#'
#' @return An object of class \code{intensity_surface}.
#' @author Jiaxun Chen, Sakis Micheas, Yuchen Wang
#' @seealso \code{\link{est_mix_damcmc}},\code{\link{est_mix_bdmcmc}}
#' @examples
#' \donttest{
#' fit <- est_mix_damcmc(pp = spatstat.data::redwood, m = 3)
#' post_intsurf <- GetPMEst(fit, burnin = 1000)
#' plot(post_intsurf)
#' fit <- est_mix_bdmcmc(pp = spatstat.data::redwood, m = 5)
#' post_intsurf <- GetPMEst(fit, num_comp = 4, burnin = 1000)
#' plot(post_intsurf)
#' post_fixed = FixLS_da(fit,approx=FALSE, plot_result = TRUE)
#' plot(GetPMEst(post_fixed))}
#'
#' @export
GetPMEst <- function(
  fit, num_comp=1, burnin = floor(fit$L / 10))
{
  fit_burnined <- drop_realization(fit, burnin)
  if (class(fit) == "bdmcmc_res")
  {
    if(missing(num_comp))
    {
      tab=GetBDTable(fit_burnined,FALSE)
      num_comp=tab$MAPcomp
    }

    if (!(num_comp %in% unique(fit_burnined$numcomp))) {
      stop(paste0("This BDMCMC chain does not contain any ", num_comp,
                "-component realizations after burn-in."))
    }
    fit_burnined <- drop_realization(fit_burnined, fit_burnined$numcomp != num_comp)
    fit_burnined$m <- num_comp
  }
  post_ps <- colMeans(fit_burnined$genps[, 1:fit_burnined$m, drop = FALSE])
  mus <- apply(fit_burnined$genmus[1:fit_burnined$m, , ,drop = FALSE], 1:2, mean)

  sigmas <- apply(fit_burnined$gensigmas[, 1:fit_burnined$m, drop = FALSE], 2,
                    function(mats) Reduce(`+`, mats) / length(mats))

  mean_lambda <- mean(fit_burnined$genlamdas)

  post_mus <- post_sigmas <- vector("list", fit_burnined$m)
  for (i in seq_along(post_mus))
  {
    post_mus[[i]] <- mus[i, ]
    post_sigmas[[i]] <- matrix(sigmas[, i], 2, 2)
  }

  normmix(post_ps, post_mus, post_sigmas, mean_lambda,
            fit$data$window, estimated = TRUE)
}

#' Drop MCMC realizations
#'
#' @description
#' The function drops realizations from a DAMCMC or BDMCMC fit and
#' returns the resulting fit object.
#'
#' For examples see
#'
#' \url{http://faculty.missouri.edu/~micheasa/sppmix/sppmix_all_examples.html
#' #drop_realization}
#'
#' @param fit Object of class \code{damcmc_res} or \code{bdmcmc_res}.
#' @param drop If one integer is provided, the function will drop the first 1:drop
#' realizations. If an integer vector is provided, it will drop these
#' iterations. If a logical vector is provided (with the same length as the chain
#' length of \code{fit}), it will be used for subsetting directly.
#'
#' @author Sakis Micheas, Yuchen Wang
#' @seealso \code{\link{est_mix_bdmcmc}}
#' @examples
#' \donttest{
#' fit <- est_mix_bdmcmc(spatstat.data::redwood, m = 5)
#' fit
#' drop_realization(fit, 500)
#' drop_realization(fit, fit$numcomp != 5)}
#'
#' @export
drop_realization <- function(fit, drop=.1*fit$L) {
  # Generalized function for burn-in or drop realizations
  if (is.numeric(drop)) {
    if (length(drop) == 1) keep <- seq_len(fit$L) > drop else keep <- -drop
  }

  if (is.logical(drop)) keep <- !drop

  newfit<-fit
  newfit$allgens_List <- fit$allgens_List[keep]
  newfit$genps <- fit$genps[keep, , drop = FALSE]
  newfit$genmus <- fit$genmus[, , keep, drop = FALSE]
  newfit$gensigmas <- fit$gensigmas[keep, , drop = FALSE]
  newfit$genzs <- fit$genzs[keep, , drop = FALSE]
  newfit$genlamdas <- fit$genlamdas[keep, , drop = FALSE]
  newfit$ApproxCompMass <- fit$ApproxCompMass[keep, , drop = FALSE]

  newfit$L <- length(newfit$allgens_List)
  if (class(fit) == "bdmcmc_res") {
    newfit$numcomp <- fit$numcomp[keep, , drop = FALSE]
    newfit$Badgen <- fit$Badgen[keep, , drop = FALSE]
  }

  return(newfit)
}

#' Retrieve the surface of MAP estimators
#'
#' @description
#' The function calculates the Maximum A Posteriori (MAP)
#' estimate of the IPPP mixture intensity surface parameters. Use function
#' \code{\link{GetPMEst}} if you want the surface of posterior means.
#'
#' For examples see
#'
#' \url{http://faculty.missouri.edu/~micheasa/sppmix/sppmix_all_examples.html
#' #GetMAPEst}
#'
#' @param fit Object of class \code{damcmc_res} or \code{bdmcmc_res}.
#' @param burnin Number of initial realizations to discard. By default, it is 1/10 of the total number of iterations.
#' @param vals Contains the density values over the point pattern and realizations in the \code{fit} object. This
#' can be obtained via a call to \code{\link{GetDensityValues}}. If this argument is missing
#' then the density values are computed herein before computing the MAP estimates.
#' @param truncate Logical variable indicating whether or not we
#' normalize the densities of the mixture components
#' to have all their mass within the window defined
#' in the point pattern \code{pp}.
#' @param priortype For different types of priors, ignored right now.
#' @param d,mu0,Sigma0,df0,sig0 Optional parameters for the prior distributions used: d are the weights
#' of the Dirichlet prior on the component probabilities.
#' mu0 and Sigma0 are the mean and covariance matrix of a bivariate normal that yields the component means.
#' df0 and sig0 are the degrees of freedom and sig0^2*Identity the parameter matrix for the Inverse Wishart prior that yields the component matrices.
#' If omitted they are set to the following values, which are the default values used in est_mix_damcmc:
#' Sigma0=cov(cbind(pp$x,pp$y)), mu0=c(mean(pp$x),mean(pp$y)),
#' sig0=1, df0=10, and d=rep(1,m).
#' @return An object of type \code{intensity_surface}.
#' @seealso \code{\link{est_mix_damcmc}},
#' \code{\link{rmixsurf}},
#' \code{\link{rsppmix}},
#' \code{\link{GetPMEst}}
#' @author Sakis Micheas
#' @examples
#' \donttest{
#' truemix_surf <- rmixsurf(m = 3, lambda=100, xlim = c(-3,3), ylim = c(-3,3))
#' plot(truemix_surf,main="True IPPP intensity surface")
#' genPPP=rsppmix(intsurf = truemix_surf, truncate = FALSE)
#' #the larger the number of realizations the better
#' fit <- est_mix_damcmc(genPPP, m = 3,L=100000)
#' MAPest=GetMAPEst(fit)
#' plot(GetPMEst(fit),main="IPPP intensity surface of posterior means")
#' plot(MAPest,main="IPPP intensity surface of MAP estimates")
#' fitBD <- est_mix_bdmcmc(pp = genPPP, m = 5)
#' MAPest=GetMAPEst(fitBD)
#' plot(MAPest,main="IPPP intensity surface of MAP estimates for MAP m")}
#'
#' @export
GetMAPEst <- function(fit,
    burnin = floor(fit$L / 10),vals,
    truncate=FALSE,priortype=1,d,mu0,Sigma0,df0,sig0)
{
  cat("\nWorking...")
  fit_burnined <- drop_realization(fit, burnin)
  if(class(fit)=="bdmcmc_res")
  {
    tab=GetBDTable(fit_burnined,FALSE)
    fit_burnined=GetBDCompfit(fit_burnined,tab$MAPcomp,burnin=0)$BDgens
  }

  m<-fit_burnined$m
  pp=fit_burnined$data
  if(truncate)
  {
    pp=spatstat.geom::ppp(pp$x,pp$y,window=pp$window)
  }
  n <- pp$n
  pattern <- cbind(pp$x,pp$y)
  if(missing(Sigma0))
    Sigma0=stats::cov(cbind(pp$x,pp$y))
  if(missing(mu0))
    mu0=c(mean(pp$x),mean(pp$y))
  if(missing(sig0))
    sig0=1
  if(missing(df0))
    df0=10
  if(missing(d))
    d=rep(1, m)
  else
    d=suppressWarnings(as.vector(matrix(data=d,nrow=m,ncol=1)))
  if(missing(vals))
  {
    vals=GetDensityValues_sppmix(pattern,
          fit_burnined,as.vector(pp$window$xrange),
          as.vector(pp$window$yrange))
    priorvals=GetPriorVals(
      pp=pattern,allgens=fit_burnined$allgens_List,priortype,d,mu0,Sigma0,df0,sig0)
  }
  else
  {
    priorvals=GetPriorVals(
      pp=pattern,allgens=fit_burnined$allgens_List,priortype,d,mu0,Sigma0,df0,sig0)
    if(length(vals$Density)!=length(priorvals))
      stop("Make sure you apply the right burnin or don't pass the density values so they're calculated again.")
  }
  logdensityvals=vals$logDensity+log(priorvals)
  MAPiter=which.max(logdensityvals)
  MAPlambda=n
  MAPps <- fit_burnined$genps[MAPiter,]
  mus <- matrix(fit_burnined$genmus[, ,
              MAPiter],m,2)
  sigmas <- fit_burnined$gensigmas[MAPiter,]

  MAPmus <- MAPsigmas <- vector("list", fit_burnined$m)
  for (i in seq_along(MAPmus)) {
    MAPmus[[i]] <- mus[i, ]
    MAPsigmas[[i]] <- matrix(unlist(sigmas[i]), 2, 2)
  }
  cat("Done\n")
  normmix(MAPps, MAPmus, MAPsigmas, MAPlambda,
          fit_burnined$data$window, estimated = TRUE)
}

#' Retrieve density values
#'
#' @description
#' This function operates on the point pattern and the realizations of
#' a DAMCMC or BDMCMC fit (object \code{damcmc_res} or \code{bdmcmc_res})
#' and returns a plethora of information about the fit. When a \code{bdmcmc_res} is passed,
#' only the realizations corresponding to the MAP number of components are used for calculations.
#'
#' For examples see
#'
#' \url{http://faculty.missouri.edu/~micheasa/sppmix/sppmix_all_examples.html
#' #GetDensityValues}
#'
#' @param fit Object of class \code{damcmc_res} or \code{bdmcmc_res}.
#' @return A list containing the following components:
#' \item{Marginal}{the value of the Marginal (approximately)}
#' \item{LogLikelihood}{the value of the LogLikelihood}
#' \item{CompDensityAtXi}{the value of the component densities across all realizations and for each data point}
#' \item{DensityAtXi}{the value of the mixture density across all realizations and for each data point}
#' \item{EntropyMAP}{an approximation of the entropy of the distribution of the component indicators}
#' \item{Density}{the joint density at each posterior realization, i.e., for each iteration of \code{ps}, \code{mus} and \code{sigmas}.}
#' @author Sakis Micheas
#' @seealso \code{\link{normmix}},
#' \code{\link{est_mix_damcmc}},
#' \code{\link{est_mix_bdmcmc}},
#' \code{\link{rsppmix}}
#' @examples
#' \donttest{
#' # create the true mixture intensity surface
#' truesurf =normmix(ps=c(.2, .6,.2), mus=list(c(0.3, 0.3), c(0.7, 0.7),
#'  c(0.5, 0.5)),sigmas=list(.01*diag(2), .01*diag(2), .01*diag(2)),
#'  lambda=100,win=spatstat.geom::square(1))
#' plot(truesurf)
#' # generate the point pattern, truncate=TRUE by default
#' genPP=rsppmix(truesurf,truncate=FALSE)
#' fit=est_mix_damcmc(pp = genPP, m = 3)
#' allvals=GetDensityValues(fit)
#' MAPest=GetMAPEst(fit,vals=allvals)
#' plot(MAPest,main="IPPP intensity surface of MAP estimates")}
#'
#' @export
GetDensityValues<- function(fit)
{
  if(class(fit)=="bdmcmc_res")
  {
    tab=GetBDTable(fit,FALSE)
    fit=GetBDCompfit(fit,tab$MAPcomp)$BDgens
#    fit=drop_realization(fit,(fit$Badgen==1))
  }
  pp=fit$data
  n <- pp$n
  pattern <- cbind(pp$x,pp$y)
  vals=GetDensityValues_sppmix(pattern,
    fit,as.vector(pp$window$xrange),
    as.vector(pp$window$yrange))
  return(vals)
}

GetPriorVals<-function(
  pp,allgens,priortype,d,mu0,Sigma0,df0,sig0)
{
  #done in c++ code now faster
  return (GetPriorVals_sppmix(pp,allgens,priortype,d,mu0,Sigma0,df0,sig0))
}

#' Retrieve the IPPP likelihood value
#'
#' @description
#' Given a point pattern this function
#' calculates the IPPP likelihood value.
#'
#' For examples see
#'
#' \url{http://faculty.missouri.edu/~micheasa/sppmix/sppmix_all_examples.html
#' #GetIPPPLikValue}
#'
#' @param pp Point pattern object of class \code{ppp}.
#' @param surf IPPP intensity surface object of class \code{intensity_surface}.
#' @param truncate Logical variable indicating whether or not we
#' normalize the densities of the mixture components
#' to have all their mass within the window defined
#' in the point pattern \code{pp}.
#' @seealso \code{\link{est_mix_damcmc}},
#' \code{\link{rmixsurf}},
#' \code{\link{rsppmix}},
#' \code{\link{GetPMEst}}
#' @author Sakis Micheas
#' @examples
#' \donttest{
#' truemix_surf <- rmixsurf(m = 3, lambda=100,xlim = c(-3,3),ylim = c(-3,3))
#' plot(truemix_surf,main="True IPPP intensity surface")
#' genPPP=rsppmix(intsurf = truemix_surf, truncate = FALSE)
#' fit <- est_mix_damcmc(genPPP, m = 3)
#' MAPest=GetMAPEst(fit)
#' GetIPPPLikValue(genPPP,MAPest)
#' GetIPPPLikValue(genPPP,GetPMEst(fit))}
#'
#' @export
GetIPPPLikValue <- function(pp,surf,truncate=FALSE)
{
  n <- pp$n
  pattern <- cbind(pp$x,pp$y)
  approxcomp=rep(1,surf$m)
  if(truncate)
    approxcomp=approx_normmix(surf,xlim= as.vector(surf$window$xrange),
                              ylim = as.vector(surf$window$yrange))
  mix=MakeMixtureListFromNormMix(surf)
  den <-densNormMix_atxy_sppmix(pattern, mix,approxcomp)
  lognfac <- log(sqrt(2*pi*surf$lambda)) + n*log(n) - n
  val<- -lognfac + n*log(surf$lambda) - surf$lambda +
    sum(log(den))
  return(exp(val))
}
