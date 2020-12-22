#' Approximate the masses of bivariate normal mixture components
#'
#' @description
#' Calculates the mass of the density of each component of a normal mixture over a given window.
#'
#' For examples see
#'
#' \url{http://faculty.missouri.edu/~micheasa/sppmix/sppmix_all_examples.html
#' #approx_normmix}
#'
#' @param mix An object of class \code{normmix}
#' @param xlim,ylim Vectors defining the x-y integration limits. A mixture component mass is estimated within
#' this window.
#'
#' @return A numerical vector with elements corresponding to the mass of each component
#' within the window. These values are required when we
#' use truncation over an observation window in order to handle edge effects.
#' @seealso \code{\link{normmix}}
#' @author Jiaxun Chen, Yuchen Wang
#' @examples
#' \donttest{
#' truemix = normmix(ps = c(.3, .7), mus = list(c(0.2, 0.2),
#'  c(.8, .8)), sigmas =list(.01*diag(2),.01*diag(2)))
#' approx_normmix(truemix,xlim= c(-2, 2), ylim = c(-2, 2))}
#'
#' @export
approx_normmix <- function(mix, xlim = c(0, 1), ylim = c(0, 1)) {

  if (!is.normmix(mix)) {
    stop("mix must be of class normmix or intensity surface")
  }

  approx <- numeric(mix$m)

  for(k in 1:mix$m)
  {
    approx[k] <- mvtnorm::pmvnorm(lower = c(xlim[1], ylim[1]),  upper = c(xlim[2], ylim[2]),
                                  mean = mix$mus[[k]],sigma = mix$sigmas[[k]])
#    approx[k]=ApproxCompMass_sppmix(xlim,ylim,mix$mus[[k]],mix$sigmas[[k]])
  }
  return(approx)
}

#' Calculate the density or intensity of a
#' normal mixture over a fine grid
#'
#' @description
#' When a \code{normmix} object is given, this
#' function calculates the mixture density over
#' a fine grid for the given window. When an \code{intensity_surface} object
#' is given, the function multiplies the density with
#'  the intensity surface parameter \code{lambda}, and returns the Poisson
#'  mixture intensity function over the grid. Used for plotting.
#'
#' For examples see
#'
#' \url{http://faculty.missouri.edu/~micheasa/sppmix/sppmix_all_examples.html
#' #dnormmix}
#'
#' @param mix An object of class \code{normmix} or \code{intensity_surface}
#' @param xlim,ylim Vectors defining the x-y integration limits. A mixture component mass is estimated within
#' this window.
#' @param L Length of the side of the square grid.
#' The density or intensity is calculated on an L * L grid.
#' The larger this value is, the slower the calculation,
#' but the better the approximation.
#' @param truncate Requests to truncate the components
#' of the mixture intensity to have all their mass
#' within the given x-y limits.
#'
#' @return An object of class \code{\link[spatstat.geom]{im}}. This is a pixel image
#'  on a grid with values corresponding to the density (or intensity surface) at that location.
#'
#' @seealso \code{\link{normmix}},
#' \code{\link{rnormmix}},
#' \code{\link{plotmix_2d}},
#' \code{\link{plot_density}},
#' \code{\link{to_int_surf}}
#'
#' @author Jiaxun Chen, Sakis Micheas, Yuchen Wang
#' @examples
#' \donttest{
#' truemix <- rnormmix(m = 3, sig0 = .1, df = 5,xlim= c(0, 5),
#'  ylim = c(0, 5))
#' normdens=dnormmix(truemix,xlim= c(0, 5), ylim = c(0, 5))
#' #2d plots
#' plot_density(as.data.frame(normdens))+ ggplot2::ggtitle(
#'  "2d mixture density plot\nWindow=[0,5]x[0,5]")
#' plot_density(as.data.frame(normdens),TRUE)+ ggplot2::ggtitle(
#'  "2d mixture contour plot\nWindow=[0,5]x[0,5]")
#' #3d plot
#' plotmix_3d(normdens)
#' #Now build an intensity surface based on the normal mixture
#' intsurf=to_int_surf(truemix,lambda = 100, win =
#'  spatstat.geom::owin( c(0, 5),c(0, 5)))
#' intsurfdens=dnormmix(intsurf,xlim= c(0, 5), ylim = c(0, 5))
#' plot_density(as.data.frame(intsurfdens))+ ggplot2::ggtitle(
#'  "2d mixture intensity plot\nWindow=[0,5]x[0,5]")
#' plot_density(as.data.frame(intsurfdens),TRUE)+ ggplot2::ggtitle(
#'  "2d mixture intensity contour plot\nWindow=[0,5]x[0,5]")
#' plotmix_3d(normdens)#3d plot
#' #For an intensity surface object we use these functions instead
#' plotmix_2d(intsurf)
#' plot(intsurf)}
#'
#' @export
dnormmix <- function(mix, xlim = c(0, 1), ylim = c(0, 1), L = 128,
                     truncate = TRUE) {

  if (!is.normmix(mix)) {
    stop("mix must be of class normmix or intensity surface")
  }

  if (is.intensity_surface(mix)) {
    xlim <- mix$window$xrange
    ylim <- mix$window$yrange
  }

  x <- seq(xlim[1], xlim[2], length.out = L)
  y <- seq(ylim[1], ylim[2], length.out = L)

  approx=rep(1,mix$m)
  if (truncate) approx <- approx_normmix(mix, xlim, ylim)

  mixlist=MakeMixtureListFromNormMix(mix)
  est=matrix(0,L,L)
  for(i in 1:L)
    for(j in 1:L)
      est[j,i]=densNormMixatx_sppmix(c(x[i],y[j]),mixlist,approx)

  if (is.intensity_surface(mix)) {
    spatstat.geom::im(est * mix$lambda, x, y)
  } else {
    spatstat.geom::im(est, x, y)
  }
}

#' Compute the Bayesian Model average
#'
#' @description
#' This function uses the posterior realizations
#' from a \code{\link{est_mix_bdmcmc}} call, to compute the
#' Bayesian Model Average across different number of components
#' and returns the fitted Poisson point process with mixture of normals intensity surface.
#'
#' For examples see
#'
#' \url{http://faculty.missouri.edu/~micheasa/sppmix/sppmix_all_examples.html
#' #GetBMA}
#'
#' @param fit Object of class \code{bdmcmc_res}.
#' @inheritParams plot.bdmcmc_res
#' @return An image as an object of class \code{\link[spatstat.geom]{im.object}}.
#' @author Sakis Micheas
#' @seealso \code{\link{est_mix_bdmcmc}},
#' \code{\link{plotmix_3d}},
#' \code{\link{plot_density}}
#' @examples
#' \donttest{
#' fit=est_mix_bdmcmc(pp = spatstat.data::redwood, m = 5)
#' BMA=GetBMA(fit)
#' burnin=.1*fit$L
#' title1 = paste("Bayesian model average of",fit$L-burnin,"posterior realizations")
#' plotmix_3d(BMA,title1=title1)
#' plot_density(as.data.frame(BMA))+ggplot2::ggtitle("Bayesian model average intensity surface")
#' plot_density(as.data.frame(BMA),TRUE)+ggplot2::ggtitle(
#'  "Contours of the Bayesian model average intensity surface")}
#'
#' @export
GetBMA<- function(fit, win = fit$data$window,
                            burnin = fit$L/10,
                            LL = 100, zlims = c(0, 0))
{
  # Start the clock!
  ptm1 <- proc.time()
  L <- fit$L
  xlims <- win$xrange
  ylims <- win$yrange
#  cat("\nTo save all open rgl graphs use Save_AllOpenRglGraphs.\n")
  distr_numcomp <- GetCompDistr_sppmix(fit$numcomp[(burnin+1):L],fit$maxnumcomp)
  zcoord <- ApproxBayesianModelAvgIntensity_sppmix(
    fit$allgens_List[(burnin+1):L],
    fit$genlamdas[(burnin+1):L],
    fit$numcomp[(burnin+1):L],
    distr_numcomp,1,fit$maxnumcomp,LL, xlims, ylims,
    fit$ApproxCompMass[(burnin+1):L,])
  # Stop the clock
  ptm<-proc.time() - ptm1
  cat(paste("\nComputation time of the Bayesian Model Average in seconds:",ptm[[1]],"\n"))
  #find the highest z
  maxz_height <- max(zcoord)
  if (zlims[1] == 0 && zlims[2] == 0) {
    zlims <- c(0, 1.1*maxz_height)
  }
  gridvals=GetGrid_sppmix(LL,xlims,ylims);
  xcoord=as.vector(gridvals[[1]]);
  ycoord=as.vector(gridvals[[2]]);
  dens_image=spatstat.geom::as.im(list(x=xcoord,y=ycoord,z=zcoord))
  return(dens_image)
}
