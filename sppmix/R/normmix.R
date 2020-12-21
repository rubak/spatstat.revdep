#' Create a 2d mixture with normal components
#'
#' @description
#' Constructor function for the \code{normmix} class. Creates a mixture in
#' two dimensions with bivariate normal components. If the parameters
#' lambda and window are set, this function will create an intensity surface
#' object of class \code{intensity_surface}.
#'
#' For examples see
#'
#' \url{http://faculty.missouri.edu/~micheasa/sppmix/sppmix_all_examples.html
#' #normmix}
#'
#' @param ps Vector of component probabilities.
#' @param mus A list where every element is a vector of length 2, defining the
#' center of each component.
#' @param sigmas A list where every element is a 2 by 2 covariance matrix,
#' defining the covariance for each component.
#' @param lambda Optional parameter denoting the average number of points over the
#' window. If set along with the \code{win} parameter,
#' the returned object will be an intensity surface.
#' @param win Optional parameter for the window of
#' observation, an object of type \code{\link[spatstat.geom]{owin}}.
#' Must be set together with \code{lambda} in
#' order to create an intensity surface.
#' @param estimated Logical variable to indicate
#' that this is an estimated mixture not the true mixture surface.
#' By default it is set to FALSE, but when using the function to
#' define an mixture based on estimates of the parameters it
#' should be set to TRUE.
#' @param ... Additional arguments for the S3 method.
#' @return An object of class "normmix" containing the following components:
#'  \item{m}{Number of components.}
#'  \item{ps}{Vector of component probabilities.}
#'  \item{mus}{List of mean vectors of the components.}
#'  \item{sigmas}{List of covariance matrices of the components.}
#'  \item{lambda}{Returned only if lambda is provided when calling
#'  \code{normmix}.}
#'  \item{window}{Returned only \code{win} is provided when calling
#'  \code{normmix}. This is an object of class \code{\link[spatstat.geom]{owin}}.}
#'  \item{estimated}{Whether the normal mixture is estimated.}
#'
#' @seealso \code{\link{rnormmix}} for
#' generating a mixture with random parameters.
#' @author Yuchen Wang, Sakis Micheas
#' @examples
#' \donttest{
#' mix1 <- normmix(ps = c(.3, .7), mus = list(c(0.2, 0.2), c(.8, .8)),
#'  sigmas = list(.01*diag(2), .01*diag(2)))
#' mix1
#' summary(mix1)
#' intsurf1 <- normmix(ps = c(.3, .7), mus = list(c(0.2, 0.2), c(.8, .8)),
#'  sigmas = list(.01*diag(2), .01*diag(2)), lambda = 100, win = spatstat.geom::square(1))
#' intsurf1
#' summary(intsurf1)}
#'
#' @export
normmix <- function(ps, mus, sigmas, lambda = NULL, win = NULL,
                    estimated = FALSE) {
  if (abs(sum(ps) - 1) > .001) {
    cat("Sum of ps=",sum(ps))
    stop("Component probabilities must sum to 1.")
  }

  if (length(ps) != length(mus) | length(ps) != length(sigmas)) {
    stop("Number of components mismatch.")
  }

  if (!is.null(lambda) & spatstat.geom::is.owin(win)) {
    # generating intensity surface
    if (lambda <= 0) stop("Intensity must be greater than 0.")
    RVAL <- list(m = length(ps), ps = ps, mus = mus, sigmas = sigmas,
                 lambda = lambda, window = win, estimated = estimated)
    class(RVAL) <- c("intensity_surface", "normmix")
  } else {
    RVAL <- list(m = length(ps), ps = ps, mus = mus, sigmas = sigmas,
                 estimated = estimated)
    class(RVAL) <- "normmix"
  }

  RVAL
}

#' @method is normmix
is.normmix <- function(x,...) {
  # check class of object, if TRUE, assuming it's created by normmix() and
  # subject to all its constraints.
  ifelse(inherits(x, "normmix"), TRUE, FALSE)
}

#' @method is intensity_surface
is.intensity_surface <- function(x,...) {
  # check class of object, if TRUE, assuming it's created by normmix() and
  # subject to all its constraints.
  ifelse(inherits(x, "intensity_surface"), TRUE, FALSE)
}

#' @rdname normmix
#' @method print normmix
#' @export
print.normmix <- function(x,...) {
  cat("Normal Mixture with", x$m, "component(s).\n")
}

#' @description
#' The \code{print} function can be used on a \code{normmix}
#' or \code{intensity_surface} object in order to display basic information.
#' @param x An object of class \code{normmix} or \code{intensity_surface}.
#' @rdname normmix
#' @method print intensity_surface
#' @export
print.intensity_surface <- function(x,...) {
  pre <- ifelse(x$estimated, "Estimated average", "Expected")
  cat(pre, "number of points over the window:", x$lambda, "\n")
  print(x$window)
  NextMethod()
}

#' @method print summary_normmix
#' @export
print.summary_normmix <- function(x,...) {
  with(x,
       for (i in unique(comp)) {
         cat("Component", i, "is centered at", "[",
             round(mu[comp == i][1], 2), ",", round(mu[comp == i][2], 2), "]",
             "\nwith probability", round(ps[comp == i][1], 4), "\n")
         prmatrix(round(x[comp == i, 4:5], 4), rowlab = rep("", 2),
                  collab = c("covariance", "matrix:"))
       }
  )
}

#' @description
#' The \code{summary} function can be used on a \code{normmix}
#' or \code{intensity_surface} object in order to display additional information.
#'
#' @param object An object of class \code{normmix} or \code{intensity_surface}.
#' @rdname normmix
#' @method summary normmix
#' @export
summary.normmix <- function(object,...) {
  mix<-object
  comps <- vector("list", mix$m)
  for (i in 1:mix$m) {
    comps[[i]] <- data.frame(comp = i, ps = mix$ps[[i]],
                             mu = mix$mus[[i]], sigma = mix$sigmas[[i]])
  }
  RVAL <- do.call(rbind, comps)
  class(RVAL) <- c("summary_normmix", oldClass(RVAL))
  RVAL
}

#' @rdname normmix
#' @method summary intensity_surface
#' @export
summary.intensity_surface <- function(object,...) {
  print(object)
  NextMethod()
}

#' Generate a mixture with normal components
#'
#' @description
#' Generates a mixture on a 2d window where the means, covariances and component probabilities
#' are chosen randomly. The number of components can be either fixed or random.
#'
#' For examples see
#'
#' \url{http://faculty.missouri.edu/~micheasa/sppmix/sppmix_all_examples.html
#' #rnormmix}
#'
#' @param m Number of components of the mixture.
#' @param sig0 Tuning parameter for generating a random matrix from an Inverse Wishart
#' distribution. If this argument is missing it is set to .1 of the minimum width/height of the window.
#' @param df Degrees of freedom for generating a random matrix from an Inverse Wishart
#' distribution. Default is 10.
#' @param dvec A vector of weights used in the Dirichlet distribution used to sample the mixture probabilities. If the dimension of dvec is not
#' the same as the number of components, then dvec is either truncated to the same dimension or repeated to have dimension m. If missing, a vector of ones is used.
#' @param mu0,Sigma0 Mean and covariance matrix for a multivariate normal distribution, used to generate all component means. If mu0 is missing the center of the window
#' is used. If Sigma0 is missing it is set to the identity matrix. If both mu0 and Sigma0 are missing, the component means are generated uniformly over the window of observation.
#' @param rand_m Request a random number of components.
#' When \code{rand_m = TRUE}, the function will
#' randomly choose a number of components from
#' \code{1:m}.
#' @param xlim,ylim Vectors defining the observation window.
#' The component means are sampled uniformly over this window.
#'
#' @return Object of class \code{normmix}.
#' @author Sakis Micheas, Yuchen Wang
#' @examples
#' \donttest{
#' mix1 <- rnormmix(m = 3, sig0 = .1, df = 5)
#' summary(mix1)
#' mix2 <- rnormmix(m = 5, sig0 = .1, df = 5, rand_m = TRUE, ylim = c(0, 5))
#' summary(mix2)}
#'
#' @export
rnormmix <- function(m, sig0, df=10, rand_m = FALSE,
                     xlim = c(0, 1), ylim = c(0, 1),
                     dvec,mu0,Sigma0)
{
  if (rand_m)
  {
    # number of components is random
    m <- sample(1:m, 1)
  }
  if(missing(dvec))
    dvec=rep(1, m)
  else
    dvec=suppressWarnings(as.vector(matrix(data=dvec,nrow=m,ncol=1)))
  #  cat(dvec)
  gen_ps <- rDirichlet_sppmix(dvec)

  if(missing(mu0))
  {
    mu0=c(mean(xlim),mean(ylim))
    if(missing(Sigma0))
      gen_mu <- cbind(x = runif(m, xlim[1], xlim[2]),
                      y = runif(m, ylim[1], ylim[2]))
    else
      gen_mu <- rnorm2_sppmix(m,mu0,Sigma0)
  }
  else
  {
    if(missing(Sigma0))
      gen_mu <- rnorm2_sppmix(m,mu0,diag(2))
    else
      gen_mu <- rnorm2_sppmix(m,mu0,Sigma0)
  }

  if(missing(sig0))
    sig0=.1*min(c(xlim[2]-xlim[1],
                  ylim[2]-ylim[1]))
  gen_sigma <- stats::rWishart(m, df, sig0^(-2) * diag(2))
  mus <- vector(mode = "list", length = m)
  sigmas <- vector(mode = "list", length = m)

  for (k in 1:m) {
    mus[[k]] <- gen_mu[k, ]
    sigmas[[k]] <- solve(gen_sigma[, , k])
  }

  normmix(gen_ps, mus, sigmas)
}

#' Generate a Poisson process surface object
#'
#' @description
#' This function creates a Poisson point process intensity
#' surface modeled as a mixture of normal components, on the given
#' 2d window. The means, covariances and component probabilities
#' are chosen randomly based on parameters passed to the
#' function. The number of components can be either fixed or random.
#'
#' For examples see
#'
#' \url{http://faculty.missouri.edu/~micheasa/sppmix/sppmix_all_examples.html
#' #rmixsurf}
#'
#' @param m Number of components of the mixture. If omitted, m is uniformly selected from
#' 1 up to 10.
#' @param lambda Average number of points over the window. If omitted
#' lambda is generated from a Gamma with shape~Unif(1,10) and
#' scale~Unif(50,100).
#' @param dvec A vector of weights used in the Dirichlet distribution used to sample the mixture probabilities. If the dimension of dvec is not
#' the same as the number of components, then dvec is either truncated to the same dimension or repeated to have dimension m. If missing, a vector of ones is used.
#' @param mu0,Sigma0 Mean and covariance matrix for a multivariate normal distribution, used to generate all component means. If mu0 is missing the center of the window
#' is used. If Sigma0 is missing it is set to the identity matrix. If both mu0 and Sigma0 are missing, the component means are generated uniformly over the window of observation.
#' @param sig0 Tuning parameter for generating a random matrix from an Inverse Wishart
#' distribution.
#' @param df Degrees of freedom for generating a random matrix from an Inverse Wishart
#' distribution.
#' @param rand_m Request a random number of components.
#' When \code{rand_m = TRUE}, the function will
#' randomly choose a number of components from
#' \code{1:m}.
#' @param xlim,ylim Vectors defining the observation window.
#' The component means are sampled uniformly over this window.
#'
#' @return Object of class \code{intensity_surface}.
#' @author Sakis Micheas
#' @seealso \code{\link{plotmix_2d}},
#' \code{\link{summary.intensity_surface}},
#' \code{\link{plot.intensity_surface}}
#' @examples
#' \donttest{
#' mixsurf1 <- rmixsurf(m = 3, lambda=100)
#' summary(mixsurf1)
#' plot(mixsurf1)
#' plotmix_2d(mixsurf1)
#' mixsurf2 <- rmixsurf(m = 5, lambda=200, rand_m = TRUE, ylim = c(-3, 3))
#' summary(mixsurf2)
#' plot(mixsurf2)
#' plotmix_2d(mixsurf2)
#' mixsurf3 <- rmixsurf(m = 5, lambda=200, rand_m = TRUE, Sigma0=.01*diag(2))
#' summary(mixsurf3)
#' plot(mixsurf3)
#' plotmix_2d(mixsurf3)}
#'
#' @export
rmixsurf <- function(m,lambda, sig0, df, rand_m = FALSE,
                     xlim , ylim,
                     dvec,mu0,Sigma0)
{
  if(missing(m))
    m=sample(1:10,1)
  if(missing(lambda))
    lambda=rgamma(1,runif(1,1,10),scale=runif(1,50,100))
  if(missing(sig0))
    sig0=.1
  if(missing(df))
    df=5
  if(missing(xlim))
    xlim = c(0, 1)
  if(missing(ylim))
    ylim = c(0, 1)

  mix <- rnormmix(m, sig0, df,rand_m,xlim,ylim,dvec,mu0,Sigma0)

  intsurf <- normmix(mix$ps, mix$mus, mix$sigmas,
                     lambda,spatstat.geom::owin(xlim,ylim),FALSE)

  return(intsurf)
}

#' Convert a normal mixture to an intensity surface
#'
#' @description
#' This function converts a \code{normmix} object into an \code{intensity_surface}
#' object. It can also be used to change
#' the parameters \code{lambda} (average number of points over the window)
#' or \code{win} (window of observation)
#' of an \code{intensity_surface} object.
#'
#' For examples see
#'
#' \url{http://faculty.missouri.edu/~micheasa/sppmix/sppmix_all_examples.html
#' #to_int_surf}
#'
#' If the class of \code{mix} is \code{normmix},
#' \code{lambda} and \code{win} are used
#' to convert \code{mix} into an intensity
#' surface class. If the class of \code{mix}
#' is \code{intensity_surface} already, \code{lambda}
#' and \code{win} are used to change the
#' original settings for these parameters.
#'
#' @param mix Object of class \code{normmix} or \code{intensity_surface}.
#' @param lambda Optional parameter treated as the average number of points over the window.
#' @param win Optional parameter of class \code{\link[spatstat.geom]{owin}}, defining the window of observation.
#' @param return_normmix Logical variable requesting to return a normal mixture (discard \code{lambda}
#' and \code{win}).
#' @seealso \code{\link{normmix}},
#' \code{\link{rnormmix}},
#' \code{\link[spatstat.geom]{square}}
#' @return Object of class \code{intensity_surface}.
#' @author Yuchen Wang
#' @examples
#' \donttest{
#' truemix <- normmix(ps=c(.4, .2,.4), mus=list(c(0.3, 0.3), c(.5,.5),c(0.7, 0.7)),
#'  sigmas = list(.02*diag(2), .05*diag(2),.01*diag(2)))
#' intsurf=to_int_surf(truemix, lambda = 100, win = spatstat.geom::square(1))
#' #plot the true mixture
#' plot(intsurf,main = "True Poisson intensity surface (mixture of normal components)")
#' # using the demo_mix normmix object
#' summary(demo_mix)
#' demo_surf1=to_int_surf(demo_mix, lambda = 100, win = spatstat.geom::square(1))
#' plot(demo_surf1)
#' # using an intensity_surface object
#' summary(demo_intsurf)
#' demo_surf2=to_int_surf(demo_intsurf, win = spatstat.geom::square(2))
#' summary(demo_surf2)
#' plot(demo_surf2)
#' demo_surf3=to_int_surf(demo_intsurf, lambda = 50)
#' plot(demo_surf3)}
#'
#' @export
to_int_surf <- function(mix, lambda = NULL, win = NULL,
                        return_normmix = FALSE) {
  intsurf <- mix

  if (is.intensity_surface(mix)) {
    # input is intensity surface
    if (!missing(lambda)) {
      stopifnot(lambda > 0)
      intsurf$lambda <- lambda
    }
    if (!missing(win)) {
      stopifnot(spatstat.geom::is.owin(win))
      intsurf$window <- win
    }

  } else if (is.normmix(mix)) {
    if (!return_normmix) {
      # input is normmix, create intensity surface
      stopifnot(!missing(lambda) & !missing(win))
      intsurf <- normmix(mix$ps, mix$mus, mix$sigmas,
                         lambda = lambda, win = win, estimated = mix$estimated)
    }
  } else {
    stop("mix must be of class normmix or intensity surface")
  }

  if (return_normmix) {
    intsurf <- normmix(intsurf$ps, intsurf$mus, intsurf$sigmas,
                       estimated = intsurf$estimated)
  }

  intsurf
}

#' Quantify the difference between two surfaces
#'
#' @description
#' This function can be used to compare two
#' intensity surfaces. In particular, if we have the true
#' surface and an estimator of the truth, then we can assess how
#' well the estimate fits the surface, i.e., how
#' close are the two surfaces. Now if we have
#' two estimates of the true surface then the estimate
#' with the smallest measure fits the truth better. We can also compare
#' two estimating surfaces this way.
#'
#' For examples see
#'
#' \url{http://faculty.missouri.edu/~micheasa/sppmix/sppmix_all_examples.html
#' #CompareSurfs}
#'
#' @param surf1,surf2 Either IPPP intensity surfaces (objects of type \code{intensity_surface})
#' or images (objects of type \code{\link[spatstat.geom]{im}}) that represent intensity surfaces.
#' @param LL Length of the side of the square grid.
#' The intensities are calculated on an L * L grid.
#' The larger this value is, the slower the calculation,
#' but the better the approximation to the difference
#' between the two surfaces.
#' @param truncate Requests to truncate the components
#' of the mixture intensities to have all their mass
#' within the observation window.
#' @details Since the two surfaces passed to the function
#' can be represented as a 2d intensity surface, any measure between
#' two images can be used for comparison purposes, provided that the
#' window is the same.
#'
#' If the two windows are different
#' the function will choose the largest one
#' and compare the two surfaces in there.
#' @return Returns a list containing all distances
#' computed and the window of comparison, an object
#' of class \code{\link[spatstat.geom]{owin}}.
#' @author Sakis Micheas
#' @seealso \code{\link{rmixsurf}},
#' \code{\link[spatstat.geom]{owin}},
#' \code{\link{rsppmix}},
#' \code{\link{est_mix_bdmcmc}},
#' \code{\link{drop_realization}},
#' \code{\link{plotmix_3d}},
#' \code{\link{plot_avgsurf}},
#' \code{\link{GetBMA}},
#' @examples
#' \donttest{
#' #compare two surfaces first
#' mixsurf1 = rmixsurf(m = 5, sig0 = .1, df = 5,xlim= c(-1,4), ylim = c(-2,1), rand_m = FALSE)
#' mixsurf2 = rmixsurf(m = 8, sig0 = .1, df = 5,xlim= c(-4,3), ylim = c(-1,2), rand_m = FALSE)
#' comp=CompareSurfs(mixsurf1,mixsurf2)
#' plot(mixsurf1,main = "First IPPP",win=comp$win)
#' plot(mixsurf2,main = "Second IPPP",win=comp$win)
#' # now we compare intensity surfaces and image objects that represent intensity surfaces
#' truemixsurf = rmixsurf(m = 5,xlim= c(-2,2), ylim = c(-2,2))
#' plot(truemixsurf,main="True IPPP surface")
#' #get a point pattern
#' genpp = rsppmix(truemixsurf,truncate=FALSE)
#' # Run BDMCMC and get posterior realizations
#' postfit=est_mix_bdmcmc(genpp,m=10,L=30000)
#' postfit=drop_realization(postfit,.1*postfit$L) #apply burnin
#' BMA=GetBMA(postfit,burnin=0)
#' title1 = paste("Bayesian model average of",postfit$L,"posterior realizations")
#' plotmix_3d(BMA,title1=title1)
#' comp=CompareSurfs(truemixsurf,BMA,LL=100)
#' #We compare the average surface and the truth for many cases below. If we pass images, we
#' #make sure it has the same dimensions or we force it to the same value by setting LL=100.
#' #We retrieve the average surfaces corresponding to MAP-1, MAP and MAP+1 components and
#' #compare them against the truth.
#' #First retrieve the frequency table and MAP estimate for number of components
#' BDtab=GetBDTable(postfit)
#' BDtab
#' MAPm=BDtab$MAPcomp
#' BDfitMAPcomp_minus1=GetBDCompfit(postfit,MAPm-1,burnin=0)
#' avgsurfMAPcomp_minus1=plot_avgsurf(BDfitMAPcomp_minus1$BDgens, LL = 100,burnin=0)
#' comp=CompareSurfs(truemixsurf,avgsurfMAPcomp_minus1,LL=100)
#' BDfitMAPcomp=GetBDCompfit(postfit,MAPm,burnin=0)
#' avgsurfMAPcomp=plot_avgsurf(BDfitMAPcomp$BDgens, LL = 100,burnin=0)
#' comp=CompareSurfs(truemixsurf,avgsurfMAPcomp,LL=100)
#' BDfitMAPcomp_plus1=GetBDCompfit(postfit,MAPm+1,burnin=0)
#' avgsurfMAPcomp_plus1=plot_avgsurf(BDfitMAPcomp_plus1$BDgens, LL = 100,burnin=0)
#' comp=CompareSurfs(truemixsurf,avgsurfMAPcomp_plus1,LL=100)}
#'
#' @export
CompareSurfs<-function(surf1,surf2,LL = 100,
                       truncate = FALSE)
{
  if(missing(surf1)|missing(surf2))
  {
    stop("Pass two intensity surfaces or images.")
  }
  if(!(is.intensity_surface(surf1)|
       spatstat.geom::is.im(surf1))&
     !(is.intensity_surface(surf2)|
       spatstat.geom::is.im(surf2)))
  {
    stop("Pass intensity_surfaces or image objects.")
  }
  minx=10000000000000
  miny=10000000000000
  maxx=-10000000000000
  maxy=-10000000000000
  if(is.intensity_surface(surf1))
  {
    minx=min(c(surf1$window$xrange[1],minx))
    maxx=max(c(surf1$window$xrange[2],maxx))
    miny=min(c(surf1$window$yrange[1],miny))
    maxy=max(c(surf1$window$yrange[2],maxy))
  }
  if(is.intensity_surface(surf2))
  {
    minx=min(c(surf2$window$xrange[1],minx))
    maxx=max(c(surf2$window$xrange[2],maxx))
    miny=min(c(surf2$window$yrange[1],miny))
    maxy=max(c(surf2$window$yrange[2],maxy))
  }
  #if an image is passed get the larger dimension
  if(spatstat.geom::is.im(surf1) &
     spatstat.geom::is.im(surf2))
    if(surf1$dim[1]!=surf2$dim[1] |
       surf1$dim[2]!=surf2$dim[2] )
      stop("Dimension mismatch in the two images.")
  if(!spatstat.geom::is.im(surf1))
    dens1=dnormmix(surf1,xlim= xlim,
        ylim =ylim,L=LL, truncate=truncate)
  else
  {
    minx=min(c(surf1$xrange[1],minx))
    maxx=max(c(surf1$xrange[2],maxx))
    miny=min(c(surf1$yrange[1],miny))
    maxy=max(c(surf1$yrange[2],maxy))
    dens1=surf1
  }
  if(!spatstat.geom::is.im(surf2))
    dens2=dnormmix(surf2,xlim= xlim,
                 ylim =ylim,L=LL, truncate=truncate)
  else
  {
    minx=min(c(surf2$xrange[1],minx))
    maxx=max(c(surf2$xrange[2],maxx))
    miny=min(c(surf2$yrange[1],miny))
    maxy=max(c(surf2$yrange[2],maxy))
    dens2=surf2
  }
  xlim=rep(0,2)
  ylim=rep(0,2)
  xlim[1]=minx
  xlim[2]=maxx
  ylim[1]=miny
  ylim[2]=maxy
  alldist=rep(0,2)
  mat1=dens1$v
  mat2=dens2$v
  mat=mat1-mat2
  logmat=log(mat1)-log(mat2)
  alldist[1]=mean(mat^2)
  cat("\n1) Average Squared Pixel distance:",
      alldist[1],"\n")
  alldist[2]=max(abs(mat))
  cat("2) Maximum Absolute Pixel distance:",
      alldist[2],"\n")
  alldist[3]=MatrixNorm(mat,2.0)
  cat("3) Frobenius norm of the matrix of Pixel distances:",
      alldist[3],"\n")
  alldist[4]=mean(logmat^2)
  cat("4) Average Squared Log Pixel distance:",
      alldist[4],"\n")
  alldist[5]=MatrixNorm(logmat,2.0)
  cat("5) Frobenius norm of the matrix of Log Pixel distances:",
      alldist[5],"\n")
  cat("Window used: [",xlim[1],",",xlim[2],
      "]x[",ylim[1],",",ylim[2],"]")
  return(list(dists=alldist,win=spatstat.geom::owin(xlim,ylim)))
}

#' Counts points in a window
#'
#' @description
#' This function counts the number of points
#' from a point pattern within a specified window.
#'
#' For examples see
#'
#' \url{http://faculty.missouri.edu/~micheasa/sppmix/sppmix_all_examples.html
#' #Count_pts}
#'
#' @param pp Object of class \code{sppmix} or \code{\link[spatstat.geom]{ppp}}.
#' @param win The window of observation as an
#' object of class \code{\link[spatstat.geom]{owin}}, defining the window of observation.
#' @seealso \code{\link{rsppmix}},
#' \code{\link{rmixsurf}},
#' \code{\link{plotmix_2d}},
#' \code{\link[spatstat.geom]{square}},
#' \code{\link[spatstat.geom]{owin}}
#' @return An integer representing the number of points from \code{pp} within the window \code{win}.
#' @author Sakis Micheas
#' @examples
#' \donttest{
#' truemix_surf=rmixsurf(m = 3,lambda=100, xlim=c(-5,5),ylim=c(-5,5))
#' genPP=rsppmix(truemix_surf)
#' plotmix_2d(truemix_surf,genPP)
#' Count_pts(genPP,spatstat.geom::square(1))
#' Count_pts(genPP,spatstat.geom::square(2))
#' Count_pts(genPP,spatstat.geom::square(3))
#' Count_pts(genPP,spatstat.geom::square(4))
#' Count_pts(genPP,spatstat.geom::square(5))
#' Count_pts(genPP,spatstat.geom::owin(c(-5,5),c(-5,5)))}
#'
#' @export
Count_pts <- function(pp,win)
{
  xlims=win$xrange
  ylims=win$yrange
  pat=cbind(pp$x,pp$y)
  checkin=CheckInWindow_sppmix(pat,xlims,ylims,TRUE,TRUE);
  data=checkin$data_inW
  return(nrow(data))
}

#' Demo objects
#'
#' @description
#' Demo objects (mixture, surface and generated pattern)
#' using the classes provided by the \code{\link{sppmix}} package.
#'
#' For examples see
#'
#' \url{http://faculty.missouri.edu/~micheasa/sppmix/sppmix_all_examples.html
#' #demo_mix}
#'
#' @author Jiaxun Chen, Sakis Micheas, Yuchen Wang
#' @examples
#' \donttest{
#' demo_mix <- normmix(ps = c(.3, .7), mus = list(c(0.2, 0.2), c(.8, .8)),
#'  sigmas = list(.01*diag(2), .01*diag(2)))
#' demo_intsurf <- normmix(ps = c(.3, .7), mus = list(c(0.2, 0.2),
#'  c(.8, .8)),sigmas = list(.01*diag(2), .01*diag(2)), lambda = 100,
#'  win = spatstat.geom::square(1))
#' demo_genPPP<-rsppmix(demo_truesurf3, truncate=FALSE)}
#'
"demo_mix"

#' @rdname demo_mix
"demo_intsurf"

#' @rdname demo_mix
"demo_genPPP"

#' @rdname demo_mix
"demo_truemix3"

#' @rdname demo_mix
"demo_truemix3comp"

#' @rdname demo_mix
"demo_truesurf3"

#' @rdname demo_mix
"demo_intsurf3comp"
