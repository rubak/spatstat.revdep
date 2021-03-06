gen_n_from_mix <- function(n, mix) {
  comp <- sample(1:mix$m, size = n, replace = TRUE, prob = mix$ps)
  spp <- vector("list", length(unique(comp)))
  for (k in unique(comp)) {
    snd <- mvtnorm::rmvnorm(sum(comp == k),
                            mix$mus[[k]], mix$sigmas[[k]])
    spp[[k]] <- cbind(snd, k)
  }
  return(do.call(rbind, spp))
}

#' Generate a point pattern from a Poisson process
#'
#' @description
#' This function generates a point pattern from a Poisson
#' point process with intensity surface modeled by
#' a mixture of normal components.
#'
#' For examples see
#'
#' \url{http://faculty.missouri.edu/~micheasa/sppmix/sppmix_all_examples.html
#' #rsppmix}
#'
#' @param intsurf Object of class \code{intensity_surface} or \code{normmix}.
#' @param truncate Logical variable indicating that the points should be
#' within the window of observation. Default is TRUE.
#' @param marks An optional vector defining the mark space. A mark value is randomly selected and attached to the generated locations.
#' Default is NULL, so that we create an unmarked point pattern.
#' @param ... Further parameters passed to \code{to_int_surf()}.
#'
#' @return A point pattern of class \code{c("sppmix", "ppp")}. The object has
#' all the traits of the \code{\link[spatstat.geom]{ppp}} object and
#' in addition, a \code{comp} member indicating from which
#' mixture component the event comes from. The object members include:
#'
#' x : a vector of x coordinates of the events,
#'
#' y : a vector of y coordinates of the events,
#'
#' n : the number of events,
#'
#' window : the window of observation (an object of class \code{\link[spatstat.geom]{owin}}),
#'
#' marks : optional vector of marks,
#'
#' comp : vector of true allocation variables.
#'
#' @details If an intensity surface is passed to \code{intsurf}, the function
#' generates a pattern from the Poisson directly.
#' We can also pass a normal mixture of class \code{normmix}
#' and specify the observation window \code{win} and the parameter
#' \code{lambda}, which is interpreted as the average number of points over the
#' window, as additional parameters.
#'
#' Even if we pass an intensity surface to \code{rsppmix()},
#' we can still overwrite the \code{lambda} and \code{win} by
#' passing them as additional parameters. See the examples for
#' specific calls to the function.
#'
#' For a given window \code{W}, the number
#' of points \code{N(W)} in \code{W} follows
#' a Poisson distribution, with intensity measure
#' \code{Lambda(W)}. The intensity surface
#' of the Poisson process is the Radon-Nikodym
#' derivative of the measure \code{Lambda} with
#' respect to the Lebesgue measure.
#'
#' The intensity surface is modeled using
#' a multiple of a mixture of normal distributions,
#' i.e., the intensity is given by
#'
#' \code{f(x|lambda,thetas)=lambda * sum(p_i*f_i(x|mu_i,sigma_i))},
#'
#' where the parameters \code{thetas} consist
#' of the mixture probabilities \code{ps},
#' normal component means \code{mus}, and
#' covariances \code{sigmas}, with \code{sum(p_i)=1}.
#'
#' When the masses of all the components \code{f_i} are within the window, then
#' the mixture \code{sum(p_i*f_i(x|mu_i,sigma_i))}
#' integrates to 1 over \code{W}, so that the
#' average number of points is given by \code{E(N(W))=lambda}.
#'
#' If \code{truncate = TRUE}, we generate points from the
#' unbounded Poisson with the given mixture intensity function
#' until there are exactly \code{n} points
#' in the window (rejection method). If \code{truncate = FALSE},
#' the function will not check if the points are inside the window.
#' @author Jiaxun Chen, Sakis Micheas, Yuchen Wang
#' @seealso \code{\link{normmix}},
#' \code{\link{rmixsurf}},
#' \code{\link[spatstat.geom]{square}},
#' \code{\link{rsppmix}},
#' \code{\link{plot.sppmix}},
#' \code{\link{plotmix_2d}},
#' \code{\link{plot2dPP}},
#' \code{\link{plot.normmix}},
#' \code{\link{rnormmix}},
#' \code{\link{plotmix_3d}}
#' @examples
#' \donttest{
#' # create the true mixture
#' truemix_surf <- normmix(ps=c(.2, .6,.2), mus=list(c(0.3, 0.3), c(0.7, 0.7), c(0.5, 0.5)),
#'  sigmas = list(.01*diag(2), .03*diag(2), .02*diag(2)), lambda=100, win=spatstat.geom::square(1))
#' plot(truemix_surf)
#' # generate the point pattern
#' genPPP1=rsppmix(truemix_surf)
#' summary(genPPP1)
#' plot2dPP(genPPP1)
#' plot2dPP(genPPP1,truemix_surf$mus)
#' plotmix_2d(truemix_surf,genPPP1)
#' # overwrite lambda or win
#' genPPP2=rsppmix(truemix_surf, lambda = 200)
#' plotmix_2d(truemix_surf,genPPP2)
#' genPPP3=rsppmix(truemix_surf, win = spatstat.geom::square(2))
#' truemix_surf$window
#' plotmix_2d(truemix_surf,genPPP3)#will not see the points outside the surface window
#' plotmix_2d(truemix_surf,genPPP3, win = spatstat.geom::square(2)) #have to pass the new window
#' #to see the points
#' #use normmix with additional parameters
#' truemix<- rnormmix(m = 3, sig0 = .1, df = 5, xlim= c(0, 3), ylim = c(0, 3))
#' plot(truemix)
#' normdens=dnormmix(truemix, xlim= c(0, 3), ylim = c(0, 3))
#' plotmix_3d(normdens)
#' genPPP4=rsppmix(truemix, lambda = 100, win = spatstat.geom::square(3))
#' # turn off truncation
#' genPPP5=rsppmix(intsurf = truemix_surf, truncate = FALSE)
#' plotmix_2d(truemix_surf,genPPP5)
#' plotmix_2d(truemix_surf,genPPP5, win = spatstat.geom::square(2))
#' plotmix_2d(truemix_surf,genPPP5,contour=TRUE)
#' intsurf6=rmixsurf(m=5,lambda=rgamma(1,shape=10,scale=5),
#'  df=5,sig0=1,rand_m=TRUE,mu0 = c(.5,.5),Sigma0 = 0.001*diag(2))
#' genPPP6=rsppmix(intsurf6,marks=1:3,truncate = FALSE)
#' plotmix_2d(intsurf6,genPPP6)
#' plot(genPPP6,showmarks=TRUE)}
#' @export
rsppmix <- function(intsurf, truncate = TRUE,
                    marks=NULL,...)
{

  intsurf <- to_int_surf(intsurf, ...)
  win <- intsurf$window

  n <- rpois(1, intsurf$lambda)

  if (n == 0) {
    stop("lambda value too small, no points generated.")
  }

  spp <- gen_n_from_mix(n, intsurf)

  if (truncate == TRUE) {
    while (sum(spatstat.geom::inside.owin(spp[, 1], spp[, 2], win)) < n) {
      spp <- rbind(spp, gen_n_from_mix(n, intsurf))
    }
    spp <- spp[spatstat.geom::inside.owin(spp[, 1], spp[, 2], win), ][1:n, ]
  } else {
    warning(paste(sum(!spatstat.geom::inside.owin(spp[, 1], spp[, 2], win)),
                  "points are outside the window."))
  }
  if(!is.null(marks))
  {
    if(!is.vector(marks))
      stop("marks argument is not a vector")
    msize=length(marks)
    tags=sample(x=1:msize,size=n,replace=TRUE)
    markvalues=rep(0,n)
    for(i in 1:n)
      markvalues[i]=marks[tags[i]]
    RVAL <- ppp(x=spp[, 1],y=spp[, 2],window=win, marks=markvalues, check = truncate)
    RVAL$markformat="vector"
  }
  else
    RVAL <- as.ppp(spp[, 1:2], W=win, check = truncate)
  RVAL$comp <- spp[, 3]
  class(RVAL) <- c("sppmix", "ppp")
  return(RVAL)
}

#' @rdname rsppmix
#' @param object A point pattern object of class \code{sppmix}.
#'
#' @export
#' @method summary sppmix
summary.sppmix <- function(object,...) {
  print(object)
  cat('\nFrequency table of the allocated number of points to each component')
  print(table(object$comp))
}
