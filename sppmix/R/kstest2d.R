#' Nonparametric Goodness-of-fit test between two point patterns
#'
#' @description
#' This function performs a two-dimensional Kolmogorov-Smirnov goodness-of-fit
#' test on two point patterns.
#'
#' For examples see
#'
#' \url{http://faculty.missouri.edu/~micheasa/sppmix/sppmix_all_examples.html
#' #kstest2d}
#'
#' @param x1,x2 Objects of class \code{\link[spatstat.geom]{ppp}}.
#' @param showinfo Logical variable. Requests to display the test conclusion based on the value of the p-value. Default is TRUE.
#' @return A list with class "htest" containing the following components:
#' \item{statistic}{Value of the KS statistic}
#' \item{p.value}{The p-value of the test}
#' \item{alternative}{A character string describing the alternative hypothesis}
#'
#' @references Peacock, J.A. (1983). Two-dimensional goodness-of-fit testing in
#' astronomy. Monthly Notices Royal Astronomy Society, 202, 615-627.
#'
#' Adapted from Matlab code by Dylan Muir.
#' @author Jiaxun Chen, Sakis Micheas
#' @seealso \code{\link{rnormmix}},\code{\link{to_int_surf}},\code{\link[spatstat.geom]{owin}}
#' @examples
#' \donttest{
#' # generate two point patterns
#' mix1 <- rnormmix(3, sig0 = .01, df = 5, xlim=c(0, 5), ylim=c(0, 5))
#' intsurf1=to_int_surf(mix1,lambda = 40, win =spatstat.geom::owin( c(0, 5),c(0, 5)))
#' mix2 <- rnormmix(8, sig0 = .01, df = 10, xlim=c(0, 5),ylim=c(0, 5))
#' intsurf2=to_int_surf(mix2,lambda = 50, win =spatstat.geom::owin( c(0, 5),c(0, 5)))
#' #generate patterns from the two different models
#' pp1 <- rsppmix(intsurf1)
#' pp2 <- rsppmix(intsurf2)
#' pp3 <- rsppmix(intsurf2)#pp3 is from the same model as pp2
#' # Test for goodness of fit, p-value should be small
#' kstest2d(pp1, pp2)
#' # Test for goodness of fit, p-value should be large
#' kstest2d(pp2, pp3)}
#'
#' @export
kstest2d <- function(x1, x2,showinfo=TRUE) {
  ecdf2d <- function(x,edge){
    count=cbind(as.numeric(x$x>=edge[1] & x$y>=edge[2]),
                as.numeric(x$x<=edge[1] & x$y>=edge[2]),
                as.numeric(x$x<=edge[1] & x$y<=edge[2]),
                as.numeric(x$x>=edge[1] & x$y<=edge[2]))
    return(count)
  }

  if((!spatstat.geom::is.ppp(x1)) | (!spatstat.geom::is.ppp(x2))){
    stop("Point pattern must be of class ppp.")
  }
  n1 <- length(x1$x)
  n2 <- length(x2$x)
  pp1 <- cbind(x1$x,x1$y)
  pp2 <- cbind(x2$x,x2$y)
  KSstat=-Inf
  for (i in 1:(n1+n2)){
    if(i <= n1) {
      edge=pp1[i,]
    } else {
      edge=pp2[i-n1,]
    }
    vfcdf1=apply(ecdf2d(x1,edge), 2, sum)/n1
    vfcdf2=apply(ecdf2d(x2,edge), 2, sum)/n2
    vfdiff=abs(vfcdf2-vfcdf1)
    fksts=max(vfdiff)
    if (fksts>KSstat) KSstat=fksts
  }

  #Peacock Z calculation and P estimation
  n <- n1*n2/(n1+n2)
  Zn <- sqrt(n) * KSstat
  Zinf <- Zn/(1-0.53 * n^(-0.9))
  pValue=2*exp(-2 * (Zinf - 0.5)^2)
  # Clip invalid values for P
  if (pValue > 1.0)
    pValue = 1.0;


  RVAL <- list(statistic=c(KSstat = KSstat),p.value = pValue,method =
                 "Two sample 2D Kolmogorov-Smirnov Test",data.name =
                 paste(sep=", ",deparse(substitute(x1)),deparse(substitute(x2))
                       ),alternative="No goodness of fit")
  class(RVAL) <- "htest"
  if(showinfo)
  {
    if (RVAL$p.value > .1)
    {
      cat("\nLarge p-value, the patterns arise from the SAME model.\n")
    }
    else
    {
      cat("\nSmall p-value, the patterns arise from DIFFERENT models.\n")
    }
  }
  return(RVAL)
}

#' Nonparametric Goodness-of-fit test for a point pattern against a surface
#'
#' @description
#' This function performs a two-dimensional Kolmogorov-Smirnov goodness-of-fit
#' test for a point pattern against a given intensity surface.
#'
#' For examples see
#'
#' \url{http://faculty.missouri.edu/~micheasa/sppmix/sppmix_all_examples.html
#' #kstest2dsurf}
#'
#' @param pp Object of class \code{\link[spatstat.geom]{ppp}}.
#' @param intsurf Object of class \code{intensity_surface}.
#' @param truncate Requests to truncate the generated point patterns
#' to be within the window of the intensity object intsurf. Default is FALSE.
#' @param iters Number of point patterns to generate and
#' compare against \code{pp}. The larger this
#' value is, the more test performed, and thus the
#' more reliable the result.
#'
#' @author Sakis Micheas
#' @seealso \code{\link{rmixsurf}},
#' \code{\link{kstest2d}},
#' \code{\link{rsppmix}},
#' \code{\link{plotmix_2d}}
#' @examples
#' \donttest{
#' # generate two intensity surfaces; assume the same window [-3,3]x[-3,3]
#' mixsurf1 <- rmixsurf(m = 3, lambda=100,xlim=c(-3,3),ylim=c(-3,3))
#' plot(mixsurf1)
#' mixsurf2 <- rmixsurf(m = 5, lambda=200,xlim=c(-3,3),ylim=c(-3,3))
#' plot(mixsurf2)
#' #generate point patterns from the two different models
#' pp1 <- rsppmix(mixsurf1, truncate=FALSE)
#' plotmix_2d(mixsurf1,pp1,colors=TRUE)
#' pp2 <- rsppmix(mixsurf2, truncate=FALSE)
#' plotmix_2d(mixsurf2,pp2,colors=TRUE)
#' # Test for goodness of fit, p-value should be small
#' kstest2d(pp1, pp2)
#' # Test each pattern for gof against both Poisson models
#' kstest2dsurf(pp1, mixsurf1)#correct model for pp1
#' kstest2dsurf(pp1, mixsurf2)#wrong model for pp1
#' kstest2dsurf(pp2, mixsurf2)#correct model for pp2
#' kstest2dsurf(pp2, mixsurf1)#wrong model for pp2}
#'
#' @export
kstest2dsurf <- function(pp, intsurf,truncate = FALSE,iters=500)
{

  if((!spatstat.geom::is.ppp(pp)) || (!is.intensity_surface(intsurf)))
  {
    stop("Bad point pattern object or intensity surface object.")
  }
  gof=0
  for(i in 1:iters)
  {
    pp1 <- suppressWarnings(rsppmix(intsurf,truncate = truncate))
    res=kstest2d(pp, pp1,showinfo=FALSE)
    if(res$p.value>.1)
      gof=gof+1
  }
  cat("K-S suggests gof",100*gof/iters,"% of the time.\n",iters,"K-S tests performed.")
}
