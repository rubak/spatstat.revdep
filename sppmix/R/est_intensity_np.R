#' Estimate the intensity surface using a non-parametric method
#'
#' @description
#' Using an Epanechnikov kernel we calculate an estimate of
#' intensity surface while accounting for edge effects.
#'
#' For examples see
#'
#' \url{http://faculty.missouri.edu/~micheasa/sppmix/sppmix_all_examples.html
#' #est_intensity_np}
#'
#' @param pattern A two-dimesional spatial point pattern in
#' the form of a \code{\link[spatstat.geom]{ppp}} object.
#' @param win Object of class \code{\link[spatstat.geom]{owin}}.
#' @param h Kernel bandwidth. \code{h} should be a positive number.
#' @param L Length of the side of the square grid.
#' @param kernel Kernel used to estimate the intensity surface.  Currently, only
#'  supports the Epanechnikov kernel.
#' @param edgecorrect Logical flag indicating whether or not
#' to use edge-correction in the estimating intensity surface. The default is TRUE.
#' @param truncate Logical flag indicating whether or not to use
#' points only within the window. The default is TRUE.
#'
#' @return An object of class \code{\link[spatstat.geom]{im}}.
#' @seealso \code{\link{rnormmix}},
#' \code{\link{to_int_surf}},
#' \code{\link{rsppmix}},
#' \code{\link{plotmix_3d}}
#' @author Jiaxun Chen, Sakis Micheas
#' @examples
#' \donttest{
#' mix1 <- rnormmix(5, sig0 = .01, df = 5, xlim=c(0, 5),ylim=c(0, 5))
#' intsurf1=to_int_surf(mix1,lambda = 40, win =spatstat.geom::owin( c(0, 5),c(0, 5)))
#' plot(intsurf1)
#' pp1 <- rsppmix(intsurf1)
#' # estimate and plot the estimated intensity surface
#' surfNP1 <- est_intensity_np(pp1, win=spatstat.geom::owin( c(0, 5),c(0, 5)), h=0.05,
#'  L=100,truncate = FALSE)
#' plotmix_3d(surfNP1,title1="Non parametric estimator of the intensity surface")
#' #truncate components to have all their mass in the window
#' surfNP2 <- est_intensity_np(pp1, win=spatstat.geom::owin( c(0, 5),c(0, 5)), h=0.5, L=100)
#' plotmix_3d(surfNP2,title1="(Truncated) Non parametric estimator of the intensity surface")}
#'
#' @export
est_intensity_np <- function(pattern, win, h, L=10, kernel=c("Epanechnikov"),
                             edgecorrect=TRUE, truncate=TRUE){
  kernel <- match.arg(kernel)
  if (truncate==TRUE) {
    pattern <- pattern[spatstat.geom::inside.owin(pattern$x,pattern$y,win)]
  }
  x <- seq(win$xrange[1],win$xrange[2],length.out = L)
  y <- seq(win$yrange[1],win$yrange[2],length.out = L)
  lenx <- length(x)
  leny <- length(y)
  intensity <- matrix(0,lenx,leny)
  loc <- expand.grid(x,y)
  pp <- cbind(pattern$x,pattern$y)
  index <- expand.grid(1:pattern$n, 1:(lenx*leny))
  xy <- loc[index[,2],]-pp[index[,1],]
  if (edgecorrect==TRUE){
    LL <- 20
    ax <- seq(win$xrange[1],win$xrange[2], length.out = LL)
    ay <- seq(win$yrange[1],win$yrange[2], length.out = LL)
    alenx <- length(ax)
    aleny <- length(ay)
    area <- (ax[2] - ax[1]) * (ay[2] - ay[1])
    centers <- expand.grid(x_center=(ax[1:(alenx-1)] + ax[2:alenx])/2,
                           y_center=(ay[1:(aleny-1)] + ay[2:aleny])/2)
    index.mass <- expand.grid(1:((alenx-1)*(aleny-1)),1:(lenx*leny))
    dist.centers <- loc[index.mass[,2],]-centers[index.mass[,1],]
    if(kernel=="Epanechnikov"){
      diff <- rowSums((dist.centers)^2)
      mass <- ifelse(diff < 1, 4*(1-diff) / (2*pi), 0)
      edgecor <- area * aggregate(mass,list(index.mass[,2]),sum) / (h^2)
      quad <- rowSums((xy/h)^2)
      val <- ifelse(quad < 1, 4*(1 - quad) / (2*pi), 0)
      intensity <- aggregate(val, list(index[, 2]), sum) / (h^2*edgecor)
    }
  } else {
    if(kernel=="Epanechnikov"){
      quad <- rowSums((xy/h)^2)
      val <- ifelse(quad < 1, 4*(1 - quad) / (2*pi), 0)
      intensity <- aggregate(val, list(index[, 2]), sum) / h^2
    }
  }
  return(spatstat.geom::im(matrix(intensity[,2],lenx,leny,byrow=T),x,y))
}

