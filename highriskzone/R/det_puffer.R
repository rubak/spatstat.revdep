#' Estimation of width of a guard region given an estimated highriskzone
#'
#' \code{det_guard_width} determines the necessary width of a guard region in which
#' the existence of additional observed bomb craters could change a intensity based estimated 
#' highriskzone within the evaluation area of interest.
#' Within the evaluation area, the high risk zone consists of all points at which the estimated 
#' intensity of unexploded bombs exceeds a certain, specified or estimated threshold c. At a given 
#' point s, the intensity of unexploded bombs is given by the sum of all evaluated bivariate normal kernels
#' centered at the observed bomb craters multiplied by a constant nxprob/1-nxprob.
#' If the estimated intensity of unexploaded bombs is zero at a point at the boarder of the evaluation area
#' an additional observation outside the area could lift the intensity only above the determined threshold if the 
#' distance to the boarder is small enough so that the density of the normal kernel (which is centered at the additional 
#' observation) is bigger than the threshold at the boarder (assuming that the estimated kernel doesn't change due to the 
#' additional observation). The function returns the biggest distance in which it is possible that the density of
#' the bivariate normal kernel of the intensity of the supplied highriskzone exceeds thresh_const times the threshold
#' of the highriskzone. If thresh_const is set to 1, the guard region is the smallest region with constant width around 
#' the evaluation area in which a single additional observation could (but not necessarily does) increase the 
#' highriskzone within the evaluation area at a point at the boarder if the intensity of unexploaded bombs was zero at this point before. 
#' If the intensity was >0 at a point at the boarder of the evaluation area, or more than 1 additional observations
#' are found nearby outside of the evaluation area, the highriskzone within the evaluation area could already expand  
#' by addditional observations with a bigger distance from the boarder. This can be considered by setting thresh_const < 1, 
#' which intuitively means that 1/thresh_const crater observation at the same point could expand the highriskzone within 
#' the evaluation area in the direction of the additional observations, or that a point the boarder becomes part of the highriskzone
#' by the observation of a single additional crater if the intensity at this point was thresh_cont times the highriskzone threshold
#' based on all crater observations within the evaluation area.
#' 
#' For more infos on the construction of guard zones see Mahling (2013, Appendix B, Approach 2)
#'  
#' 
#' @param highriskzone the estimated highriskzone for the evaluation area
#' @param thresh_const the constant multiplied with the determined threshold, 0 < thresh_const < 1. 
#' @return The constant width of the guard region.
#' @importFrom mvtnorm dmvnorm
#' @export
#' @examples
#' ## change npixel to 1000 to obtain nicer plots
#' spatstat.geom::spatstat.options(npixel=100)
#' data(craterA)
#' # reduce number of observations for faster computation
#' thin.craterA <- craterA[1:50]
#' hrzi1 <- det_hrz(thin.craterA, type = "intens", criterion = "area", cutoff = 100000, nxprob = 0.1)
#' det_guard_width(hrzi1, thresh_const = .25)

det_guard_width <- function(highriskzone, thresh_const = .5) {
  if (class(highriskzone)[1] != "highriskzone")
    stop("highriskzone has to be of class highriskzone!")
  if (thresh_const < 0 | thresh_const >1)
    stop("thres_const has to be in [0,1]")
  if(is.null(highriskzone$covmatrix)) 
    stop("highriskzone has to be estimated based on an intensity based approach and information on (estimated) kernel covariance is necessary")
  
  cov_hat <- highriskzone$covmatrix
  eigen_cov <- eigen(cov_hat)
  
  cutval <- thresh_const * (1-highriskzone$nxprob) / highriskzone$nxprob * highriskzone$threshold
  if (dmvnorm(c(0,0), c(0,0), cov_hat) < cutval) {
    buffer_length <- 0
  } else {
    # for a given distance of center d, highest density of bivariate normal 
    # is in direction of major axis, determined by first eigenvector of covmat
    # biggest distance of center and a coordinate with density y, is therefore 
    # euclidean distance of center and coordinate on major axis
    # with density y
    maj_axis_dir <- sqrt(eigen_cov$values[1]) * eigen_cov$vectors[,1]
    mult_fac <- uniroot(function(x) dmvnorm(maj_axis_dir * x, 
                                            c(0,0), cov_hat) - cutval, 
                        c(0, 2*max(cov_hat)))
    buffer_length <- sqrt(sum((mult_fac$root * maj_axis_dir)^2))
  }
buffer_length
}
