#' Determination of the high-risk zone.
#'
#' \code{det_hrz_restr} determines the high-risk zone through the method of fixed radius
#' (type = "dist" and criterion = "direct"), the quantile-based method (type = "dist" and
#' criterion = "area"/"indirect") and the intensity-based method (type =  "intens").
#' Restriction areas can be taken into account.
#'
#' Used in functions eval_method, sim_clintens, sim_intens.
#' \cr
#' This function contains the same functionalities as \code{\link[highriskzone]{det_hrz}}.
#' In addition, it offers the possibility to take into account so-called restriction areas. This is relevant in
#' situations where the observed point pattern \code{ppdata} is incomplete. If it is known that no observations
#' can be made in a certain area (for example because of water expanses),
#' this can be accounted for by integrating a hole in the observation window.
#' The shape and location of the hole is given by \code{hole}, whereas \code{integratehole} is used to state
#' whether the hole is to become part of the resulting high-risk zone.
#' This may also be a reasonable approach if only few observations could be made in a certain area.
#' Another approach consists in weighting the observed events with their reciprocal observation probability when
#' estimating the intensity. To do so, the observation probability can be specified by using \code{obsprobs} (value of the
#' observation probability for each event) or \code{obsprobsimage} (image of the observation probability). Note that the
#' observation probability may vary in space.
#' \cr
#' If there are no restriction areas in the observation window, \code{\link[highriskzone]{det_hrz}}
#' can be used instead.
#' \cr
#' Note that for \code{criterion = "area"}, \code{cutoff} specifies the area of the high-risk zone outside the hole. If
#' \code{integratehole = TRUE}, the area of the resulting high-risk zone will exceed \code{cutoff}.
#' 
#'
#' For further information, Mahling et al. (2013) and Mahling (2013), Chapters 4 and 8 and Appendix A (References).
#'
#' @param ppdata  Observed spatial point process of class ppp.
#' @param type  Method to use, can be one of \code{"dist"}(method of fixed radius or quantile-based method), or
#'              \code{"intens"}(intensity based method)
#' @param criterion  criterion to limit the high-risk zone, can be one of
#'        \code{"area"} (giving size of hrz), \code{"indirect"} (giving quantile/alpha depending on type),
#'        or \code{"direct"} (giving radius/threshold c depending on type)
#' @param cutoff  Value of criterion (area, radius, quantile, alpha or threshold).
#'                 Depending on criterion and type.
#' @param hole  (optional) an object of class \code{owin} representing a region inside the observation window of
#'               the \code{ppdata} where no observations were possible.
#' @param integratehole    Should the \code{hole} be part of the resulting high-risk zone? Defaults to \code{TRUE}.
#' @param obsprobs  (optional)  Vector of observation probabilities associated with the observations contained in \code{ppdata}.
#'                              Must be given in the same order as the coordinates of the observations. Only meaningful
#'                              for the intensity-based method if some observations are located in areas where not all
#'                              events can actually be observed. For example, if only one third of the events in a specific region
#'                              could be observed, the observation probability of the corresponding observations
#'                              is 1/3.
#' @param obsprobimage  (optional) an object of class \code{im} giving the observation probabilities inside the
#'                        observation window. Ranges of the coordinates must equal those of \code{ppdata}.
#'                        Only used if \code{obsprobs} is not given.
#' @param distancemap  (optional) distance map: distance of every pixel to the nearest observation
#'                     of the point pattern; only needed for \code{type="dist"}. If not given,
#'                     it will be computed by \code{\link[spatstat.geom]{distmap}}.
#' @param intens  (optional) estimated intensity of the observed process (object of class "im",
#'                see \code{\link[spatstat.core]{density.ppp}}), only needed for type="intens". If not given,
#'                it will be estimated.
#' @param nxprob  Probability of having unobserved events.
#'                Default value is 0.1.
#' @param covmatrix  (optional) Covariance matrix of the kernel of a normal distribution, only needed for
#'                    \code{type="intens"} if no intensity is given. If not given, it will be estimated.
#' @param returnintens  Should the image of the estimated intensity be returned? Defaults to \code{TRUE}.
#' @export
#' @return An object of class "\code{highriskzone}", which is a list of
#'    \item{ typehrz, criterion, cutoff, nxprob}{ see arguments}
#'    \item{ zone }{ Determined high-risk zone: Object of class "owin" based on a binary mask.
#'                   See \code{\link[spatstat.geom]{owin}}. }
#'    \item{ threshold }{ determined threshold. If type = "dist" and criterion = "direct" it is the specified radius.
#'    If criterion = "indirect" or "area" the determined radius used to construct a risk zone fulfilling the specified criterion 
#'    and cutoff. If type = "dist" it is the specified or calculated threshold c, the maximum intensitiy of unexploded bombs 
#'    outside the risk zone.}
#'    \item{ calccutoff }{ determined cutoff-value. For type="dist" and criterion="area", this is the
#' quantile of the nearest-neighbour distance. For type="intens" and criterion="area" or "direct", it is the failure
#' probability alpha. For all other criterions it is NA.}
#'    \item{ covmatrix }{ If not given (and \code{type="intens"}), it is estimated. See \code{\link[ks]{Hscv}}.}
#'    \item{ estint }{ Estimated intensity. See \code{\link[spatstat.core]{density.ppp}}.}
#' @seealso \code{\link[spatstat.geom]{distmap}}, \code{\link[spatstat.geom]{eval.im}}, \code{\link[spatstat.geom]{owin}}
#' @examples
#' set.seed(1211515)
#' data(craterA)
#' #change npixel = 100 to 1000 to get a nicer picture
#' spatstat.geom::spatstat.options(npixel=100)
#' # reduce number of observations for faster computation
#' craterA <- craterA[sample(1:craterA$n, 150)]
#' # define restriction area
#' restrwin <- spatstat.geom::owin(xrange=craterA$window$xrange, yrange=craterA$window$yrange,
#'                            poly=list(x=c(1500, 1500, 2000, 2000), y=c(2000, 1500, 1500, 2000)))
#'
#' # create image of observation probability (30% inside restriction area)
#' wim <- spatstat.geom::as.im(craterA$window, value=1)
#' rim <- spatstat.geom::as.im(restrwin, xy=list(x=wim$xcol, y=wim$yrow))
#' rim$v[is.na(rim$v)] <- 0
#' oim1 <- spatstat.geom::eval.im(wim - 0.7 * rim)
#' # determine high-risk zone by weighting the observations
#' hrzi1 <- det_hrz_restr(ppdata=craterA, type = "intens", criterion = "indirect",
#'                        cutoff = 0.4, hole=NULL, obsprobs=NULL, obsprobimage=oim1, nxprob = 0.1)
#' 
#' # determine high-risk zone by accounting for a hole
#' hrzi2 <- det_hrz_restr(ppdata=craterA, type = "intens", criterion = "indirect",
#'                        cutoff = 0.4, hole=restrwin, obsprobs=NULL, obsprobimage=NULL, nxprob = 0.1)
 

det_hrz_restr <- function(ppdata, type, criterion, cutoff, hole=NULL, integratehole=TRUE, 
                          obsprobs=NULL, obsprobimage=NULL, distancemap=NULL, intens=NULL, 
                          nxprob=0.1, covmatrix=NULL, returnintens=TRUE){
  
  win <- ppdata$window
  calccutoff <- NA
  type <- match.arg(type, choices=c("intens", "dist"))  
  criterion <- match.arg(criterion, choices=c("area", "indirect", "direct"))
  check_det_hrz_restr_input(ppdata = ppdata, type = type, criterion = criterion, cutoff = cutoff,
                            hole=hole, integratehole=integratehole, obsprobs=obsprobs, obsprobimage=obsprobimage,
                      distancemap = distancemap, intens = intens, nxprob = nxprob, covmatrix = covmatrix)
  
  if (!is.null(hole)){
    winminus <- setminus.owin(A=ppdata$window, B=hole)
    ppdata <- ppp(x=ppdata$x, y=ppdata$y, window=winminus, marks=ppdata$marks)
  }
  
  if (!is.null(obsprobimage)){
    if (!is.null(hole)){
      holeim <- as.im(X=winminus, W=ppdata$window, value=0)
      obsprobimage <- eval.im(obsprobimage + holeim)
    }
    sumim <- summary(obsprobimage)
    if(sumim$min <= 0){warning("obsprobimage contains value 0 or negative values")}
  }

  # weights which are used for intensity estimation later on (their reciprocal value)
  # if neither obsprobs nor obsprobimage is given, all marks (weights) equal 1
  if(is.null(obsprobs)){
    if(!is.null(obsprobimage)){
      # Markierung mittels obsprobimage - hier noch nicht Kehrwert!
      obsprobs <- obsprobimage[ppdata]
    }else{obsprobs <- rep(1, times=ppdata$n)}
  } 
  #ppdata$marks <- obsprobs
  invweight <- obsprobs
  
  if(min(invweight) <= 0){stop("obsprobs contains value 0 or negative values")}
  
  # set the right values
  if(type=="dist"){
    
    if(is.null(distancemap)){distancemap <- distmap(ppdata)}
    
    if(criterion=="area"){
      
      res_det_radius <- det_radius(ppdata=ppdata, distancemap=distancemap, areahrz=cutoff, win=win)
      threshold <- res_det_radius$thresh
      calccutoff <- res_det_radius$cutoffdist
      
    }else{
      
      threshold <- ifelse(criterion=="indirect", quantile(nndist(ppdata), p=cutoff, type=8), cutoff)
      
    }
    
    HRZimage <- eval.im(distancemap < threshold)
    
  }
  
  if(type=="intens"){
    
    if(is.null(intens)){
      
      estim <- est_intens_weight(ppdata, covmatrix=covmatrix, weights=1/invweight)
      intens <- estim$intensest
      covmatrix <- estim$covmatrix
      
    }
    
    if(criterion=="area"){
      
      res_det_thresholdfromarea <- det_thresholdfromarea_rest(win=win, intens=intens, nxprob=nxprob, 
                                                              areahrz=cutoff, hole = hole, 
                                                              integratehole = integratehole)
      threshold <- res_det_thresholdfromarea$thresh
      calccutoff <- res_det_thresholdfromarea$calccutoff
      
    }else{
      
      threshold <- ifelse(criterion=="indirect", det_threshold(intens=intens, alpha=cutoff, nxprob=nxprob), 
                          (1-nxprob)/nxprob * cutoff)
      
    }
    
    HRZimage <- eval.im(intens > threshold)
    
    if(criterion == "direct")
      calccutoff <- det_alpha(intens, threshold = threshold, nxprob = nxprob)
    
    threshold <- nxprob / (1-nxprob) * threshold
  }
  
  if (is.null(hole) | !integratehole){
    totalim <- HRZimage
  }else{
    him <- as.im(hole)
    totalim <- eval.im(him | HRZimage)
    
  }

  Rwindow <- owin(xrange=win$xrange, yrange=win$yrange, mask=as.matrix(totalim))
  if(returnintens){
    result <- list(typehrz=type, criterion=criterion, cutoff=cutoff, nxprob = nxprob, zone=Rwindow, threshold=threshold, calccutoff=calccutoff, covmatrix=covmatrix, estint=intens)
  }else{
    result <- list(typehrz=type, criterion=criterion, cutoff=cutoff, nxprob = nxprob, zone=Rwindow, threshold=threshold, calccutoff=calccutoff, covmatrix=covmatrix)
  }
  class(result) <- "highriskzone"
  return(result)
  
}

#' Checks the arguments of det_hrz_restr
#' 
#' For each argument it is checked if it is of a correct value or class.
#' @inheritParams det_hrz_restr
#' 
check_det_hrz_restr_input <- function(ppdata, type, criterion, cutoff, hole, integratehole, 
                                      obsprobs, obsprobimage, distancemap, intens, 
                                      nxprob, covmatrix, returnintens){
  
  #errors
  #check if arguments have correct values
  if ( !is.ppp(ppdata) ) {
    stop("data is not of class ppp")
  }
  
  #this is for every criterion
  stopifnot( cutoff > 0 | length(cutoff) != 1 ) 
  #special: alpha and the quantile can only be in [0,1]
  if ( criterion == "indirect" ) stopifnot( cutoff < 1 ) 
  
  if ( !is.im(distancemap) & !is.null(distancemap) ) stop("distancemap must be of class im (see distmap) or NULL.")
  
  if(!is.null(hole) & !is.owin(hole)) stop("hole must be of class owin")
  stopifnot(is.logical(integratehole))
  if(!is.null(obsprobs))
    stopifnot(length(obsprobs) == ppdata$n | any(obsprobs>1 | obsprobs<0))
  
  if(!is.null(obsprobimage)) {
    stopifnot(is.im(obsprobimage))
    if(!is.null(obsprobs))
      warning("obsprobs and obsprobimage specified, obsprobs are used")
  }
  
  if ( !is.im(intens) & !is.null(intens) ) stop("wrong input for intens. The intensity must be of class im (see density.ppp) or NULL.")
  
  if ( !(nxprob > 0 & nxprob < 1) ) stop("nxprob is a probability and therefore it has to be in the interval [0, 1]")
  
  if ( !is.null(covmatrix) && (!isSymmetric(covmatrix) | !all(eigen(covmatrix)$values > 0)) ){
    stop("covmatrix has to be symmetric and positive semidefinit")
  
  stopifnot(is.logical(returnintens))  
  }
  
}
