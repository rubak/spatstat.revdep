#' Filters a set of points based on a probability map.
#'
#' This function rarefies a set of points using a raster map of probabilities of a point persisting in that location.
#'
#' @import raster
#' @importFrom classInt classIntervals
#' @importFrom stats runif
#'
#'
#' @param points data frame containing columns 'y' and 'x' (spatial coordinates on projected coordinate system).
#' @param MAP object of class RasterLayer. Map of probability [0;1] of a point to persist in any pixel.
#' @param Nclass optional number of classes into which MAP is going to be binned (quantile-based). This can be useful if the probability distribution of MAP is highly skewed.
#'
#' @return data frame with same structure as 'points' but with less rows (points).
#'
#' @author Luca Butikofer
#'
#' @export


spatFilter<-function(points, MAP, Nclass=FALSE){

  if(Nclass!=F){
    MAPVal<-extract(MAP,extent(MAP))
    Breaks<-classIntervals(MAPVal,n=Nclass,na.rm=T)$brks
    MAPVal<-.bincode(MAPVal, breaks=Breaks, right = F, include.lowest = F)
    MAPAt<-extract(MAP,data.frame(x=points$x,y=points$y),
                   cellnumbers=T,df=T)
    valAt<-MAPVal[MAPAt[,2]]
  }else{
    valAt<-extract(MAP,data.frame(x=points$x,y=points$y))
  }
  stay<-valAt>runif(length(valAt),min(valAt,na.rm=T),max(valAt,na.rm=T))
  points<-subset(points,stay=="TRUE")
  rownames(points)<-NULL
  return(points)
}
