#' Generates random points within a geographic boundary.
#'
#' This function generates random points within a geographic area set as either sp::SpatialPolygons or raster::RasterLayer.
#'
#' @importFrom sp spsample coordinates
#'
#' @param rpopn total number of random points to be generated.
#' @param boundary either a RasterLayer object (0: inside; >= 1: outside) or a SpatialPolygon(s) object.
#' @param SP name of the species being simulated for labelling the output.
#' @param year year name (used by wrapper function simulacro() to lable consecutive time steps in its simulate time series).
#'
#' @return data frame of random points.
#'
#' @author Luca Butikofer
#'
#' @export


RPG<-function(rpopn,boundary,year="anyYear",SP){


  if(class(boundary)=="RasterLayer"){
    boundary<-as.data.frame(boundary,xy=T)
    random<-boundary[boundary[,3]==0,][
      sample(nrow(subset(boundary,boundary[,3]==0)),rpopn),]
    random2<-as.numeric(row.names(random))
    random3<-boundary[random2,2:1]
  } else if (class(boundary)=='SpatialPolygons'){
    random<-spsample(boundary,rpopn,type='random')
    random2<-data.frame(coordinates(random))
    random3<-random2[,2:1]
  }
  random3$year<-rep(year,rpopn)
  random3$species<-rep(SP,rpopn)
  random3$Pnat<-rep(0,rpopn)
  random3$Dist<-rep("artificial",rpopn)
  random3<-random3[,c(3,1,2,4,5,6)]
  rownames(random3)<-NULL
  return(random3)
}
