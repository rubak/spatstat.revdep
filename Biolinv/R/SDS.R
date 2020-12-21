#' Generates a list of dispersal vectors.
#'
#' This function generates a set of spatial vectors whose length is sampled form a user-defined dispersal kernel and whose direction is sampled randomly from 1 degree to 360 degrees.
#'
#' @param probDist Probability Distribution of npop.
#' @param DIST vector of possible distances [m] to be sampled for dispersal.
#' @param npop number of distance-angle couples that need to be produced.
#'
#' @return data frame of y and x shifts.
#'
#' @author Luca Butikofer
#'
#' @export
#'
#' @examples
#' dist<- seq(.1, 30, .1)
#' prob<- fx(x=dist, a=7, c=2)
#' deltas<- SDS(DIST=dist, probDist=prob)



SDS<-function(probDist,DIST, npop=10000){


  dist<-1000*sample(DIST,npop,replace=T,prob=probDist)
  angles<-sample(seq(0,360,1),npop,replace=T)
  longitude<-cosd(angles)*dist
  latitude<-sind(angles)*dist
  deltas<-data.frame(y=latitude,x=longitude)
  return(deltas)
}
