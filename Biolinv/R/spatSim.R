#' Computes spatial dissimilarity of two point processes.
#'
#' This function uses Rippley's K-function (see Details) to compute spatial dissimilarity of two point processes.
#'
#' @importFrom spatstat ppp owin Kest
#'
#' @param MOD0 data frame containing 'y' and 'x' coordinates (projected coordinate system) to compare with MOD2.
#' @param MOD2 data frame containing 'y' and 'x' coordinates (projected coordinate system) to compare with MOD0.
#' @param WINDOW window of observation of the point patterns (MOD0 and MOD2)(see ?spatstat.geom::owin).Must be an object of class 'owin'.
#' @param R numeric vector of searc distances for the K-function.
#'
#' @return squared sum of distances between K-functions. This is a measure of spatial dissimilarity.
#'
#' @author Luca Butikofer
#'
#' @export
#'
#' @examples
#' ran<- data.frame('y'=sample(1000),'x'=sample(1000))
#' nor<- data.frame('y'=rnorm(1000,sd=150,mean=500),'x'=rnorm(1000,sd=150,mean=500))
#' window<-  spatstat.geom::owin(xrange=c(0,1000),yrange=c(0,1000))
#' spatSim(ran, nor, WINDOW= window, R=0:200)


spatSim<-function(MOD0,MOD2,WINDOW,R){

  mod0<-ppp(MOD0$x,MOD0$y,window=WINDOW)    #MEASURED DATA
  mod2<-ppp(MOD2$x,MOD2$y,window=WINDOW)    #SIMULATED DATA

  K0<-Kest(mod0,r=R)
  K2<-Kest(mod2,r=R)
  DK0<-K0$iso-K2$iso    #Delta K
  sqDK0<-sum(DK0)^2    #Squared sum of Delta K

  return(sqDK0)
}
