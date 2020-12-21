#' Runs the EM algorithm.
#'
#' Attributes to the populations in an invasion time series a probability value of being of natural origin, as opposite of anthropogenic origin.
#'
#' @param dataset the data frame to be analised (WGS84, colums order should be "year","lat","long","origin")
#' @param sigma starting value for the standard deviation of the natural dispersal kernel (assumed to be a half normal)
#' @param pi starting value for the proportion of natural points in the dataset
#' @param randompoints data frame of 'y' and 'x' coordinates of random points (projected coordinate system)
#'
#' @return dataset argument with two additional colums. 'dist': distance from nearest point of natural origin or nearest anchor point (see Details); 'Pnat': probability of being of natural origin.
#'
#' @author Beatrix Jones, Luca Butikofer
#'
#' @export
#'
#' @examples
#' data('nzp')
#' data('frogs')
#' randp<- RPG(rpopn=1000, boundary=nzp, SP= 'random_frog')
#' frogsEM<- EM(dataset= frogs, randompoints= randp, sigma=6, pi=0.5)


EM<-function(dataset, randompoints, sigma, pi){

  dati<-dataset
  realp<-nrow(subset(dati,dati$Pnat==1))/nrow(dati)

  input.dist<-distances(dati)
  output.weights<-distance.by.year(input.dist)
  Kerns<-random.distances(dati,randompoints)
  update<-EM.max(output.weights[[1]], output.weights[[2]],sigma,pi,Kerns)

  dataset<-cbind(dataset, unlist(update[[1]]), unlist(output.weights[[1]]))
  names(dataset)[(ncol(dataset)-1):(ncol(dataset))]<-c("Pnat","Dist")

  return(dataset)
}
