\name{test.intersection}
\alias{test.auto.intersection}
\alias{test.intersection}
\alias{test.auto.intersection.p}
\alias{test.intersection.p}
\title{
Test Overlapping of Polygons Against Random Rotation
}
\description{
These functions test the overlaping surface area of a colection of polygons against a a null model of random rotation. \code{test.auto.intersection} test for overlapping between polygons of the same type (i.e., the same species) whereas \code{test.intersection} test for overlapping between polygons of different types (i.e., between different species).  \code{test.auto.intersection.p} and \code{test.intersection.p} functions are parallelized versions which might significantly  reduce computing times.
}
\usage{
test.auto.intersection(ventana1, nsim = 199, win = NULL,
    prop=1, centroides1= NULL,    diametros1 =NULL)
test.intersection(ventana1, ventana2, nsim, prop = 1, win = NULL,
   centroides1 = NULL, diametros1 = NULL, centroides2 = NULL, 
   diametros2 = NULL)
test.auto.intersection.p(ventana1, nsim = 199, win = NULL,
    prop=1, centroides1= NULL,    diametros1 =NULL, ncl=2)
test.intersection.p(ventana1, ventana2, nsim, prop = 1, win = NULL,
   centroides1 = NULL, diametros1 = NULL, centroides2 = NULL, 
   diametros2 = NULL, ncl=2)
}
\arguments{
  \item{ventana1}{
A multiple-polygon window with the format \code{\link[spatstat]{owin}} of \pkg{spatstat}. The "accessory" species in the case of \code{test.intersection}.
}
  \item{ventana2}{
A multiple-polygon window with the format \code{\link[spatstat]{owin}} of \pkg{spatstat}. The species whose polygons will be rotated (i.e., the "focal" species) in the case of \code{test.intersection}.
}
  \item{nsim}{
Number of simulations for the Monte Carlo test.
}
  \item{prop}{
Proportion of the diameter of each polygon which will be employed to select potential close neighbours. See details.
}
  \item{win}{
Observation window, employed to control for edge effects. An objetct with the format \code{\link[spatstat]{owin}} of \pkg{spatstat}. If not provided, it will be estimated by the \emph{x} and \emph{y} ranges from \code{ventana}.
}
  \item{centroides1}{
Matrix or \code{data.frame} with the coordinates of the centroids of the polygons of \code{ventana1}.
}
  \item{diametros1}{
Vector with the diameters of the polygons of \code{ventana1}.
}
  \item{centroides2}{
Matrix or \code{data.frame} with the coordinates of the centroids of the polygons of \code{ventana2}.
}
  \item{diametros2}{
Vector with the diameters of the polygons of \code{ventana2}.
}
  \item{ncl}{
Number of clusters for the parallel implementation.
}
}
\details{
The summary statistic employed in the test is the overall sum of the overlapping areas of intersecting polygons. To reduce the computing burden, area intersections are only computed for polygon pairs for which the distance between their centoids is smaller than the sum of their largest diameters (weighted by the argument \code{prop}). To avoid edge effects, the test only considers polygons whose centroids are separated from the border of the observation window by a distance larger than their largest diameter.
In the pairwise tests, the polygons of \code{ventana1} are kept fixed in their original positions (i.e., the accessory species) and the polygons of \code{ventana2} (i.e., the focal species) are rotated.
}
\value{
Both \code{test.auto.intersection} and \code{test.intersection} return a vector of length \code{nsim}+1, with the sum of observed overlaping areas in the first position and subsequently with the sum of overlapping areas in each the simulated (i.e., randomly rotated) realizations of the null model.
}
\author{
Marcelino de la Cruz Rot
}

\examples{

data(Astragalus)
data(Sesleria)
# Test overlapping between Astragalus individuals
  # Ideally nsim should be at least 199
   Astragalus.test<- test.auto.intersection(Astragalus, nsim=19)
   # Observed overlapping area
   Astragalus.test[1]
   # p-value (negative value indicates that the observed overlapping is smaller 
   # than expected)
   pval(Astragalus.test)
   
# Test overlapping between Astragalus  and Sesleria individuals.
# Here, Sesleria is the accesory species (its individuals are kept fixed during the
# test) and Astragalus  the focal one (its individuals are rotated)
# Ideally nsim should be at least 199
   Sesleria.Astragalus.test<- test.intersection(ventana1= Sesleria, 
                                                   ventana2= Astragalus, nsim=19)
   # Observed overlapping area
   Sesleria.Astragalus.test[1]
   # p-value (negative value indicates that the observed overlapping is smaller
   # than expected)
   pval(Sesleria.Astragalus.test)
   
\donttest{ 
# Reducing computing burden when making repetitive testing
# First, put all the polygonal regions in a list, i.e.
 owins<- list(Astragalus, Sesleria)
 
 # compute diameters and centroids of the individual polygons 
 # in each polygonal region
 
centroids<- list()
diams<- list()
  for ( i in 1: length(owins)){
    cd<- centroidiam(owins[[i]])
    centroids[[i]] <- cd$centroids
    diams[[i]] <- cd$diams
}

# set the number of simulations for each test
nsim<-199

# create an array to store the results
result <- array(NA, dim=c(length(owins),length(owins),nsim+1))

t0<-Sys.time()
for ( i in 1: length(owins)){
   for ( j in 1: length(owins)){
     cat(i,j,"\n")
       if(j!=i) result[i,j,] <- test.intersection(owins[[i]], owins[[j]], nsim=nsim,
		centroides1=centroids[[i]], diametros1=diams[[i]],
		centroides2=centroids[[j]], diametros2=diams[[j]]) else
		result[i,j,] <-  test.auto.intersection(owins[[i]], nsim=nsim,
		centroides1=centroids[[i]], diametros1=diams[[i]])
   }
}
 Sys.time()-t0
 
 # observed values (focal species in columns)
 (observed<-  t(result[,,1]))
 # p-values
 tabla.p<- apply(result,c(1,2),pval)
 (p_values <- t(tabla.p))
 
 # Compare with parallelized versions:
  
# create an array to store the result.ps
result.p<- array(NA, dim=c(length(owins),length(owins),nsim+1))

t0<-Sys.time()
for ( i in 1: length(owins)){
   for ( j in 1: length(owins)){
     cat(i,j,"\n")
       if(j!=i) result.p[i,j,] <- test.intersection.p(owins[[i]], owins[[j]], nsim=nsim,
		centroides1=centroids[[i]], diametros1=diams[[i]],
		centroides2=centroids[[j]], diametros2=diams[[j]]) else
		result.p[i,j,] <-  test.auto.intersection.p(owins[[i]], nsim=nsim,
		centroides1=centroids[[i]], diametros1=diams[[i]])
   }
}
 Sys.time()-t0
 
 # observed values (focal species in columns)
 (observed.p<-  t(result.p[,,1]))
 # p-values
 tabla.p.p<- apply(result.p,c(1,2),pval)
 (p_values.p <- t(tabla.p.p))
 
 
}
 
}
\keyword{ spatial }
