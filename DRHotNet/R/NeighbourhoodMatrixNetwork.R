#' Creates the neighbourhood structure of a linear network
#' 
#' Given a linear network structure, this function creates the neighbourhood matrix ("queen" criterion) associated to it. Two segments of the network are neighbours if they share a vertex   
#' 
#' @param network - A \code{linnet} object representing a linear network structure
#' @return Returns a \code{listw} object in \code{"W"} style
#' @examples 
#' library(DRHotNet)
#' library(spatstat.core)
#' library(spdep)
#' library(raster)
#' library(maptools)
#' chicago_neighbourhood <- NeighbourhoodMatrixNetwork(chicago$domain)
#' class(chicago_neighbourhood)
#' chicago_neighbourhood$neighbours[[1]]
#' @export
NeighbourhoodMatrixNetwork <- function(network){
  aux=SpatialLines2PolySet(as.SpatialLines.psp(as.psp(network)))
  aux=PolySet2SpatialPolygons(aux)
  queen=poly2nb(aux, queen=TRUE)
  W=nb2listw(queen, style="W", zero.policy=TRUE)
  return(W)
}

