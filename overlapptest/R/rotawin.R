rotawin <-
function(ventana){
   deg2rad <- function(deg) deg*pi/180
   for( i in 1:length(ventana$bdry)){
      ventana$bdry[[i]] <- rotate(owin(poly=ventana$bdry[[i]]), centre="centroid", angle=deg2rad(runif(1, 0,360)))$bdry[[1]]
   }
   return(ventana)
}
