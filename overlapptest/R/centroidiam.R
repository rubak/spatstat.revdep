centroidiam <-
function(ventana1){
     centroides1<- NULL
     diametros1 <- NULL
     for( i in 1:length(ventana1$bdry)){
           progressreport(i, length(ventana1$bdry))
           pol <-owin(poly=ventana1$bdry[[i]])
           centroides1 <- rbind(centroides1, unlist(centroid.owin(pol)))
           diametros1<- c(diametros1, diameter(pol))
        }
        return(list(diams=diametros1, centroids=centroides1))
}
