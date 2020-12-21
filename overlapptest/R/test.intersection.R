test.intersection <-
function(ventana1, ventana2, nsim, prop=1, win=NULL, 
                              centroides1=NULL, diametros1=NULL, centroides2=NULL, diametros2=NULL){

    deg2rad <- function(deg) deg*pi/180
    
    result.obs<- NULL
    result.sim <- matrix(0, nrow=length(ventana1$bdry), ncol=nsim)
    
    #primero , calculo de los centroides y diametros de los subpoligonos,
    # para calcular la interseccion solo con los que se encuentren a una distancia menor o igual que la suma de sus diametros mayores
    # (se podria afinar todavia mas y usar el radio de giracion)
    
    if( is.null(centroides1) |is.null(diametros1)){
        cat("computing centroids, etc, for ventana1 \n\n")
        for( i in 1:length(ventana1$bdry)){
           progressreport(i, length(ventana1$bdry))
           pol <-owin(poly=ventana1$bdry[[i]])
           centroides1 <- rbind(centroides1, unlist(centroid.owin(pol)))
           diametros1<- c(diametros1, diameter(pol))
        }
    }
    if( is.null(centroides2) |is.null(diametros2)){
        cat("computing centroids, etc, for ventana2 \n\n")
        for( i in 1:length(ventana2$bdry)){
            progressreport(i, length(ventana2$bdry))
            pol <-owin(poly=ventana2$bdry[[i]])
            centroides2 <- rbind(centroides2, unlist(centroid.owin(pol)))
            diametros2<- c(diametros2, diameter(pol))
        }
   }
    
    # correccion borde (ventana1)
    if(is.null(win)) win<- owin(xrange=ventana1$xrange, yrange=ventana1$yrange)
    centroid1.ppp<-ppp(centroides1[,1], centroides1[,2], window=win)
    bordeok<- which(bdist.points(centroid1.ppp) > diametros1/2)
       # APlica la correccion
        ventana1$bdry <- ventana1$bdry[bordeok]
        centroides1<- centroides1[bordeok,]
	diametros1<-diametros1[bordeok]
	
    # distancia entre los centroides de la ventana 1 (filas) con los centroides de la ventana2 (columnas)
     centroidist<-crossdist(centroides1[,1], centroides1[,2],centroides2[,1],centroides2[,2])
     # suma de diametros ( o porporcion de diametros) de cada poligono de veentana1 con cadapoligono de ventana2
     sumdiam<- outer(diametros1*prop, diametros2*prop, "+")
     
     # parejas a tener en cuenta
     okpair<-sumdiam>centroidist
     # individuos de ventana1 que habria que rotar
     ok1<- rowSums(okpair)>0
     cuales1<-which(ok1)
     
     # individuos de ventana2 que habria que tener en cuenta
     ok2<- colSums(okpair)>0
     cuales2<-which(ok2)
     
     ventana1$bdry <- ventana1$bdry[cuales1]
     ventana2$bdry <- ventana2$bdry[cuales2]
     
intersecto <- intersect.owin(ventana1, ventana2, fatal = FALSE)
   if(!is.null(intersecto))  obs<- area.owin(intersecto) else obs<- 0   
   
   #if there are not polygons left after controling for distances, put all at 0
   if( length(cuales1)<1 | length(cuales2)<1) simu<- rep(0, nsim)   else {   
      simu<- NULL
      for ( i in 1:nsim){
         progressreport(i,nsim)
         intersecto <- intersect.owin( ventana1, rotawin(ventana2), fatal =FALSE)
         if(!is.null(intersecto))  simu<- c(simu, area.owin(intersecto)) else simu<- c(simu,0)
      }
   }
   result<- c(obs, simu)
   #class(result) <- c("intest", class(result))
   return(result)
}
     