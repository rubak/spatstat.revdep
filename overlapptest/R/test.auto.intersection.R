test.auto.intersection<- function(ventana1, nsim=199, win=NULL,  
                prop=1, centroides1= NULL,    diametros1 =NULL){
    deg2rad <- function(deg) deg*pi/180
    
    result.obs<- rep(0, length(ventana1$bdry))
    result.sim <- matrix(0, nrow=length(ventana1$bdry), ncol=nsim)
    
    #primero , calculo de los centroides1 y diametros1 de los subpoligonos,
    # para calcular la interseccion solo con los que se encuentren a una distancia menor o igual que la suma de sus diametros1 mayores
    # (se podria afinar todavia mas y usar el radio de giracion)
    if( is.null(centroides1) |is.null(diametros1)){
       cat("computing centroids, etc \n\n")
      for( i in 1:length(ventana1$bdry)){
         progressreport(i, length(ventana1$bdry))
         pol <-owin(poly=ventana1$bdry[[i]])
        centroides1 <- rbind(centroides1, unlist(centroid.owin(pol)))
        diametros1<- c(diametros1, diameter(pol))
       }
    }
    
    
      # correccion borde (ventana1)
    if(is.null(win)) win<- owin(xrange=ventana1$xrange, yrange=ventana1$yrange)
    centroid.ppp<-ppp(centroides1[,1], centroides1[,2], window=win)
    bordeok<- which(bdist.points(centroid.ppp) > diametros1/2)
       # APlica la correccion
        #ventana1$bdry <- ventana1$bdry[bordeok]
        #centroides1<- centroides1[bordeok,]
	#diametros1<-diametros1[bordeok]

    
    
    
    # distancia entre los centroides1 de cada subpoligono
    centroidist<- as.matrix(dist(centroides1))
    # calculo de la interseccion entre los poligonos proximos
    
    cat("\n\n computing intersections, etc \n\n")
    #for( i in 1:length(ventana1$bdry)){
    for( i in bordeok){
      progressreport(i, length(ventana1$bdry))
      
      # suma de diametros1 del par de poligonos (focal-resto)
      #sumdiam<- diametros1+diametros1[i]
      sumdiam<- diametros1*prop+diametros1[i]*prop
      
      # Seleccion de poligonos que estan cerca del focal
      cualesok<- centroidist[i,] <= sumdiam
      #quitamos ademas el focal
      cualesok[i] <- FALSE
      ventanasin<- ventana1
      ventanasin$bdry <- ventanasin$bdry[cualesok]
      # definicion del poligono focal
      pol <-owin(poly=ventana1$bdry[[i]])
      
      # initerseccion observada
      intersecto.obs <- intersect.owin(pol, ventanasin, fatal = FALSE)
      #if(!is.null(intersecto.obs))  result.obs<- c(result.obs, area.owin(intersecto.obs)) else result.obs<- c(result.obs, 0)
      if(!is.null(intersecto.obs))  result.obs[i]<- area.owin(intersecto.obs) 
      
      
      
      # intersecciones simuladas bajo modelo nulo de rotacion
      for (sim in 1:nsim){
        pol<-rotate(pol, centre="centroid", angle=deg2rad(runif(1,
                                                                0,360)))
        intersecto.sim <- intersect.owin(pol, ventanasin, fatal = FALSE)
        if(!is.null(intersecto.sim))  result.sim[i,sim]<-
          area.owin(intersecto.sim)
      }
    }
    cat("\n\n")
    return(colSums(cbind(result.obs, result.sim)))
}
