check.ventana <-
function(ventana){
    malos1 <- NULL
    malos2<- NULL
    for( i in 1:length(ventana$bdry)){
       test <- try(owin(poly=ventana$bdry[[i]]), silent=TRUE)
       if(class(test)=="try-error"){
            malos1<- c(malos1, i) # anota en que subpoligono hay un error
            # intenta darle la vuelta a las coordenadas a ver si eso lo resuleve
            test2 <- try( owin(poly=lapply(ventana$bdry[[i]], rev), silent=TRUE))
              if(class(test2)!="try-error"){
                   # si al dar la vuelta a las cooredenadas se arregla, grabamos las nuevas coordenadas
                   ventana$bdry[[i]] <-lapply(ventana$bdry[[i]], rev)
	        } else malos2<- c(malos2, i)
	}
     }
   attr(ventana, "corrected")  <-malos1
   attr(ventana, "not.corrected") <- malos2
   if(!is.null(malos1)) cat( length(malos1), "problematic polygon(s) detected \n \n")
   if(!is.null(malos1) & is.null(malos2)) cat("all problematic polygons have been repared\n \n")
   if(!is.null(malos2)) cat(length(malos2), "problematic polygon(s) could not been repared\n \n")
   return(ventana)
}
