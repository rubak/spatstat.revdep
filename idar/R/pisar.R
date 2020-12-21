
controldis <- function(d, m, mimark) {
         
	 #function to control and the distance object and
	# to reasure that it match the species in the ppp
	
       # TODO: control that d is a named/unnamed vector

        # meter opcion para que saque la media de los NA en el vector final.
	# (deberia ser un argumento tambien de pisar y de risar)

        
      if(is.null(d)){
	      d<- rep(1, ncol(m))
	      return(d)
      }
   
     #if(class(d)=="dist"){  # !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!   CAMBIADO 04/12/2019 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     if(inherits (d, what="dist")){ 
         d<- as.matrix(d)
     }
      
      if(is.matrix(d)){
         cuales <- colnames(d)%in%colnames(m)
	 if(sum(cuales)!=ncol(m)) stop("Names of species in the ppp and in the 'd'-distance matrix dont match")
	 if(!mimark%in%rownames(d)) stop("Name of the focal species not found  in the 'd'-distance matrix")
		 
	 # extract the distances to focal species	 
	 d<- d[mimark, cuales]
	 #match names ofspecies in cosamt and in the distance vector and re-order
	 names.ok<-match (colnames(m),names(d))
	 d <- d[names.ok]
       }
  
        return(d)
}

#---------------------------------------------------------------------------------------------------





pisar<-
function (mippp, mippp.sp = NULL, mimark = NULL, namesmark = NULL, d=NULL, 
    r = NULL, buffer = 0, bfw = NULL) 
{
    
    if (!is.null(namesmark)) 
        mippp$marks <- factor((mippp$marks[namesmark][[1]]))
    if (!is.null(mippp.sp) & !is.ppp(mippp.sp)) {
        mippp.sp <- NULL
        mimark <- mippp.sp
    }
    if (!is.null(mimark)) 
        if (mimark %in% levels(mippp$marks) == FALSE) {
            stop(paste(mimark, " can't be recognized as a mark\n\n\n     have you indicated in which column of thedataframe are the species\n      marks? (argument 'namesmark'\n\n"))
        }
    if (is.null(mippp.sp)) 
        mippp.sp <- mippp[mippp$marks == mimark]
    if (buffer != "adapt") {
        if (is.null(bfw)) 
            bfw <- erosion(mippp$window, buffer)
        mippp.sp <- mippp.sp[inside.owin(mippp.sp, w = bfw)]
        npoints <- rep(mippp.sp$n, length(r))
        names(npoints) <- r
    }
    cosamt <- mitable(mippp.sp, mippp, r)
    if (buffer == "adapt") {
        bdp <- bdist.points(mippp.sp)
        for (i in 1:length(r)) cosamt[[i]] <- cosamt[[i]][bdp >= 
            r[i], ]
        npoints <- sapply(cosamt, function(x) dim(x)[1])
    }
    
    #--------------------------------------------------------------------
    # vector de distancias de la focal al resto.
      d<- controldis(d, cosamt[[1]], mimark)
    #--------------------------------------------------------------------
    
    Ptj <- function(x) sum(x == 0)/length(x)
    pisar.r <- function(x) {
        result <- sum(d* apply(x, 2, function(x) 1 - Ptj(x))) # cambio:multiplicamos por el vector d
        return(result)
    }
    pisar <- sapply(cosamt, pisar.r)
    result <- data.frame(r = r, pisar = pisar)
    result <- fv(result, argu = "r", ylab = substitute(PISAR(r), 
        NULL), valu = "pisar", fmla = pisar ~ r, alim = c(min(r), 
        max(r)), labl = c("r", "%s(r)"), desc = c("radius of circle", 
        "%s"), fname = "PISAR")
    return(result)
}

#-----------------------------------------------------------------------------------------------------------------------------------------------

risar<- 
function (mippp, mippp.sp = NULL, mimark = NULL, namesmark = NULL, d=NULL, d0=NULL,
    r = NULL, buffer = 0, bfw = NULL) 
{
    
    if (!is.null(namesmark)) 
        mippp$marks <- factor((mippp$marks[namesmark][[1]]))
    if (!is.null(mippp.sp) & !is.ppp(mippp.sp)) {
        mippp.sp <- NULL
        mimark <- mippp.sp
    }
    if (!is.null(mimark)) 
        if (mimark %in% levels(mippp$marks) == FALSE) {
            stop(paste(mimark, " can't be recognized as a mark\n\n\n     have you indicated in which column of thedataframe are the species\n      marks? (argument 'namesmark'\n\n"))
        }
    if (is.null(mippp.sp)) 
        mippp.sp <- mippp[mippp$marks == mimark]
    if (buffer != "adapt") {
        if (is.null(bfw)) 
            bfw <- erosion(mippp$window, buffer)
        mippp.sp <- mippp.sp[inside.owin(mippp.sp, w = bfw)]
        npoints <- rep(mippp.sp$n, length(r))
        names(npoints) <- r
    }
    cosamt <- mitable(mippp.sp, mippp, r)
    if (buffer == "adapt") {
        bdp <- bdist.points(mippp.sp)
        for (i in 1:length(r)) cosamt[[i]] <- cosamt[[i]][bdp >= 
            r[i], ]
        npoints <- sapply(cosamt, function(x) dim(x)[1])
    }
    
    #--------------------------------------------------------------------
    # vector de distancias de la focal al resto.
    # La entrada "d" puede ser una matriz, una dist o un vector.
    # en todos los casos con "names" identicos (y al final en el mismo orden) que los de cosamt 
    # que salen de los levels de las marcas delppp
    # la salida para calcular risar es un vector de longitud igual al numero de especies (columnas) en cosamt
    
     
     d<- controldis(d=d, m=cosamt[[1]], mimark)
     d0<- controldis(d=d0, m=cosamt[[1]], mimark)
     
    
    #--------------------------------------------------------------------
    
    Ptj <- function(x) sum(x == 0)/length(x)
    pisar.r <- function(x) {
        result <- sum(d* apply(x, 2, function(x) 1 - Ptj(x))) # cambio:multiplicamos por el vector d
        return(result)
    }
    isar.r <- function(x) {
        result <- sum(d0* apply(x, 2, function(x) 1 - Ptj(x))) # cambio:multiplicamos por el vector d
        return(result)
    }
    pisar <- sapply(cosamt, pisar.r)
    isar <- sapply(cosamt, isar.r)
    risar<-pisar/isar
    result <- data.frame(r = r, risar = risar)
    result <- fv(result, argu = "r", ylab = substitute(rISAR(r), 
        NULL), valu = "risar", fmla = risar ~ r, alim = c(min(r), 
        max(r)), labl = c("r", "%s(r)"), desc = c("radius of circle", 
        "%s"), fname = "rISAR")
    return(result)
}

