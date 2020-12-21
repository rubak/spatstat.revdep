# con ncl=3, reduce a la mitad el tiempo del serial con 199 simulaciones para Sesleria

test.auto.intersection.p <-
function (ventana1, nsim = 199, win = NULL, prop = 1, centroides1 = NULL, 
    diametros1 = NULL, ncl=2) 
{
    deg2rad <- function(deg) deg * pi/180
    intersect.autofunct<- function(x, p, vs){
	   p <- rotate(p, centre = "centroid",
  	                      angle = deg2rad(runif(1,  0, 360)))
	   intersecto.sim <- intersect.owin(p, vs, fatal = FALSE)
	   if (is.null(intersecto.sim)) return(0) else return(area.owin(intersecto.sim))
        }
	
        todo<- function(i, diametros1,centroidist, ventana1,nsim){
             sumdiam <- diametros1 * prop + diametros1[i] * prop
            cualesok <- centroidist[i, ] <= sumdiam
            cualesok[i] <- FALSE
            ventanasin <- ventana1
            ventanasin$bdry <- ventanasin$bdry[cualesok]
            pol <- owin(poly = ventana1$bdry[[i]])
            intersecto.obs <- intersect.owin(pol, ventanasin, fatal = FALSE)
	    obs<-0
            if (!is.null(intersecto.obs)) obs<- area.owin(intersecto.obs) 
	    res<-c(obs, sapply(1:nsim, intersect.autofunct, p=pol, vs=ventanasin))
	    return(res)
        }
	
	
	
    result.obs <- rep(0, length(ventana1$bdry))
    result.sim <- matrix(0, nrow = length(ventana1$bdry), ncol = nsim)
    if (is.null(centroides1) | is.null(diametros1)) {
        cat("computing centroids, etc \n\n")
        for (i in 1:length(ventana1$bdry)) {
            progressreport(i, length(ventana1$bdry))
            pol <- owin(poly = ventana1$bdry[[i]])
            centroides1 <- rbind(centroides1, unlist(centroid.owin(pol)))
            diametros1 <- c(diametros1, diameter(pol))
        }
    }
    if (is.null(win)) 
        win <- owin(xrange = ventana1$xrange, yrange = ventana1$yrange)
    centroid.ppp <- ppp(centroides1[, 1], centroides1[, 2], window = win)
    bordeok <- which(bdist.points(centroid.ppp) > diametros1/2)
    centroidist <- as.matrix(dist(centroides1))
    cat("\n\n computing intersections, etc \n\n")
     #require(parallel)
        e<-environment()
        cl <- makeCluster(ncl)
	   clusterExport(cl, c("todo","deg2rad","intersect.autofunct"), envir=e) 
	   clusterExport(cl, c("rotate", "intersect.owin", "area.owin","owin"))
           result<-parSapplyLB(cl, bordeok, todo,diametros1=diametros1,
                               centroidist=centroidist, ventana1=ventana1,nsim=nsim)
           stopCluster(cl)
    
    cat("\n\n")
    #print(dim(result))
    return(rowSums(result)) 
}