# parallel implementation
#ncl=number of clusters
test.intersection.p<-
function (ventana1, ventana2, nsim, prop = 1, win = NULL, centroides1 = NULL, 
    diametros1 = NULL, centroides2 = NULL, diametros2 = NULL,
    ncl=2) 
{
    deg2rad <- function(deg) deg * pi/180
    intersectfunct<- function(x,v1,v2){ 
	       #focal<-rotawin(v2)
	       intersecto <- intersect.owin(v1, rotawin(v2),   fatal = FALSE)
               if (!is.null(intersecto)) return(area.owin(intersecto)) else return(0)
            }
    result.obs <- NULL
    result.sim <- matrix(0, nrow = length(ventana1$bdry), ncol = nsim)
    if (is.null(centroides1) | is.null(diametros1)) {
        cat("computing centroids, etc, for ventana1 \n\n")
        for (i in 1:length(ventana1$bdry)) {
            progressreport(i, length(ventana1$bdry))
            pol <- owin(poly = ventana1$bdry[[i]])
            centroides1 <- rbind(centroides1, unlist(centroid.owin(pol)))
            diametros1 <- c(diametros1, diameter(pol))
        }
    }
    if (is.null(centroides2) | is.null(diametros2)) {
        cat("computing centroids, etc, for ventana2 \n\n")
        for (i in 1:length(ventana2$bdry)) {
            progressreport(i, length(ventana2$bdry))
            pol <- owin(poly = ventana2$bdry[[i]])
            centroides2 <- rbind(centroides2, unlist(centroid.owin(pol)))
            diametros2 <- c(diametros2, diameter(pol))
        }
    }
    if (is.null(win)) 
        win <- owin(xrange = ventana1$xrange, yrange = ventana1$yrange)
    centroid1.ppp <- ppp(centroides1[, 1], centroides1[, 2], 
        window = win)
    bordeok <- which(bdist.points(centroid1.ppp) > diametros1/2)
    ventana1$bdry <- ventana1$bdry[bordeok]
    centroides1 <- centroides1[bordeok, ]
    diametros1 <- diametros1[bordeok]
    centroidist <- crossdist(centroides1[, 1], centroides1[, 
        2], centroides2[, 1], centroides2[, 2])
    sumdiam <- outer(diametros1 * prop, diametros2 * prop, "+")
    okpair <- sumdiam > centroidist
    ok1 <- rowSums(okpair) > 0
    cuales1 <- which(ok1)
    ok2 <- colSums(okpair) > 0
    cuales2 <- which(ok2)
    ventana1$bdry <- ventana1$bdry[cuales1]
    ventana2$bdry <- ventana2$bdry[cuales2]
    intersecto <- intersect.owin(ventana1, ventana2, fatal = FALSE)
    if (!is.null(intersecto)) 
        obs <- area.owin(intersecto)
    else obs <- 0
    if (length(cuales1) < 1 | length(cuales2) < 1) 
        simu <- rep(0, nsim)
    else {
        #require(parallel)
	 # define environment created by the function
	    e<-environment()
        cl <- makeCluster(ncl)
	   clusterExport(cl, c( "ventana1", "ventana2",  "intersectfunct"), envir=e)
	   clusterExport(cl, c("rotawin", "intersect.owin", "area.owin"))
	   simu<-parSapplyLB(cl, 1:nsim, intersectfunct, v1=ventana1, v2=ventana2)
	stopCluster(cl)
	
   }
   result <- c(obs, simu)
    return(result)
}
