K1K2<-
function (X, i, j, nsim = 99, nrank = 1, r = NULL, correction = "isotropic") 
{
    marx <- marks(X)
    I <- (marx == i)
    J <- (marx == j)
    k12 = Kmulti.ls(X, I, J, r, corre = correction)
    r = k12$r
    k1 = Kest(split(X)[names(split(X)) == i][[1]], r = r, correction = correction)
    k2 = Kest(split(X)[names(split(X)) == j][[1]], r = r, correction = correction)
    k1k2.o = k1[[3]] - k2[[3]]
    k112 = k1[[3]] - k12[[3]]
    k212 = k2[[3]] - k12[[3]]
    k1k2.s = NULL
    k112.s = NULL
    k212.s = NULL
    for (n in 1:nsim) {
        progressreport(n, nsim)
        X$marks[I | J] = sample(X$marks[I | J])
        marx <- marks(X)
        I <- (marx == i)
        J <- (marx == j)
        k12.s = Kmulti.ls(X, I, J, r = r, corre = correction)
        k1.s = Kest(split(X)[names(split(X)) == i][[1]], r = r, 
            correction = correction)
        k2.s = Kest(split(X)[names(split(X)) == j][[1]], r = r, 
            correction = correction)
        k1k2.s = cbind(k1k2.s, k1.s[[3]] - k2.s[[3]])
        k112.s = cbind(k112.s, k1.s[[3]] - k12.s[[3]])
        k212.s = cbind(k212.s, k2.s[[3]] - k12.s[[3]])
    }
    orderstat <- function(x, n) sort(x)[n]
    k1k2.lo <- apply(k1k2.s, 1, orderstat, n = nrank)
    k1k2.hi <- apply(k1k2.s, 1, orderstat, n = nsim - nrank + 
        1)
    k112.lo <- apply(k112.s, 1, orderstat, n = nrank)
    k112.hi <- apply(k112.s, 1, orderstat, n = nsim - nrank + 
        1)
    k212.lo <- apply(k212.s, 1, orderstat, n = nrank)
    k212.hi <- apply(k212.s, 1, orderstat, n = nsim - nrank + 
        1)
    k1k2 = k12
    k1k2[[2]] = k1k2.hi
    k1k2[[3]] = k1k2.o
    k1k2[[4]] = k1k2.lo
   
   attributes(k1k2)$fname <- "D" ## function difference-------------------
    attributes(k1k2)$labl = c( "r", "hat(%s)[hi](r)", "hat(%s)[obs](r)" ,"hat(%s)[lo](r)" ) ## ------------------- 
    attributes(k1k2)$ylab = "D(r) = K[1](r) - K[2](r)"                   ## ------------------- 
    attributes(k1k2)$yexp = expression('D(r) = ' * K[1](r) - K[2](r))## ------------------- 
    
    attributes(k1k2)$names = c("r", "hi",  "D", "lo")
    attributes(k1k2)$desc = c("distance argument r" , "upper pointwise envelope of %s  from simulations", 
         "observed value of %s for data pattern", 
        "lower pointwise envelope of %s  from simulations") ## ------------------- 
	
	attributes(k1k2)$dotnames <- c("hi", "D", "lo") ## ------------------- 
       #attributes(k1k2)$shade <- c("lo", "hi") ## ------------------- 

	
   
    k1k12 = k12
    k1k12[[2]] = k112.hi
    k1k12[[3]] = k112
    k1k12[[4]] = k112.lo
  
     attributes(k1k12)$fname <- "D" ## function difference-------------------
    attributes(k1k12)$labl = c( "r", "hat(%s)[hi](r)", "hat(%s)[obs](r)" ,"hat(%s)[lo](r)" ) ## ------------------- 
    attributes(k1k12)$ylab = "D(r) = K[1](r) - K[12](r)"                   ## ------------------- 
    attributes(k1k12)$yexp = expression('D(r) = ' * K[1]*'(r)' - K[12]^"*" *'(r)')## ------------------- 
    
    attributes(k1k12)$names = c("r", "hi",  "D", "lo") ## ------------------- 
    attributes(k1k12)$desc = c("distance argument r" , "upper pointwise envelope of %s  from simulations", 
         "observed value of %s for data pattern", 
        "lower pointwise envelope of %s  from simulations") ## ------------------- 
	
	attributes(k1k12)$dotnames <- c("hi", "D", "lo") ## ------------------- 
       #attributes(k1k12)$shade <- c("lo", "hi") ## ------------------- 


    k2k12 = k12
    k2k12[[2]] = k212.hi
    k2k12[[3]] = k212
    k2k12[[4]] = k212.lo
    
    attributes(k2k12)$fname <- "D" ## function difference-------------------
    attributes(k2k12)$labl = c( "r", "hat(%s)[hi](r)", "hat(%s)[obs](r)" ,"hat(%s)[lo](r)" ) ## ------------------- 
    attributes(k2k12)$ylab = "D(r) = K[2](r) - K[12](r)"                   ## ------------------- 
    attributes(k2k12)$yexp = expression('D(r) = ' *K[2]*'(r)' - K[12]^"*" *'(r)')## ------------------- 
    
    attributes(k2k12)$names = c("r", "hi",  "D", "lo") ## ------------------- 
    attributes(k2k12)$desc = c("distance argument r" , "upper pointwise envelope of %s  from simulations", 
         "observed value of %s for data pattern", 
        "lower pointwise envelope of %s  from simulations") ## ------------------- 
	
	attributes(k2k12)$dotnames <- c("hi", "D", "lo") ## ------------------- 
       #attributes(k2k12)$shade <- c("lo", "hi") ## ------------------- 
    return(list(k1k2 = k1k2, k1k12 = k1k12, k2k12 = k2k12))
}
