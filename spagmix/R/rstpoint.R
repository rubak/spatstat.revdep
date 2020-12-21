rstpoint <- function(n, f, W = NULL, correction = 1.5, maxpass = 50) {
  if(!inherits(f,"stim")&&!inherits(f,"stden")) stop("'f' must be of class 'stim' or 'stden'")

  if(inherits(f,"stden")){
    g <- list(v=f$z,tlay=f$tgrid,a=abind(lapply(f$z,as.matrix),along=3))
    f <- g
  }

  if(is.null(W)) W <- as.polygonal(Window(f$v[[1]]))
  if(!is.owin(W)) stop("'W' must be of spatstat class 'owin'")
  if(is.null(intersect.owin(Window(f$v[[1]]),W,fatal=FALSE))) stop("'W' must overlap the spatial domain of 'f'")

  nacc <- 0
  dummyim <- f$v[[1]]
  dummyim[] <- 1
  xacc <- yacc <- tacc <- c()
  ngen <- ceiling(n*correction)
  tlay <- f$tlay
  tstep <- tlay[2]-tlay[1]
  tgrid <- seq(min(tlay)-0.5*tstep,max(tlay)+0.5*tstep,length=length(tlay)+1)
  pass <- 0
  while(nacc<n&&pass<maxpass) {
    cand <- rimpoly(ngen,dummyim,W)
    cand.t.raw <- runif(ngen,min(tgrid),max(tgrid)) #Randomly generated times
    cand.t.int <- findInterval(cand.t.raw,tgrid,all.inside=TRUE) #Interval number for each time

    cand.dens <- rep(NA,ngen)
    ut <- unique(cand.t.int)
    for(i in 1:length(ut)){
      cand.ti <- which(cand.t.int==ut[i])
      cand.dens[cand.ti] <- safelookup(f$v[[ut[i]]],cand[cand.ti])
    }
    cand.probs <- cand.dens/max(f$a,na.rm=TRUE)
    cand.acc <- which(runif(ngen)<cand.probs) #Indices of the accepted points

    xacc <- c(xacc,cand$x[cand.acc])
    yacc <- c(yacc,cand$y[cand.acc])
    tacc <- c(tacc,cand.t.raw[cand.acc])
    nacc <- nacc + length(cand.acc)
    pass <- pass + 1
  }

  if(nacc<n) warning(paste("'maxpass' reached with only",nacc,"points accepted"))

  result <- ppp(c(),c(),window=W)
  if(nacc>0){
    accs <- 1:min(n,nacc)
    result <- ppp(x=xacc[accs],yacc[accs],marks=tacc[accs],window=W)
  }

  return(result)
}
