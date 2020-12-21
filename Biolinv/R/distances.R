#' Hidden functions for EM algorithm.
#'
#' These functions are used by the EM() function and were contained in the file 'distances.R'.
#'
#' @keywords internal
#'
#' @importFrom fields rdist
#'
#' @author Beatrix Jones, Luca Butikofer
#'
#' @export

distances<-function(dati){ #dati should have (at least) 3 columns; column 2 should be lattitude and 3 longitude; one column should be called "year"
  #print("inside distances")
  years<- unique(sort(dati$year))
  distance.pastYrs<-rep(NA, dim(dati)[[1]])
  distance.currentYr<-rep(NA, dim(dati)[[1]])
  for( i in 2:length(years)){
    if(sum(dati$year==years[i])>0 & sum(dati$year< years[i])>0){
      past <- subset(dati, dati$year<years[i])
      pres<-  subset(dati, dati$year==years[i])
      dist <- 0.001*rdist(pres[,3:2],past[,3:2])
      dist.min<-apply(dist,1,min)
      distance.pastYrs[dati$year==years[i]]<-dist.min
      if(sum(dati$year==years[i])>1){
        dist<-0.001*rdist(pres[,3:2], pres[,3:2])
        diag(dist)<-NA
        dist.min<-apply(dist, 1, min, na.rm=T)
        distance.currentYr[dati$year==years[i]]<-dist.min}
    }}
  return(data.frame(distance.pastYrs, distance.currentYr, year=dati$year))
}

distance.by.year<-function( obs.dist) {
  #print("inside distance by year")
  #obs.dist should be the output of distances
  year1<-min(obs.dist$year)
  distances<-apply(obs.dist[obs.dist$year>year1,1:2],1,min, na.rm=T)  #shortest distance in current or past year.
  distances<-c(rep(NA, sum(obs.dist$year==year1)),distances)

  output<-split(distances,obs.dist$year)   #output will be a list with entry for each year
  #print(output)

  past.dist<-split(obs.dist$distance.pastYrs, obs.dist$year)  #same list for past years


  for(i in 2:length(output)){
    if(sum(past.dist[[i]]==output[[i]])==0){  #if the closest  point is never a point from past years
      index=which.min(past.dist[[i]])           #pick the closest past year distance point
      output[[i]][index]<-past.dist[[i]][index]  #replace that distance with the longer, past year distance
    }
  }
  weights<-list()
  for(i in 1:length(output)){
    weights[[i]]<-rep(0.5, length(output[[i]]))}
  output[[1]]<-rep(NA, length(output[[1]]))
  return(list(output,weights))
}

random.distances<-function(dati, rdm10000){
  min.dist <- list()
  kern<-list()
  years<-unique(sort(dati$year))
  for (i in 1:length(years)){
    pres <- subset(dati,dati$year< years[i]+1)
    dist <- 0.001*rdist(rdm10000[,3:2],pres[,3:2])
    out.min <- apply(dist,1,min)
    min.dist[[i]]<-out.min
    kern[[i]]<-density(min.dist[[i]],from=0.001, to=max(min.dist[[i]])+5, n=2^13)
    bar.width<-kern[[i]]$x[2]-kern[[i]]$x[1]
    nc<-sum(kern[[i]]$y*bar.width)
    kern[[i]]$y<-kern[[i]]$y/nc
  }

  return(kern)
}

g.x<-function(x,Kern){    #Kern should be the kernel pertaining to just that year.

  a<-rep(NA, length(x))
  for(i in 1:length(x)){
    a[i]<-which.min(abs(Kern$x-x[i]))}
  b<-Kern$y[a]
  return(b)
}

ddispersq<-function(x,alpha,c){

  answer<-c/(alpha*gamma(1/c))*exp(-1*(x/alpha)^c)
  answer[x<=0]<-0
  return(answer);
}

f2.x<-function(x, sigma, year.index, xrange){  #xrange should be the xvalues used in the kernel for that year
  Kern<-ddispersq(xrange, sigma*sqrt(2),2)
  bar.width<-xrange[2]-xrange[1]
  nc<-sum(Kern*bar.width)
  Kern<-Kern/nc
  a<-rep(NA, length(x))
  for( i in 1:length(x)){
    a[i]<-which.min(abs(xrange-x[i]))}
  b<-Kern[a]
  return(b)
}

EM.update<-function(output,Weights,Sigma,Pi, Kerns){
  new.weights<-Weights
  new.sigma2<-0
  n<-sapply(output,length)
  for(i in 2:length(n)){
    if(n[i]>0){
      for(j in 1:n[i]){
        new.weights[[i]][j]<-Pi*f2.x(output[[i]][j],Sigma,i, Kerns[[i]]$x)/(Pi*f2.x(output[[i]][j],Sigma,i,Kerns[[i]]$x)+(1-Pi)*g.x(output[[i]][j],Kerns[[i]]))
      }
      #print(c(n[i],length(new.weights[[i]]), length(output[[i]])))
      new.sigma2<-new.sigma2+sum(new.weights[[i]]*output[[i]]^2)
    }}
  new.sigma2<-new.sigma2/sum(sapply(new.weights,sum))
  new.sigma2<-new.sigma2^(1/2)

  new.pi<-sum(sapply(new.weights, sum))/sum(n)
  return(list(new.weights, new.sigma2, new.pi))
}

EM.max<-function(output, Weights, Sigma, Pi,  Kerns){
  EM.update(output, Weights,Sigma,Pi,Kerns)->updateA
  EM.update(output, updateA[[1]], updateA[[2]],updateA[[3]],Kerns)->updateB

  while(abs(updateA[[3]]-updateB[[3]]) >0.00001){
    updateA=updateB
    updateB=EM.update(output, updateA[[1]], updateA[[2]], updateA[[3]], Kerns)
    #print(c(updateB[[3]], updateB[[2]]))
  }
  return(updateB)
}
