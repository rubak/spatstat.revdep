# function to compute percentage acumulators, neutralas and repellers from a list of envelope objects to idar functions
proportion.idar <- function(envlist, alfa = 0.05){

   # envlist: list of envelope objects (one for each focal species)

   allsim.ok <- sapply(envlist, function(x) !is.null(attributes(x)$simfuns))
   # number of usefull species (i.e. with simulated functions)
   
   nsp <-sum(allsim.ok)
   if( nsp <length(envlist)) warning("some of the envelope objects don't wear the simulated functions")

# use  only envelope objects with simulated functions
   envlist <- envlist[allsim.ok]

  nsims<- sapply(envlist, function(x) dim(attributes(x)$simfuns)[2]-1)
  nsim<- nsims[1] #TODO: Check that all envelope objects have the same numbre of simulations !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  r <- envlist[[1]]$r  #TODO: Check that all envelope objects have the same  sequence of r values; use only common r values  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   # get rank of observed function among the simulated functions at each scale r
   scores <- sapply(envlist,
                      function(sp.env) {apply(data.frame(sp.env$obs, data.frame(attributes(sp.env)$simfuns)[,-1]), 1, function(x) rank(x)[1])}
                )
 
      quant <- c(0,1)+c(alfa/2,-alfa/2)
      dataquants <- quantile(0:(nsim+1), probs=quant)

       accumulators <- apply(scores, 2, function(x) x>=dataquants[2])
       repellers    <- apply(scores, 2,  function(x) x<=dataquants[1])
       neutrals         <- apply( scores, 2,  function(x) x>dataquants[1] & x<dataquants[2])

       behaviour <- as.data.frame(neutrals)
       behaviour[accumulators]<- "A"
       behaviour[neutrals]<- "."
       behaviour[repellers]<- "R"
    

       
       percentage <- data.frame(
           p.accumulators = apply(accumulators,1,sum)*100/nsp,
           p.repellers = apply(repellers,1,sum)*100/nsp,
           p.neutrals = apply(neutrals,1,sum)*100/nsp
            )
       
       result <- list(percentage=percentage, nsp=nsp, nsim=nsim, alfa=alfa, r=r,  behaviour= behaviour)
       class(result) = c("pidar", class(result))
       return(result)
       
}

plot.pidar<-function(x,cols=c(1,2,3),type=c("l","o","o"), pch=c(NA,19,19),
                            lty=c(1,1,1),  legend=TRUE, p.legend="topleft",...){

   plot(x$r, x$percentage$p.neutrals, type=type[1], pch=pch[1],col=cols[1], lty=lty[1],
                 ylim=c(0,100), xlab="Radius r of circle (m)", ylab="Percentage of cases",...)
   lines(x$r, x$percentage$p.accumulators, type=type[2], pch=pch[2],col=cols[2],lty=lty[2],...)
   lines(x$r, x$percentage$p.repellers, type=type[3], pch=pch[3], col=cols[3],lty=lty[3],...)
   
   if(legend==TRUE){
        if(length(p.legend)>1){
           xleg<- p.legend[1]
           yleg<- p.legend[2]
           legend(x=xleg,y=yleg, legend=c("neutrals", "accumulators","repellers"),
                                           col=cols, pch=pch, lty=lty)
        }
         if(length(p.legend)==1){
             legend(x=p.legend, legend=c("neutrals", "accumulators","repellers"),
                                           col=cols, pch=pch, lty=lty)
        }
   }
}


