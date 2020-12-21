 ipsimlist<- function(pp, mimark,listsim){
   # inhomogeneous simulation of a type within a
   # multivariate point pattern
   # listsim: list with simulated pp from simulador2

   listsim <- lapply(listsim, function(x, pp2=pp, mimark2=mimark)   {
                    # First split the multivariate pp
                     u <- split(pp2)
	             # Second, put the fit inhomogenous  simulated IPP
                      u[[mimark2]] <- unmark(x)
	         # recompose back the splited as a multivariate pp
                 split(pp2) <- u
                 return(pp2)
                 }
	        )
   return(listsim)
 }
 
 
 
 
 
