# es el que se llamaba ipsim4CORR.R


 ipsim<- function(pp, mimark,sigma=0, lambda=NULL, namesmark=NULL){
    # check that in case of multimarked (data.frame) pp, namesmark is provided
    # if mippp has a data,frame of marks, ensure that the name of the column with the labels for species is set in argument "namesmark"
    dimarks <- dim(pp$marks)
    if(!is.null(dimarks) & is.null(namesmark)) stop(paste("For the simulations you should indicate to ipsim argument 'namesmark' which of these ('",paste(names (pp$marks), collapse="' or '"),"') wear the labels for species in the community pp", sep=""))
     if(!is.null(namesmark)){
             nmarks<- names(pp$marks)
            u <- split(pp, f=namesmark)
            # randomize (IPP)) the point pattern of the label
            
            if(!is.null(lambda)){
                  um <-  rpoispp(lambda)
            }
            if(is.null(lambda)){
                if (sigma==0){
                     lambda=u[[mimark]]$n/area.owin(pp$window)# agnadido 29041017 $window
                      um <-  rpoispp(lambda, win=pp$window)
                 }

                if (sigma>0){
                     lambda <- density.ppp(unmark(u[[mimark]]), sigma=sigma)
                     um <-  rpoispp(lambda)
                 }
                     
            } 
            

            um$marks=data.frame(matrix(NA, nrow=um$n,ncol= length(nmarks)))
            names(um$marks) <- nmarks
            um$marks[,namesmark] <- rep(mimark, um$n)
            u[[mimark]] <- um
            split(pp, f=namesmark) <- u
            namesmarks <-names(pp$marks)
      }
     if(is.null(namesmark)){
         # split the multivariate pp
          u <- split(pp)
          # randomize (IPP)) the point pattern of the label

            if(!is.null(lambda)){
                  um <-  rpoispp(lambda)
            }

          if(is.null(lambda)){
                if (sigma==0){
                     lambda=u[[mimark]]$n/area.owin(pp$window) # agnadido 29041017 $window
                      um <-  rpoispp(lambda, win=pp$window)
                 }

                if (sigma>0){
                     lambda <- density.ppp(unmark(u[[mimark]]), sigma=sigma)
                     um <-  rpoispp(lambda)
                 }
                     
            } 

          u[[mimark]] <- um

          # recompose back the splited as a multivariate pp
          split(pp) <- u
      }
      return(pp)
 }