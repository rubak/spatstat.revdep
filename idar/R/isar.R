###############################################################################
#
#    ISAR.R
#
#   Individual Species Area Relationship
#
#  $Version: 0.1 $  $Date: 2013/02/17 23:48:03 $
#
#   function to compute Indivdual Species Area Relationships (i.e., species
#   accumulation curves around individual species) in a community.
#
#    Value: an object of class isar, bassicaly a list with the following
#    elemnts:
#     r: vector of radii at which ISAR(r) has been estimated
#     isar: vector with the ISAR(r) values
#     npoints: number of focal points employed to estimate isar(r) at each 
#              radius r
#
#    Arguments:
#     mippp:    Multivariate point pattern (ppp of spatstat)
#     mippp.sp: Focal (unmarked) point pattern (ppp of spatstat) or character
#                string indicating one of the marks in the multivariate mippp 
#     mimark:   character string indicating one of the marks in the
#                multivariate mippp
#     namesark: name of the column with species names in the data.frame of
#                marks (only for multimarked ppp's)
#     r:        vector with the sequence of radii (>0) at which estimate
#                ISAR(r)
#     buffer:   a number indicating the width around the window of the focal 
#                point pattern that will be excluded from the computations of
#                ISAR(r) to control de edge effect, or the string "adapt", 
#                indicating that an adaptive border (of width ri) will be 
#                excluded for the computation of each ISAR(ri)
#     bfw:      an owin object indicating a subset of the focal point pattern 
#                that will be employed to compute ISAR(r)
#
#   Details:
#     btw can be employed to select different habitats and compute ISAR((R) i 
#      each habitat

 




 
isar<- function(mippp, mippp.sp=NULL, mimark=NULL,  namesmark=NULL, r=NULL,
	buffer=0, bfw=NULL) {


  #If species marks are within a data.frame of marks, get this column and 
  # discard the rest
  if(!is.null(namesmark)) mippp$marks <- factor((mippp$marks[namesmark][[1]]))  # OJO: que representa el [[1]] del final: es para coger los valores. si no se pone da problemas !!!!!!!!!!!!!!!!!!!!!!!!!!!
  #!!!!!!!!!!!!!!!!!!
  
  
  #If mippp.sp has not been provided as a ppp, build it from the multivariate
  # mippp and the provided mark
  if(!is.null(mippp.sp) & !is.ppp(mippp.sp)) {
   mippp.sp <- NULL
   mimark <- mippp.sp
  }

  if(!is.null(mimark)) if(mimark%in%levels(mippp$marks)==FALSE){
     stop(paste(mimark, " can't be recognized as a mark\n\n
     have you indicated in which column of thedataframe are the species
      marks? (argument 'namesmark'\n\n"))
      }
  if (is.null(mippp.sp)) mippp.sp<- mippp[mippp$marks==mimark]
       
  #If indicated, correct edge effect with a fixed buffer of size "buffer"
  # or with the window "bfw"
  if(buffer!="adapt"){	
    if(is.null(bfw)) bfw<- erosion( mippp$window, buffer)
    mippp.sp <- mippp.sp[inside.owin(mippp.sp, w=bfw)]
    npoints <- rep(mippp.sp$n, length(r))
    names(npoints) <- r
  }
       

  #Tabulate the number of individuals of each species (i.e. types) around each 
  # individual of the focal species within circles of each of the r radii
  
  
  cosamt <- mitable(mippp.sp, mippp, r)
	
  #If indicated correct edge effect with an adaptive buffer
  if(buffer=="adapt"){
    bdp <- bdist.points(mippp.sp)
    for (i in 1:length(r)) cosamt[[i]]<-cosamt[[i]][bdp>=r[i],]
    npoints <- sapply(cosamt, function(x) dim(x)[1])
  }
  	
  #Define bivariate emptines probability Ptj(0,r)
  # i.e., the probability of not finding species j within circles of radius r 
  # around the individuals of focal species.
  # x will be the vector of counts of specis j within circles of radius r 
  # around each individual of the focal species.
  Ptj <- function(x) sum(x==0)/length(x)

  #Define function to compute ISAR
  isar.r <- function(x){
    result <-sum(apply(x,2,function(x) 1-Ptj(x)))
    return(result)
  }
  
  #Compute ISAR(r) for each radius r
  isar <- sapply(cosamt, isar.r)

  # return results
  # result <- list(r=r, isar=isar, npoints=npoints)
   result <- data.frame(r=r, isar=isar)
  
   result <- fv(result, argu="r", ylab=substitute(ISAR(r), NULL),
              valu="isar", fmla=isar ~ r,
              alim=c(min(r),max(r)),
              labl=c("r", "%s(r)"),
              desc=c("radius of circle",
                     "%s"),
              fname="ISAR")

  
  
  
  #class(result)<-c("isar", class(result))
  return(result)
}
###############################################################################