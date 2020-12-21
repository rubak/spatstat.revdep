###############################################################################
#
#    IPSRAR.R
#
#   Individual Phylogenetic Species Richness -Area Relationship
#
#  $Version: 0.1 $  $Date: 2013/02/17 23:48:03 $
#
#   function to compute Indivdual Phylogenetic Species Richness Area Relationships (i.e., acumulation of PSR curves
#    around individual species) in a community.
#
#    Value: an object of class fv , bassicaly a data.frame with the following
#    elemnts:
#     r: vector of radii at which IPSRAR(r) has been estimated
#     ipsv: vector with the IPSRAR(r) values
#     ###npoints: number of focal points employed to estimate isar(r) at each 
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
#     tree:    A phylogenetic tree in "phylo" format (ape) or phylogenetic covariance matrix 
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

 




 
ipsrar<- function(mippp, mippp.sp=NULL, mimark=NULL,  namesmark=NULL, tree=NULL,
                      r=NULL, buffer=0, bfw=NULL, correct.phylo="mean") {

   
 
  #Check mippp
  #If species marks are within a data.frame of marks, get this column and 
  # discard the rest
  if(!is.null(namesmark)) mippp$marks <- factor((mippp$marks[namesmark][[1]]))  # OJO: que representa el [[1]] del final: es para coger los valores. si no se pone da problemas !!!!!!!!!!!!!!!!!!!!!!!!!!!
  #!!!!!!!!!!!!!!!!!!
  
    #Check tree  
    # this is for checktree
    idar <-"ipsrar"
    if(is.null(tree)) stop("you should provide a phylogenetic tree to compute ipscar")
    #if(class(tree)=="phylo") tree<- vcv.phylo(tree, corr = TRUE)  #!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  CAMBIADO 04/12/2019 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    if(inherits(tree, what="phylo")) tree<- vcv.phylo(tree, corr = TRUE)
    if(!is.matrix(tree)) stop("you should provide a phylogenetic tree or a phylogenetic covariance matrix to compute ipscar")
    tree <- checktree(tree=tree,  mippp=mippp, idar=idar, correct.phylo=correct.phylo)



  # Check mippp.sp
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
   
      # Check that there are enough differnt individuasl to compute psv (and avoid error!). Frequent problem for short r's
    cosamt.ok<- sapply(cosamt, function(y) any(apply(y,1,function(x) sum(x>0))>1)) #!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 17/02/2015
    
    
   
	
  #If indicated correct edge effect with an adaptive buffer
  if(buffer=="adapt"){
    bdp <- bdist.points(mippp.sp)
    for (i in 1:length(r)) cosamt[[i]]<-cosamt[[i]][bdp>=r[i],]
    npoints <- sapply(cosamt, function(x) dim(x)[1])
  }
  	
   
  #Compute IPSVAR(r) for each radius r
   # preset result to NA
    ipsrar<- rep(NA, length(r))
  
  
  
  # first for each individual  (for those individuals which have enough neighbors)) #!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 17/02/2015
  #ipsrar <- lapply(cosamt, function(x) psr(x, tree=tree, compute.var=F)$PSR)
  ipsrar[cosamt.ok] <- lapply(cosamt[cosamt.ok], function(x) psr(x, tree=tree, compute.var=F)$PSR)
  
  
  # then average ipsv for each radius
  ipsrar <- sapply(ipsrar, mean, na.rm=T) # TODO: optional na.rm. At short distances lot of NA's and maybe unreliable means

  # return results
  
   result <- data.frame(r=r, ipsrar=ipsrar)
  
   result <- fv(result, argu="r", ylab=substitute(IPSRAR(r), NULL),
              valu="ipsrar", fmla=ipsrar ~ r,
              alim=c(min(r),max(r)),
              labl=c("r", "%s(r)"),
              desc=c("radius of circle",
                     "%s"),
              fname="IPSRAR")

  
  
  
  #class(result)<-c("isar", class(result))
  return(result)
}
###############################################################################