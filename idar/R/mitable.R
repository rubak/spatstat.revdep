
##############################################################################
#
#   mitable.R
#   function to tabulate number of (occurence of) different species (i.e.
#    types) of a multivariate point pattern ppp2 within circles of radius r
#    around each individuals of a (unmarked) point pattern ppp1
#
#
#    Value: a list of length r; each element of the list a data.frame with
#            nrows = number of points in ppp1 and ncols = unique(ppp2$marks),
#            each cell being the counts of each type of pp2 points within a
#            radius r of the corresponding point of ppp1
#    Arguments:
#     ppp1: Focal (unmarked) point pattern (ppp of spatstat) or data.frame
#           with columns named "x" an "y"
#     ppp2: Multivariate point pattern (ppp of spatstat) or data.frame whith 
#           columns named "x", "y" and "marks"
#     r:    vector with the sequence of radii (>0) at which tabulate counts
#
#



mitable <- function(ppp1,ppp2,r){
    
    if(is.null(ppp2$marks)) stop("ppp2 should be a multivariate ppp or a data.frame with a '$marks' column")
    
    x1 <- ppp1$x
    y1 <- ppp1$y
    
    x2 <- ppp2$x
    y2 <- ppp2$y
    
    l1<- length(x1)
    l2<- length(x2)

    sp<- factor(ppp2$marks)
    especies <-levels (sp)

    nombresp <-  as.numeric(sp)
    nomsp <- unique(nombresp)
    nsp <-length(nomsp)
    nr <-length(r)
    ltab<-l1*nsp #length of a the matrix ppp1*nsp
 
    
    # mitable.f renamed mitablee.f to avoid troubles with registrations etc 27.04.2017
    ans <- .Fortran('mtb', x1=as.double(x1),y1=as.double(y1),
                    x2=as.double(x2),y2=as.double(y2),
                    l1=as.integer(l1), l2=as.integer(l2),
                    nombresp=as.integer(nombresp), nsp=as.integer(nsp),
                    r=as.double(r), nr=as.integer(nr),ltab=as.integer(ltab),
                    abu=double(ltab*nr),
                    PACKAGE="idar")
    
       
      
      result=vector(mode="list", length=nr)
      names(result)=r
      for(i in 1:nr){
        result[[i]]<-matrix(ans$abu[(((i-1)*ltab)+1):(i*ltab)],l1,nsp,
	                        byrow=TRUE,dimnames=list(NULL, 
				     especies=especies))
      }    
      return(result)
}