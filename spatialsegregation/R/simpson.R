# last change: 090715
###############################################################################

#' Spatial Simpson index
#' 
#' Compute the spatial and non-spatial Simpson index for a given multitype point pattern.
#' 
#' @param X Multitype point pattern of class \code{ppp} (see package 'spatstat')
#' @param r Vector of sizes for neighbourhoods, e.g. \code{geometric} graph with different ranges.
#' @param spatial If FALSE, return the classical aspatial index value.
#' @param ... Further parametes for the function \code{\link{segregationFun}}.
#' 
#' @details 
#' 
#' The form of Simpson index is  \var{S = 1 - sum pi_tau}, where the sum is over the types of the pattern, and \var{pi_tau} is like
#' in Shimatani\& Kubota 2004.
#' The function \code{simpsonF} is the main calculation function. Uses function \code{\link{segregationFun}}.
#' 
#' The function \code{simpson.index} is a shortcut to get a single value for the pattern using 4-nearest neighbours graph by default. 
#'
#' @return
#' If spatial, returns an \code{fv}-object, see \code{spatstat} for more information. Otherwise a numeric value.
#'
#' @references 
#' Rajala, Illian: A family of spatial biodiversity measures based on graphs, Env. Ecol. Stat. 2012
#' 
#' Shimatani, Kubota: Quantitative assesment of multispecies spatial pattern with high species diversity. Ecological Research, 19, 2004.
#'   
#' @aliases simpson.index
#' @export

simpsonF<-function(X, r=NULL, ...)
{
	# check that X is multitype ppp-object
	verifyclass(X, "ppp")
	if(length(levels(X$marks))<2) stop("Use only on a multitype point pattern (data.frame-marks not yet supported).")
	
	# the main calc function: returns the typewise mean of (deg_i(o)/deg)^2
	res<-segregationFun(X, r=r, fun="simpson", ...)
	
	# calc the non-spatial (global) value
	aspat<-simpson.index(X,spatial=FALSE)
	
	#TODO: the CSR values: possibly not right, the mean degree is not included
	theo<-rep(aspat, length(res$parvec))
	
	# calc the spatial value
	S <- 1 - unname(rowSums(res$v))
	
	# create the fv-object
	simpson.final<-fv(data.frame(theo=theo, par=res$parvec, S=S),
			argu="par",
			alim=range(res$parvec),
			ylab=substitute(S,NULL),
			desc=c("CSR values","Parameter values","Spatial Simpson index"),
			valu="S",
			fmla=".~par",
			unitname=res$unitname,
			fname="Spatial Simpson index"
	)
	
	
	
	# include also the typewise values of which the index is a summary
	attr(simpson.final,"typewise")<- -res$v
	
	# add the global index value too
	attr(simpson.final,"Aspatial Simpson index")<-aspat
	
	# add the note about neighbourhood definition
	attr(simpson.final,"note")<-res$note
	
	# return
	simpson.final
}

####################################################################################################
#' @export
#' @describeIn simpsonF The Spatial Simpson Index 
simpson.index<-function(X, spatial=FALSE, ...)
{
	#the traditional aspatial Simpson index 1-D
	if(!spatial)
	{			
		sum0<-summary(X)
		ints<-sum0$marks[,3]
		m<-union(X$marks,NULL)
		pii<-ints/sum(ints)
		N<-X$n
		n<-pii*N		
		S<- 1 - sum( n*(n-1)  )/(N*(N-1))#Simpson index of diversity
		names(S)<-"Non-spatial Simpson index"
	}
	if(spatial)
	{			   #spatial Simpson index for a set of edgelists	
		S<-simpsonF(X, ...)$S
		names(S)<-"Spatial Simpson index"
	}
	S
}
