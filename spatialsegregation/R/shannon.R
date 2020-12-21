#' Spatial Shannon index
#' 
#' Compute the spatial and aspatial Shannon index for a given multitype point pattern.
#' 
#' @param X Multitype point pattern of class \code{ppp} (see package 'spatstat')
#' @param r Vector of sizes for neighbourhoods, e.g. \code{geometric} graph with different ranges.
#' @param spatial If FALSE, return the classical aspatial index value.
#' @param v2 If TRUE, use the real number of types in neighbourhoods as the log-base instead of total population type count.
#' @param ... Further parametes for the function \code{\link{segregationFun}}.
#' 
#' 
#' @details 
#' 
#' The form of Shannon index is  \var{H = 1 - E(o)/E(N)}, where \var{E(N)} is the global entropy and \var{E(o)} is the local entropy calculated as \var{E(o)= - sum pi_tau log(pi_tau)}, where the sum is over the different types present in the pattern, and \var{pi_tau} is the expected frequency of type \var{tau} points in a neighbourhood of a typical point of the pattern. 
#' 
#' The function \code{shannonF} is the calculation function. Uses function \code{\link{segregationFun}}. 
#' 
#' The function \code{shannon.index} is a shortcut to get the non-spatial Shannon index. 
#' 
#' @return 
#'  Returns an \code{fv}-object, see \code{spatstat} for more information. The index returns a scalar.
#'  
#' @references 
#' Rajala, Illian: A family of spatial biodiversity measures based on graphs, Env. Ecol. Stat. 2012
#'
#' Reardon, O'sullivan: Measures of spatial segregation. Sociological methodology, 34:121-162, 2004.
#' @export

shannonF<-function(X, r=NULL, v2=FALSE, ...)
#Shannon index for graphs, with possibly various range-parameters
{
	# check that X is multitype ppp-object
	verifyclass(X, "ppp")
	if(length(levels(X$marks))<2) stop("Use only on a multitype point pattern (data.frame-marks not yet supported).")
	
	# the main calc function: it calculates only the pi_tau-vector
	res<-segregationFun(X=X, r=r, fun="shannon", funpars=ifelse(v2,1,0), ...) 
	
	# calc the aspatial (global) entropy i.e. shannon index
	eglobal <- -shannon.index(X,spatial=FALSE)
	
	# a function to calculate  -1*sum( pi_tau*log(pi_tau) )
	f<-function(pvec) 
	{
		E1<-0
		pvec<-pvec/sum(pvec)
		S<-length(pvec)
		ok<-pvec>0
		E1 <- sum(pvec[ok]*log(pvec[ok],base=S))
		(1-E1/eglobal)
	}
	
	# if we take the log base as the individual degree instead of total S -
	#   then the result is the mean
	if(v2)
	{
		H<-apply(-res$v,1,mean)
		desc<-"Spatial Shannon index, v2 (check documentation)"
	}
	# the base is S: calc the E(o)/E
	else
	{
		H<-unname(apply(res$v,1,f))
		desc<-"Spatial Shannon index"
	}
		
	# the CSR values
	theo<-rep(0, length(res$parvec))
	
	# create the fv-object
	shannon.final<-fv(data.frame(theo=theo, par=res$parvec, H=H),
			          argu="par",
					  alim=range(res$parvec),
					  ylab=substitute(H,NULL),
					  desc=c("CSR values","Parameter values",desc),
					  valu="H",
					  fmla=".~par",
			      unitname=res$unitname,
					 fname=desc
					 )
	
	# add the typewise pii_tau values of which the index is a summary
	attr(shannon.final,"typewise")<- res$v
	
	# add the global index value too
	attr(shannon.final,"Aspatial Shannon index")<--eglobal
	
	# add the note about neighbourhood definition
	attr(shannon.final,"note")<-res$note
	
	# return
	shannon.final
}


####################################################################################################
#' @export 
#' @describeIn shannonF Traditional index.
shannon.index<-function(X, spatial=FALSE, ...)
{
	#the traditional aspatial Shannon index
	if(!spatial)
	{  
		sum0<-summary(X)
		ints<-sum0$marks[,3]
		m<-union(X$marks,NULL)
		pii<-ints/sum(ints)
		H<- (-1)*sum( pii * log(pii, length(m)),na.rm=TRUE)
		names(H)<-"Non-spatial Shannon index"
	}
	if(spatial)
	{			   #spatial Simpson index for a set of edgelists	
		H<-shannonF(X=X, ...)$H
		names(H)<-"Spatial Shannon index"
	}
	H
}
