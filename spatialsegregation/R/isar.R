# Last update: 090715
###############################################################################

#' Individual Species Area Relationship
#' 
#' Compute the Individual Species Area Relationship ( ISAR ) or Local Species
#' Richness, for a given multitype point pattern.
#' @param X Multitype point pattern of class \code{ppp} (see package 'spatstat')
#' @param r Vector of sizes for neighbourhoods, e.g. \code{geometric} graph with
#'   different ranges.
#' @param target Default NULL. Calculate only for target type. If NULL computes
#'   for each type + mean over all types.
#' @param v2 Logical. Estimate species-to-neighbours-ratio instead of just total
#'   number of species.
#' @param v3 Logical. Instead of summing number 1 for each species present, sum
#'   the average X$mass of each species present.
#' @param v4 Logical. Estimate ISAR using empty space probabilities instead of
#'   direct counts (equals the normal version in all my tests)
#' @param ntype Sets the n'hood type to \code{knn} by default in isar.index.
#' @param ... Further parameters for the function \code{segregationFun}.
#'   
#' @details Extension of ISAR-function introduced in WGGH07. In effect
#' calculates the expected amount of different types present in the
#' neighbourhood of a point in the pattern.
#' 
#' The function \code{isarF} is the calculation function for different
#' neighbourhoods. Uses function \code{\link{segregationFun}}.
#' 
#' The function \code{isar.index} is a shortcut to get a single value for the
#' pattern. Uses 4-nn graph by default.
#' 
#' @references 
#' Rajala, Illian: A family of spatial biodiversity measures based on graphs, Env. Ecol. Stat. 2012
#' 
#' Wiegand, Gunatilleke, Gunatilleke, Huth: How individual species structure diversity in tropical forests. PNAS, nov 16, 2007. 
#' 
#' @export
isarF<-function(X, r=NULL, target=NULL, v2=FALSE, v3=FALSE, v4=FALSE, ... )
{
	# check that X is ppp-object
	verifyclass(X, "ppp")
	if(length(levels(X$marks))<2) stop("Use only on a multitype point pattern (data.frame-marks not yet supported).")
	# if no target given, calculate for all types
	if(is.null(target))
	{
		targeti <- 0
		valu   <- "ISARmean"
	}
	# else convert to an integer
	else
	{
		if(!is.factor(X$marks))warning("Marks of X are not in factor form. Transforming.")
		X$marks<-as.factor(X$marks)
		targeti<- which( levels(X$marks)  == target)
#		targeti<-which( union(X$marks, NULL) == target)
		if(length(targeti)!=1) stop("Target type not one of pattern types.")
	}
	
	#v2 logical if a degree weighted version should be calculated
	if(v2) funtype <- "Neighbour-count-weighted-ISAR"
	if(v3){v2<-2; funtype <- "Mass-weighted ISAR"; }
	if(v4){v2<-3; funtype <- "eISAR"; }
	else funtype <- "ISAR"
	
		
	# use the main calc function
	res<-segregationFun(X=X, fun="isar", r, funpars=c(targeti, as.integer(v2)), ...)
	
	# theoretical values in CSR: depends on the neighbourhood type
	ntype<-res$ntype
	mdeg<-function(l,k)c( pi*l*k^2, k, 4, 6)[charmatch(ntype, kGraphs)]
	
	# get the intensities
	sum0<-summary(X)
	omarks<-order(union(X$marks,NULL)) # right order of marks
	l<-sum0$marks[,3][omarks]
	#    calc the theoretical values, also for the degree weighted version
	theo<-NULL
	for(para in res$parvec)theo<-
				c(theo, 
                 sum(1-exp(-mdeg(sum(l),para)*l/sum(l) )) / ifelse(v2,mdeg(sum(l), para),1))
										
		
	# create the fv-object
	isar.final<-fv(data.frame(theo=theo,par=res$parvec), 
			       argu="par",
				   alim=range(res$parvec),
				   ylab=substitute(ISAR, NULL),
			       desc=c("CSR values","Parameter values"),
				   valu="theo",
				   fmla=".~par",
			   unitname=res$unitname,
	   			  fname=funtype
			           )
					   
	# add all typewise values if no target type given
	if(targeti==0)
	{
		# the values from calculation
		tw<-res$v
		# set the names right, and don't forget to check inclusion (might drop some types off)
		colnames(tw)<-union(marks(X[res$included]),NULL)
		isar.final<-bind.fv(x=isar.final,
						    y=tw,
					  	 desc=paste("Typewise",funtype,"for type",colnames(tw)),
					     labl=colnames(tw)
					         )
					
		isar.final<-bind.fv(x=isar.final,
				            y=data.frame("ISARmean"=apply(res$v,1,mean,na.rm=TRUE)),
				         desc=paste("Mean",funtype,"over types"),
						 labl="ISARmean",
						 preferred="ISARmean"
				 	     	 )		
		# a frequency weighted mean instead of just a mean, w=freqs/sum(freqs)
		#Iw=apply(res$v,2,weighted.mean,w=w,na.rm=TRUE), 
	}
	
	# if target type given add the values for the target type
	else
	{
		isar.final<-bind.fv(x=isar.final,
				            y=data.frame("ISAR"=res$v[,1]),
						 desc=paste(funtype,"for type",target),
						 labl="ISAR",
					preferred="ISAR"
								)
	}

	# attach the frequencies too
	attr(isar.final,"frequencies")<-freqs(X[res$included])
	
	# and some notes
	attr(isar.final,"neighbourhoodType")<-res$ntype
	attr(isar.final,"note")<-res$note
	# add pointwise values
	attr(isar.final,"point.values")<-res$point.values
	# return 
	isar.final
}

###############################################################################
#' @export
#' @describeIn isarF Shortcut for 4-nearest neighbour value.
isar.index<-function(X, r=4, ntype="knn", ...)
{
	if(length(r)>1)stop("Use isarF for vector of parameter values.")
	I0<-isarF(X, r=r, ntype=ntype, ...)
	data.frame(meanISAR=I0$I, par=I0$par)
}