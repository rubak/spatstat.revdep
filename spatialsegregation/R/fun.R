# last change : 090715
#################################################

#constant: the supported neighbourhoods
kGraphs<-c("geometric","knn","gabriel","delaunay","bgeometric")



#' Main computer function for spatialsegregation
#' 
#' Compute the spatial exposure (segregation vs. mingling) features from a given multitype point pattern. Usage of  shortcuts \code{minglingF}, \code{isarF}, \code{shannonF}, \code{simpsonF} etc. highly recommended.
#' 
#' @param X Multitype point pattern of class \code{ppp} (see package 'spatstat')
#' @param fun Default "isar". Takes "isar","mingling","shannon", "simpson", "mci" and "biomass", see below.
#' @param r Vector for the neighbourhood defining graph, e.g. "geometric" graph with different ranges. See below.
#' @param ntype Default "geometric". Type of the neighbourhood graph. Accepts: "knn", "geometric", "delauney", "gabriel".
#' @param funpars Default NULL. Parameter(s) for the measure. Mingling: c(i,j), where i= only for type i (0 for all), j=1 -> ratio version. ISAR: i, i=type (integer). Shannon: 0 or 1, see \code{v2} in \code{shannonF}. Simpson: none.
#' @param toroidal Default FALSE. If TRUE, use a toroidal correction in distance calculation. Works at the moment only for rectangular windows and "geometric" or "knn" graph.
#' @param minusRange If TRUE, adaptive minus-sampling is employed. Overrides \code{included}-vector. If given as a positive number, \code{included}-vector is created with points with distance atleast minusRange from the border.
#' @param included boolean-vector of length |pp|. included[i]==TRUE => pp[i] included in calculations. Used for minus-sampling border correction.
#' @param dbg Default FALSE. Print additional runtime texts.
#' @param doDists Default TRUE. Precalculate distances for speed. Be aware of memory requirement n*(n-1)/2
#' @param prepRange Default 0. If >0, shrink the search space for neighbourhoods by searching only points within distance R i.e. precalculates a geometric graph.
#' @param prepGraph Precalculated graph for the point pattern. If given, The \code{prepRange}, \code{dodists} and \code{toroidal} are ignored and calculations are carried using the prepGraph as a starting point. Useful for huge datasets.
#' @param prepGraphIsTarget If TRUE, precalculated graph \code{prepGraph} is used to calculate a single function value directly, all other neighbourhood parameters are ignored.
#' @param weightMatrix See \code{isarF} for this.
#' @param translate Use translation correction (see e.g. documentation of \code{spatstat::Kest} for details). Used only in mingling index.
#' @param ... ignored.
#' @details 
#' This is the general function for computing the spatial exposure (segregation/mingling) features. Used by \link{minglingF}, \link{shannonF}, \link{simpsonF}, \link{isarF}, \link{mciF} and \link{biomassF}, which should be preferred for better (and nicer) outcome. 
#' 
#' Possible neighbourhood relations include geometric, k-nearest neighbours, Delaunay, and Gabriel. Delaunay and Gabriel are parameter free, so given \code{r} has no meaning. In geometric graph, \code{r} is a vector of distances (sizes of the surrounding 'disc') and for k-nn \code{r} is the vector of neighbourhood abundances (so r is k in k-nn) for each point to consider in the calculation of the spatial exposure measures. In general, the basic type of spatial summary in literature uses 'geometric' graph  with several ranges. 
#'
#'For \code{geometric} and \code{knn}, the calculations are done by shrinking the graph given by the largest value of \code{r}. If dealing with large datasets, it is advisable to give preprocessing range, \code{prepRange}. The algorithm first calculates a geometric graph with parameter \code{prepRange}, and uses this as basis for finding the needed neighbourhoods. Speeds up calculations, but make sure \code{prepRange} is large enough e.g. in \code{geometric}, \code{prepRange}>max(\code{r}). \code{prepGraph}, if given, works as the preprocessed geometric graph (overrides prepRange), and can be computed using the package \code{spatgraphs}; useful for huge datasets, where the dominating graph needs to be computed only one. The \code{doDists} option speeds up calculations by precomputing the pairwise distances but requires  approx. n*(n-1)/2 * 32 bytes of memory. 
#'
#' For border correction, use \code{minusRange} for reduced border correction (for rectangular windows only). If using \code{geometric} or \code{knn} neighbourhoods, the option \code{toroidal} for toroidal correction is also available. The vector \code{included} can be given for more specific \code{minus}-correction, only those points with TRUE (or 1) value are used in calculation. However, the neighbourhoods will include all points.
#' 
#' @return 
#' 	Returns an object of class \code{fv}, see spatstat for more details. Basically a list with the computed values and parameter values.
#' 	
#' @export
#' @useDynLib spatialsegregation
#' @import spatstat

segregationFun<-function(X, fun="isar", r=NULL, ntype="geometric", funpars=NULL, 
		                 toroidal=FALSE, minusRange=TRUE,  included=NULL, dbg=FALSE, 
						 doDists=FALSE, prepRange=0.0, prepGraph=NULL, prepGraphIsTarget=FALSE, 
						 weightMatrix=NULL, translate=FALSE, ...)
# function types:
#	1 mingling
#   2 shannon
#   3 simpson
#   4 ISAR
#	5 MCI
#   6 biomass sum
{
	# a note about the neighbourhoods
	note<-NULL
	# turn the neighbourhood type into an integer
	if(ntype=="delauney")ntype<-"delaunay" # stupid mistake in the past, now corrected.
	ntypei<-charmatch(ntype, kGraphs) - 1 #minus1 for c-side...
	# check 
	if(is.na(ntypei))
		stop(paste("Error segregationFun: Wrong neighbourhood type, try one of:", 
						paste(kGraphs, collapse=",", sep=""),"."))
	if(ntypei==4) ntypei<-12
	# init neighbourhood parameters if not given: geometric range taken from Kest in spatstat
	# r not given
	if(is.null(r))
	{
		# snip
		if(ntypei==0)
		{
			W <- X$window
			npoints <-X$n
			lambda <- npoints/area.owin(W)
			rmaxdefault <- rmax.rule("G", W, lambda)
			breaks <- handle.r.b.args(r, NULL, W, rmaxdefault=rmaxdefault)
			rvals <- breaks$r
			rmax  <- breaks$max
			parvec<-rvals
		}
		# end snip
		# defaults for knn, gabriel and delauney. TODO: k-nn not yet clear what is a good range
		else
			parvec <- switch(ntype, "knn"=1:20, 0)
	}
	else parvec <- sort(r) # r given. make sure in ascending order.
	
	# TODO: for the c++-module coherence with 'spatgraphs', skip mass geometric
	if(ntypei>1) ntypei <- ntypei+1  
	
	# if a precalculated Graph is given, check its from spatgraphs
	if(class(prepGraph)!="sg" & !is.null(prepGraph) )stop("Error segregationFun: Prepared graph is not of class 'sg'.")
	else if(!is.null(prepGraph)) note<-c(note,paste("PrepGraph given, type ",prepGraph$type,", par ",prepGraph$par,sep="")) #add a note about prepGraph
	else note<-""
	 
	
	# turn the wanted function type into an integer
	funi <- charmatch(fun, (kFuns<-c("mingling","shannon","simpson","isar","mci","biomass")))
	if(is.na(funi)) stop("Error segregationFun: wrong function type.")
	
	# check the weightMatrix
	if(!is.null(weightMatrix)){
		if(sum(weightMatrix<0)>0)stop("Can't handle negative weights.")
		if(diff(dim(weightMatrix))!=0)	stop("Error segregationFun: weightMatrix not square matrix.")
		if(!isTRUE(all.equal(colnames(weightMatrix),rownames(weightMatrix))) |
				!isTRUE(all.equal(colnames(weightMatrix),levels(as.factor(X$marks))))){
			stop("Error segregationFun: weightMatrix row and column names and order should match X's mark levels.")
		}
	}
	else weightMatrix<-diag(0, ncol=length(levels(as.factor(X$marks))))-1 # no type-to-type weights given
	
	# modify the pp: convert factor labels to integers, add window area, add border distances. 	
	X<-sg.modify.pp(X)
	
	# if a minus sampling border correction is to be used, compute the ones to exclude
	# TRUE + geometric graph => adaptive
	# positive + any graph => fixed
	if(is.logical(minusRange)){
		if(minusRange){
			if(ntypei!=0) stop("Error segregationFun: adaptive minus-sampling only for geometric neighbourhood.")
			included<-rep(-1, X$n)
			note<-c(note, "adaptive minus-sampling border correction")
		}
	}
	else if(minusRange>0)
	{
		included<- X$bdist > minusRange
		note<-c(note,paste("Minus correction, fixed radius=", minusRange, ";", sep=""))
	}
	# if we accept that all points are good for computation, included vector is all 1's
	if(is.null(included) | length(included)!=X$n) included<-rep(1,X$n)
	
	# if translation corrected: Only for mingling and Simpson
	if(translate){
		note<-c(note,paste("translation correction", sep=""))
		if(funi!=1 & funi!=3) stop("Translation correction available only for Mingling and Simpson functions and indices.")
	}
	
	# check if the prepGraph is given when used as the target neighbourhood configuration
	if(prepGraphIsTarget && is.null(prepGraph)) stop("Error segregationFun: prepGraph not given but needed for calculation.")
		
	# this is for the c-side check to know if we are giving a prepGraph
	prepGraph$'isnull'<- as.integer(is.null(prepGraph))
	
	# set parameters according to prepGraph
	if(prepGraphIsTarget)
	{
		parvec<-prepGraph$parameters
		ntype<-prepGraph$type
	}
	
	# the main call 
	res<-.External("fun_c", as.integer(dbg), X, as.numeric(funpars), 
					as.integer(ntypei), as.numeric(parvec), 
					as.integer(funi), as.integer(toroidal), 
					as.numeric(prepRange), as.integer(doDists),
					as.integer(translate),
					as.integer(included), prepGraph, 
					as.integer(prepGraphIsTarget),
					as.numeric(weightMatrix),
					PACKAGE="spatialsegregation")

	
	# turn the result list into a matrix: 1 col per valuetype, 1 row per parameter
	a<-t(matrix(unlist(res),ncol=length(parvec)))

	# name the rows par1, par2, ...
	#rownames(a)<-paste("par",1:length(parvec),sep="")#no need for this anymore
 
	# the unit name: different from Kest etc, reports the meaning of the unit depending on -
	#   the neighbourhood type
	unitname<-switch(ntypei,"1"="range","2"="neighbour",NULL)
	list(v=a, included=as.logical(included), 
		 parvec=parvec, unitname=unitname, 
		 ntype=ntype, note=note, point.values=X$mass2 # point.values2=X$mass2
        )
}


#####################
#' Prepare the input data pattern object
#' 
#' Prepares given point pattern object for computations.
#' 
#' @param pp Point pattern object
#' 
#' @export
sg.modify.pp<-function(pp)
{
	n<-length(pp[["x"]])
	pp[["n"]]<-n
	
	if(is.null(pp[["window"]])) pp[["window"]]<-list(xrange=range(pp[["x"]]), yrange=range(pp[["y"]]), type="rectangle")
	if(!is.null(dim(pp$marks)))stop("Error modify.pp: data.frame-marks not yet supported. Put factor marks to $marks or $types, and continuous marks to $mass.")
	if(length(pp[["mass"]]) != n ) # set the masses
	{
		if(length(pp[["marks"]])< n | !is.numeric(pp[["marks"]])) pp$mass<-rep(1.0,n)
		else pp$mass<-pp$marks
	}
	
	
	if(length(pp[["types"]]) < n) # set the types
	{
		if( (is.factor(pp$marks) | (is.integer(pp$marks)) & length(pp[["marks"]])==n) ) x<-as.factor(pp$marks) 
		else x<-as.factor(rep(1,n))
	}
	else x<-as.factor(pp$types)
	pp$mass<-as.numeric(pp$mass)
	pp$types<-as.integer(x)
	pp$marks<-NULL
	
	if(is.null(pp[["z"]]) || length(pp[["z"]])!=length(pp[["x"]])) pp$z<-rep(0.0,n) # if 2D only
	if(is.null(pp[["window"]][["z"]])) pp$window$z<-as.numeric(c(0.0,1.0)) # if 2D only
	
	pp$window$x<-as.numeric(pp$window$x)
	pp$window$y<-as.numeric(pp$window$y)
	pp$window$z<-as.numeric(pp$window$z)
	pp$x<-as.numeric(pp$x)
	pp$y<-as.numeric(pp$y)
	pp$z<-as.numeric(pp$z)	
	
	## add the distance to border
	pp$bdist<-bdist.points(pp)
	# include area
	pp$area<-as.numeric(area.owin(pp$window))
	
	# add a container for pointwise values
	pp$mass2 <- rep(0, pp$n)
	#
	pp
}

#####################

