\name{GriegSmith}
\alias{GriegSmith}
\alias{plot.GriegSmith}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Grieg-Smith Calculation
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
This function accepts a point process object, a two column matrix of x-y coordinate pairs or a three column matrix containing x-y coordinates and quadrat counts in the third column. If the data contains quadrat counts, then the counts=TRUE option must be selected.
The function returns a GriegSmith object which is a matrix with block sizes, sum of squares for each block size as well as mean sums of squares. Simulation envelopes are produced as well through randombly permuting the quadrat counts. The 5th and 95th percentiles of the permutations create the simulation envelope.

Ploting the GriegSmith object produces a plot of the MSr as well as the simulation envelopes. 

}
\usage{
GriegSmith(datapoints, counts = FALSE, env = 100)
\method{plot}{GriegSmith}(x,main, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{datapoints}{
%%     ~~Describe \code{datapoints} here~~
datapoints is either a point process object (ppp), a two column matrix of x-y coordinates, or a three column matrix of x-y quadrat coords with a third column of quadrat counts.

}
  \item{counts}{
%%     ~~Describe \code{counts} here~~

If datapoints is a three column matrix with quadrat counts, then set counts=T

}
  \item{env}{
%%     ~~Describe \code{env} here~~
How many permuatations should be used to create the simulation envelopes default=100.
}
\item{x}{A GriegSmith object created with GriegSmith()}
\item{main}{the graph title}
\item{...}{other parameters passed to the plot function}

}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
%% ~put references to the literature/web site here ~
Statistical Methods for Spatial Data Analysis. Oliver Schabenberger and Carol A. Gotway . Boca Raton, FL: Chapman & Hall/CRC, 2005.

Greig-Smith, P. 1952. The use of random and contiguous quadrats in the study of structure in plant communities. Annals of Botany 16:293-316.
}
\author{
%%  ~~who you are~~
Brian McGuire
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{

data(amacrine,package="spatstat")
GS_ama<-GriegSmith(amacrine);
plot(GS_ama)




## The function is currently defined as
function(datapoints,counts=FALSE,env=100){


	if(counts==FALSE){

	if(is.ppp(datapoints)){

		xmin<-datapoints$window$xrange[1]
		xmax<-datapoints$window$xrange[2]

		ymin<-datapoints$window$yrange[1]
		ymax<-datapoints$window$yrange[2]


		datapoints<-cbind(datapoints$x,datapoints$y);


		
	}
	else{

		xmax<-max(datapoints[,1]);
		xmin<-min(datapoints[,1]);

		ymax<-max(datapoints[,2]);
		ymin<-min(datapoints[,2]);





	}


		
	
	numpts<-length(datapoints[,1]);
	startingdim<-ceiling(log(numpts)/(2*log(2)));
	counts<-sums(datapoints,2^startingdim,xmin,xmax,ymin,ymax)
		

	}
	else {
		if (max(datapoints[,1]) != max(datapoints[,2])) stop("Your count data must have equal dimensions")


		datapoints<-datapoints[order(datapoints[,2],datapoints[,1]),]
		numpts<-sum(datapoints[,3])
		startingdim<-ceiling(log(max(datapoints[,1]))/log(2))
		counts<-matrix(datapoints[,3],nrow=2^startingdim,byrow=TRUE);
	
		
	}	


	actual<-iterate(counts,startingdim);
	sims<-envelopes(env,counts,startingdim);
	final<-cbind(actual,sims);
	

	colnames(final)<-c("blocksize","SSr","MSr","MSr.05","MSr.95");

	
	class(final) <- "GriegSmith"


	final;

	


  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
