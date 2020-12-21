\name{sumofsquares}
\alias{sumofsquares}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{sumofsquares
%%  ~~function to do ... ~~
}
\description{Used by GriegSmith function
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
sumofsquares(sizematrix, singlecounts)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{sizematrix}{
%%     ~~Describe \code{sizematrix} here~~
}
  \item{singlecounts}{
%%     ~~Describe \code{singlecounts} here~~
}
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
}
\author{
%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{

## The function is currently defined as
function(sizematrix,singlecounts){

#	print(sizematrix);

	xsublength<-sizematrix[1];
	ysublength<-sizematrix[2];


	xsize<-length(singlecounts[1,]);
	ysize<-length(singlecounts[,1]);


	
	xmin<-rep(seq(from=1,to=xsize,by=xsublength),ysize/ysublength);	
	ymin<-sort(rep(seq(from=1,to=ysize,by=ysublength),xsize/xsublength));

	
	submatricies<-cbind(xmin,ymin);
	
	squaredsums<-sum(apply(submatricies,1,addem,data=singlecounts,xmatlen=xsublength,ymatlen=ysublength));


#	print(squaredsums);

	## sum up all the numbers in each matrix, square those numbers and add them;
	


  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
