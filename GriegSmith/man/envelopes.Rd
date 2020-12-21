\name{envelopes}
\alias{envelopes}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{envelopes
%%  ~~function to do ... ~~
}
\description{Used by GriegSmith function
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
envelopes(iterations = 100, countdata, dimention)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{iterations}{
%%     ~~Describe \code{iterations} here~~
}
  \item{countdata}{
%%     ~~Describe \code{countdata} here~~
}
  \item{dimention}{
%%     ~~Describe \code{dimention} here~~
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
function(iterations=100,countdata,dimention){



	## Randobly arrange the counts;



	GSpermprev<-array();

	for (i in 1:iterations){
	
		perm<-matrix(sample(as.vector(countdata),size=2^(dimention*2),replace=FALSE),nrow=2^dimention,byrow=TRUE);
		GSperm<-iterate(perm,dimention);

		
		GSpermprev<-cbind(GSpermprev,GSperm[,3])
	


	}


	ret.val<-cbind(apply(GSpermprev[,-1],1,quantile,probs=c(.05)),apply(GSpermprev[,-1],1,quantile,probs=c(.95)))


  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
