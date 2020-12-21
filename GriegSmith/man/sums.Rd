\name{sums}
\alias{sums}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{sums
%%  ~~function to do ... ~~
}
\description{Used by GriegSmith function
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
sums(coords, dim, xmin = min(coords[, 1]), xmax = max(coords[, 1]), ymin = min(coords[, 2]), ymax = max(coords[, 2]))
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{coords}{
%%     ~~Describe \code{coords} here~~
}
  \item{dim}{
%%     ~~Describe \code{dim} here~~
}
  \item{xmin}{
%%     ~~Describe \code{xmin} here~~
}
  \item{xmax}{
%%     ~~Describe \code{xmax} here~~
}
  \item{ymin}{
%%     ~~Describe \code{ymin} here~~
}
  \item{ymax}{
%%     ~~Describe \code{ymax} here~~
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
function (coords,dim,xmin=min(coords[,1]),xmax=max(coords[,1]),ymin=min(coords[,2]),ymax=max(coords[,2])){


	xints<-((xmax-xmin)/dim);
	yints<-((ymax-ymin)/dim);

	xbins<-seq(from=xmin, to=xmax-xints, by=xints);
	ybins<-seq(from=ymin, to=ymax-yints, by=yints);



	bins<-cbind(c(sapply(xbins,rep,dim)), rep(ybins,dim));
	cnts<-matrix(apply(bins,1,belongtoint,vect=coords,int.x=xints,int.y=yints),nrow=dim,byrow=TRUE);


  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
