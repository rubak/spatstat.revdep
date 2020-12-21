\name{addem}
\alias{addem}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{addem
%%  ~~function to do ... ~~
}
\description{Used by GriegSmith function
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
addem(startingvals, xmatlen, ymatlen, data)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{startingvals}{
%%     ~~Describe \code{startingvals} here~~
}
  \item{xmatlen}{
%%     ~~Describe \code{xmatlen} here~~
}
  \item{ymatlen}{
%%     ~~Describe \code{ymatlen} here~~
}
  \item{data}{
%%     ~~Describe \code{data} here~~
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
function(startingvals,xmatlen,ymatlen,data){

	x<-startingvals[1];
	y<-startingvals[2];


	#print("here");
	#print(data[x:(x+xmatlen-1),y:(y+ymatlen-1)]);

	sum(data[x:(x+xmatlen-1),y:(y+ymatlen-1)])^2;




  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
