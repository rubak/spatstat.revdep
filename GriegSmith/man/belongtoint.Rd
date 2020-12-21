\name{belongtoint}
\alias{belongtoint}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{belongtoint
%%  ~~function to do ... ~~
}
\description{Used by GriegSmith function
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
belongtoint(bin.vals.vect, int.x, int.y, vect)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{bin.vals.vect}{
%%     ~~Describe \code{bin.vals.vect} here~~
}
  \item{int.x}{
%%     ~~Describe \code{int.x} here~~
}
  \item{int.y}{
%%     ~~Describe \code{int.y} here~~
}
  \item{vect}{
%%     ~~Describe \code{vect} here~~
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
function(bin.vals.vect,int.x,int.y,vect){

	xval<-bin.vals.vect[1];
	yval<-bin.vals.vect[2];	



	sum(vect[,1] >= xval & vect[,1] < (xval+int.x) & vect[,2] >= yval & vect[,2] < (yval+int.y));
 
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
