% Copyright 2001-4 by Roger S. Bivand
\name{leglabs}
\alias{leglabs}
\title{Make legend labels}
\description{
\code{leglabs} makes character strings from the same break points. The \code{plot.polylist()} function may be used as a generic S3 method. }
\usage{
leglabs(vec, under="under", over="over", between="-", reverse=FALSE)
}

\arguments{
  \item{vec}{vector of break values}
  \item{under}{character value for under}
  \item{over}{character value for over}
  \item{between}{character value for between}
  \item{reverse}{flag to reverse order of values, you will also need to reorder colours, see example}
%  \item{y}{numeric vector to be broken into class intervals}
%  \item{rightmost.closed}{as in findInterval()}
%  \item{all.inside}{as in findInterval()}
}
\author{Roger Bivand \email{Roger.Bivand@nhh.no}}

\seealso{\code{\link{findInterval}}}

\examples{
mappolys <- readShapeSpatial(system.file("shapes/columbus.shp", package="maptools")[1], ID="NEIGNO")
brks <- round(quantile(mappolys$CRIME, probs=seq(0,1,0.2)), digits=2)
colours <- c("salmon1", "salmon2", "red3", "brown", "black")
plot(mappolys, col=colours[findInterval(mappolys$CRIME, brks,
 all.inside=TRUE)])
legend(x=c(5.8, 7.1), y=c(13, 14.5), legend=leglabs(brks),
  fill=colours, bty="n")
title(main=paste("Columbus OH: residential burglaries and vehicle",
 "thefts per thousand households, 1980", sep="\n"))
#legend with reversed order
plot(mappolys, col=colours[findInterval(mappolys$CRIME, brks,
 all.inside=TRUE)])
legend(x=c(5.8, 7.1), y=c(13, 14.5), legend=leglabs(brks, reverse = TRUE),
  fill=rev(colours), bty="n")
title(main=paste("Columbus OH: residential burglaries and vehicle",
 "thefts per thousand households, 1980 (reversed legend)", sep="\n"))
}
\keyword{spatial}
