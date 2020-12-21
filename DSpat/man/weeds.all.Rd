\name{weeds.all}
\alias{weeds.all}
\docType{data}
\title{Dubbo weed data with constructed y-coordinate}
\description{
Locations of devils claw in a farming paddock.  Locations to all weeds are given
as x,y coordinates and are contained in one of eight 150m wide transects (75m each side). The
weeds seen by observers are specified as Seen=1.}
\usage{data(weeds.all)}
\format{
  A data frame with 742 observations on the following 5 variables.
  \describe{
    \item{\code{label}}{label of the transect 1 to 8}
    \item{\code{x}}{x coordinate along horizontal (east-west) of 1200m x 1200m paddock}
    \item{\code{y}}{y coordinate along vertical (north-south) of 1200m x 1200m paddock}
    \item{\code{distance}}{absolute perpendicular distance from line}
    \item{\code{Seen}}{weed was seen if 1 and 0 if missed}
  }
}
\details{
The data provided from Melville and Welsh did not have the y-coordinate.  We have constructed
y-coordinates by drawing randomly from a uniform distribution in the y-direction such that no two weeds are at
the exact same location.  The code used to create \code{weeds.all} was as follows:

\preformatted{
data(weeds)
# Fudge the data ever so slightly to appease spatstat so the transects don't abut
# and all points are contained within the strips.
weeds$SignedDistance[weeds$SignedDistance==75]=74.99
weeds$SignedDistance[weeds$SignedDistance==-75]=-74.99
weeds.all=data.frame(label=weeds$Transect, x=(weeds$Transect-1)*150+75
               +weeds$SignedDistance, y=floor(runif(dim(weeds)[1])*1200),
               distance=weeds$Distance, Seen=weeds$Seen)
while(any(duplicated(data.frame(x=weeds.all$x,y=weeds.all$y))))
{
  npts=sum(as.numeric(any(duplicated(data.frame(x=weeds.all$x,y=weeds.all$y)))))
  weeds.all$y[duplicated(data.frame(x=weeds.all$x,y=weeds.all$y))]=
                                                               runif(npts)*1200
}
save(weeds.all,file="weeds.all.rda")
}

See \code{\link{weeds}} for more details.
}
\references{
Melville, G. J., and A. H. Welsh. 2001. Line transect sampling in small regions. Biometrics 57:1130-1137.
}
\keyword{datasets}
