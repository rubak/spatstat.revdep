\name{sample.points}
\alias{sample.points}
\alias{hndetfct}
\title{Sample points within each transect and filter with specified detection function}
\description{
Create a dataframe of observations by simulating distance sampling of a point process
with a systematic set of lines over a rectangular grid.  The \code{transects},
\code{lines} and point process(\code{points.ppp}) are input arguments.  Detection
of observations is specified with a user-defined detection function which takes
a distance vector and set of parameters \code{det.par} as its arguments.
}
\usage{
sample.points(transects,lines,points.ppp,detfct=NULL,det.par=NULL,
                      det.formula=~1,covariates=NULL)
hndetfct(x,scale)
}
\arguments{
  \item{transects}{list of transect polygons}
  \item{lines}{dataframe of lines}
  \item{points.ppp}{simulated point process}
  \item{detfct}{detection function name}
  \item{det.par}{parameters for the detection function}
  \item{det.formula}{formula of covariates to use for scale of distance
               if det.formula=~-1, uses a strip transect}
  \item{covariates}{a matrix with columns x,y and any number of covariates
                  x and y are the mid points of the grid cells; the order
                  of the rows must match the formulation for function im}
  \item{x}{perpendicular distance for detection function}
  \item{scale}{scale for detection function}
  
}

\value{observation dataframe with fields label,x,y,distance for line label, x,y coordinates 
of the observation location and its perpendicular distance from the line}
\details{
Definition for half-normal detection function (hndetfct) is \code{exp(-(x^2/(2*exp(scale)^2)))}
}
\author{Jeff Laake}
\seealso{ \code{\link{simCovariates}},\code{\link{simPts}},\code{\link{create.lines}}}
\examples{
study.area=owin(xrange=c(0,100),yrange=c(0,100))
hab.range=30
probs=c(1/3,2/3)
covariates = simCovariates(hab.range, probs)
xlines=create.lines(study.area,nlines=10,width=5,angle=45)
ls=lines_to_strips(xlines,study.area)
plot(ls$lines,lty=2)
plot(owin(poly=ls$transects),add=TRUE)
xpp=simPts(covariates=covariates,int.formula=~factor(habitat),int.par=c(0,1,2),EN=1000)
obs=sample.points(transects=ls$transects,lines=xlines,points.ppp=xpp,
                    hndetfct,c(1),covariates=covariates)
plot(ppp(x=obs$x,y=obs$y,window=study.area),add=TRUE,pch=20)
}

