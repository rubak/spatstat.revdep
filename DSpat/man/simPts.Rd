\name{simPts}
\alias{simPts}
\title{Simulates point process on a rectangular grid}
\description{
Generates a set of points from either a homogeneous or inhomogeneous Poisson process
with optional clustering.  This is a wrapper function for \code{rpoispp} from
the \code{spatstat} package. The intensity is defined by covariates on a grid,
an intensity formula and parameters. The correlation structure is defined by the model
and the correlation parameters.}
\usage{
simPts(covariates,int.formula=~1, int.par=c(1), EN=100,
                    model, showplot=FALSE, showpts=FALSE, ...)
}
\arguments{
  \item{covariates}{a matrix with columns x,y and any number of covariates
                  x and y are the mid points of the grid cells; the order
                  of the rows must match the formulation for function im}
  \item{int.formula}{formula for deriving expected intensity from covariates}
  \item{int.par}{parameters for intensity formula}
  \item{EN}{expected number of points}
  \item{model}{a spatial model for suimulation. See \code{\link{RFsimulate}} for further specification.}
  \item{showplot}{if TRUE, plot intensity and point process}
  \item{showpts}{if TRUE show points on the plot}
  \item{...}{parameters, if any, passed to plot}
}
\value{ppp object of point locations
}
\author{ Devin Johnson; Jeff Laake}
\seealso{ \code{\link{simCovariates}}}
\examples{
hab.range=30
probs=c(1/3,2/3)
covariates = simCovariates(hab.range, probs)
xpp=simPts(covariates=covariates,int.formula=~factor(habitat),int.par=c(0,1,2))
plot(xpp)
}

