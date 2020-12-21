\name{simCovariates}
\alias{simCovariates}
\title{Simulates covariates for an example in DSpat}
\description{
Create a set of covariates in a 100x100 world with a vertical linear feature
 and discrete habitats.}
\usage{
simCovariates(hab.range=30, probs=c(1/3,2/3), river.loc=50)
}
\arguments{
  \item{hab.range}{habitat range that controls patchiness}
  \item{probs}{ordered probablities that define habitat cutoffs}
  \item{river.loc}{x coordinate for north-south river location}
}
\value{  dataframe with columns x,y,river and habitat
}
\details{
The number of habitat types is the length of \code{probs} plus 1.  The
habitats are stored as a numeric from 1 to the number of types, but should
be fitted with habitat as a factor variable.  The distance to the river is
a scaled distance from 0 to 1.
}
\author{ Devin Johnson; Jeff Laake}
\seealso{ \code{\link{simPts}}}
\examples{
covariates = simCovariates(hab.range=50, probs=c(1/3,2/3,7/8))
}

