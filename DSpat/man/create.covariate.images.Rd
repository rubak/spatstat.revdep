\name{create.covariate.images}
\alias{create.covariate.images}
\title{Create a list of covariate images}
\description{Creates a list of covariates images from a dataframe of covariates
 defined on a grid for the study area.  Images are created for
 variables contained in vector of names \code{varnames} and the values of the covariates are
 in the \code{covariates} dataframe.
}
\usage{
create.covariate.images(covariates, varnames)
}
\arguments{
  \item{covariates}{covariate dataframe (see DSpat for structure)}
  \item{varnames}{vector of names of fields contained in covariates that will be used for fitted model}
}
\value{covariate.im - list of covariate images (class im) }
\author{Jeff Laake}

