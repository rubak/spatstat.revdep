% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/getFeatures.R
\name{getFeatures}
\alias{getFeatures}
\alias{getFeatures,ANY-method}
\alias{getFeatures,geom-method}
\alias{getFeatures,Spatial-method}
\alias{getFeatures,sf-method}
\alias{getFeatures,ppp-method}
\alias{getFeatures,Raster-method}
\alias{getFeatures,matrix-method}
\title{Get the table of feature attributes}
\usage{
\S4method{getFeatures}{ANY}(x, ...)

\S4method{getFeatures}{geom}(x, ...)

\S4method{getFeatures}{Spatial}(x, ...)

\S4method{getFeatures}{sf}(x, ...)

\S4method{getFeatures}{ppp}(x, ...)

\S4method{getFeatures}{Raster}(x)

\S4method{getFeatures}{matrix}(x)
}
\arguments{
\item{x}{the object from which to derive the attribute table.}

\item{...}{subset based on logical predicates defined in terms of the
columns in \code{x} or a vector of booleans. Multiple conditions are
combined with \code{&}. Only rows where the condition evaluates to TRUE are kept.}
}
\value{
A table of the feature attributes of \code{x} or an object where the
  features table has been subsetted.
}
\description{
Get tabular information of the attributes of features.
}
\examples{
getFeatures(x = gtGeoms$polygon)

# get a subset of an sf-object
getFeatures(x = gtSF$multilinestring, a == 1)

# get the values of a RasterLayer
getFeatures(x = gtRasters$continuous)
}
\seealso{
Other getters: 
\code{\link{getCRS}()},
\code{\link{getExtent}()},
\code{\link{getGroups}()},
\code{\link{getHistory}()},
\code{\link{getLayer}()},
\code{\link{getName}()},
\code{\link{getPoints}()},
\code{\link{getRes}()},
\code{\link{getType}()},
\code{\link{getWindow}()}
}
\concept{getters}
