\name{vfenv}
\alias{vfenv}
\docType{data}
\title{environment with the current normative values used \code{visualFields}}
\description{
environment containing the current version of normative values to be used with visualFields. The default used is \code{\link{nvsapdefault}}
}
\usage{data(vfenv)}
\format{It has as many element as the version of normative values used plus one with the name of the version used. To set a structure as the new dataset use \code{\link{setnv}}. The structure has to be similar to \code{\link{nvsapdefault}}.
  \describe{
    \item{\code{...}}{same variables as in the corresponding normative values version used. See default value of \code{nv$nvname}}
    \item{\code{nv$nvname}}{char. Name of the variable with the normative values to use. Default is \code{"nvsapdefault"}}
  }
}

\seealso{\code{\link{getnv}},\code{\link{setnv}}}

\keyword{datasets}

\author{Ivan Marin-Franch}