\name{vfsort}
\alias{vfsort}
\title{sort vf-objects}
\description{sorts a vf-object by id, eye tested, and date and time of test}
\usage{vfsort( vf, decreasing = FALSE )}
\arguments{
  \item{vf}{a vf-object}
  \item{decreasing}{logical. Should the sort order be increasing or decreasing?. Default is \code{FALSE}}
}
\value{
return a sorted vf
}

\examples{
  vfsort( vf91016left, decreasing = TRUE )
}

\author{Ivan Marin-Franch}