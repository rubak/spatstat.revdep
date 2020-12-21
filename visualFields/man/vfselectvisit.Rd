\name{vfselectvisit}
\alias{vfselectvisit}
\title{select visits per subject from a vf-object}
\description{select a determined amount of visits a vf-object. This is done if the interest is to select the last n visits, or the first n visits, or visits within a range, etc}
\usage{vfselectvisit( vf, sel = "last", numTests = 1,
               beginDate = NA, endDate = NA )}
\arguments{
  \item{vf}{a vf-object}
  \item{sel}{Type of selection, do we want visits from \code{last}, from first \code{first}, within a \code{date} range, or specific visit numbers. Default is \code{last}}
  \item{numTests}{number of tests to select from \code{last}, \code{first}, or a \code{date} range. If \code{sel} is an array of indices or set to a \code{date} range, \code{numTests} is overruled. Default is \code{1}}
  \item{beginDate}{when \code{sel = "date"} is date from (inclusive). If it is \code{NA}, then the data is set to \code{1900-01-01}. Default is \code{NA}}
  \item{endDate}{when \code{sel = "date"} is date to (inclusive). If it is \code{NA}, then the data is set to today. Default is \code{NA}}
}
\value{
returns the subselected visits per subject. For those subject for which at least \code{numTests} visits could not be selected are completely removed
}

\examples{
  vfselectvisit( vf91016left )
}

\author{Ivan Marin-Franch}