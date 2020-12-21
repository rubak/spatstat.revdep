\name{agecalc}
\alias{agecalc}
\title{Calculates age}
\description{Tool to calculate ages from date of birth and date of test, etc}
\usage{agecalc( from, to, daysyear = NULL )}
\arguments{
  \item{from}{date from}
  \item{to}{date to}
  \item{daysyear}{if \code{NULL}, calculates real age. If informed, then calculates year based on the number of days per year (e.g. 362.25)}
}

\value{returns age in years}

\examples{
  agecalc( "1977-01-31", "2014-01-30" )
}

\author{Ivan Marin-Franch}