\name{helper functions}
\alias{clean.up.data}
\alias{freqs}
\alias{minusID}
\alias{shake}


\title{Functions for the aid of segregation measures}

\description{
Small functions included in package \code{spatialsegregation}, used for manipulation of forest datasets which
have dbh-values (pp with an element $dbh).
}

\usage{
 clean.up.data(pp, dbh = 10, atleast = 10)
 freqs(pp)
 minusID(pp, minusR, dbh, atleast=0)
 shake(pp, a = 0.001)
}
\arguments{
\item{pp}{Multitype point pattern (see package 'spatstat')}
\item{atleast}{Include specii with abundance atleast \code{atleast}.}
\item{dbh}{Include only those points with dbh atleast \code{dbh}.}
\item{minusR}{Range from the border withing which to exluce points (used for correction of estimates).}
\item{a}{Size of displacement: x+Unif(-a,a), y+Unif(-a,a). }
}

\details{
Small functions to manipulate multitype point patterns.

\code{clean.up.data}: Returns a subsample fullfilling the given constrains.

\code{freqs}: Returns the abundance vector.

\code{minusID}: Returns a 0-1-vector indicating inclusion in a simple minus-correction.

\code{shake}: Shakes the pattern, i.e. adds a random displacement shift to each point.
}

\author{
Tuomas Rajala
tuomas.rajala@iki.fi
}

\keyword{ spatial }
