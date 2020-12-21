\encoding{Latin1}
\name{fdis}
\alias{fdis}
\title{Average Functional Dispersion}
\description{
  Computes average Functional Dispersion for several communities.
}

\usage{
  fdis(x, traits)
}


\arguments{
  \item{x}{A community data matrix containing the abundances of the species in the different communities. Rows are sites and species are columns.}
   \item{traits}{A distance matrix among species in \code{\link[base]{matrix}} format, i.e. obtained using  \code{\link[base]{as.matrix}}(\emph{d}), where \emph{d} is a  \code{dist} object from \code{\link[FD]{gowdis}}, \code{\link[stats]{dist}} or similar functions.}
}

\details{
This function is a wrap to \code{\link[FD]{fdisp}} in package \pkg{FD}. It manages some of the possible problems that could appear when computing automatically functional dispersion for local communities in \link{ifdar} (mainly "empty" communities, common when computing ifdar(r) for very small \emph{r}'s). It is a kind of internal function that wouldn't be usually called by the user.
}

\value{
  Numeric. The average functional dispersion of the communities in matrix \code{x}
}

\seealso{
      \code{\link[FD]{fdisp}} for a description of the functional dispersion measure (FDis) of \enc{Laliberté and Legendre (2010)}{Laliberte and Legendre (2010)}.
}


\references{

   \enc{Laliberté, E. and  Legendre, P. (2010) A distance-based framework for measuring functional diversity from multiple traits. \emph{Ecology} 91, 299-305.}{Laliberte, E. and  Legendre, P. (2010) A distance-based framework for measuring functional diversity from multiple traits. \emph{Ecology} 91, 299-305.} 

}

\author{ Marcelino de la Cruz \email{marcelino.delacruz@urjc.es}}

\keyword{ spatial }
