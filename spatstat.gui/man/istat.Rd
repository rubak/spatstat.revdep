\name{istat}
\alias{istat}
\title{Point and Click Interface for Exploratory Analysis of Point Pattern}
\description{
  Compute various summary functions for a point pattern
  using a point-and-click interface.
}
\usage{
 istat(x, xname)
}
\arguments{
  \item{x}{
    The spatial point pattern to be analysed.
    An object of class \code{"ppp"}.
  }
  \item{xname}{
    Optional. Character string to use as the title of the dataset.
  }
}
\value{
  \code{NULL}.
}
\details{
  This command launches an interactive (point-and-click) interface 
  which offers a choice of spatial summary functions
  that can be applied to the point pattern \code{x}.

  The selected summary function is computed for the point pattern
  \code{x} and plotted in a popup window.

  The selection of functions includes
  \code{\link[spatstat.core]{Kest}}, \code{\link[spatstat.core]{Lest}},
  \code{\link[spatstat.core]{pcf}},
  \code{\link[spatstat.core]{Fest}} ,\code{\link[spatstat.core]{Gest}}
  and \code{\link[spatstat.core]{Jest}}.
  For the function \code{\link[spatstat.core]{pcf}} it is possible to control
  the bandwidth parameter \code{bw}.
  
  There is also an option to show simulation envelopes of
  the summary function.
}
\section{Note}{
   Before adjusting the bandwidth parameter \code{bw},
   it is advisable to select \emph{No simulation envelopes}
   to save a lot of computation time.
}
\seealso{
  \code{\link{iplot}}
}
\examples{
   if(interactive()) {
      istat(swedishpines)
   }
}
\author{
  \spatstatAuthors.
}
\keyword{spatial}
\keyword{hplot}
