\encoding{Latin1}
\name{mitable}
\alias{mitable}

\title{Tabulate Marks in Neighbourhood of Every Point in a Point Pattern}
\description{
  For a sequence of radii defining different neighbourhood sizes, visit each point in a focal point pattern, find the neighbouring points in a target point pattern, and compile a frequency table of the marks of these neighbour points for each radii.
}

\usage{
  mitable(ppp1,ppp2,r)
  }

\arguments{
  \item{ppp1}{"Focal" point pattern (an \code{\link[spatstat]{ppp}} object of \pkg{spatstat}) or a \link[base]{data.frame} with columns named "x" and "y"}
  \item{ppp2}{"Target" multivariate point pattern (an \code{\link[spatstat]{ppp}} object of \pkg{spatstat}) or a \link[base]{data.frame} with columns named "x", "y", and "marks". The column "marks" wears the (species) name of each point}
   \item{r}{Vector with the sequence of radii (\emph{r>0}) that define different neiborhood sizes.}
}

\details{
Given both a "focal" and a "target" point patterns, \code{mitable} visits each point in the focal point pattern, finds the neighbouring points in the target point pattern, and compile a frequency table of the marks of these neighbour points, where the neibourhood is defined by circles orf radius \emph{r} around the focal points. From an ecological point of view, it provides a \enc{\emph{relevé}}{\emph{releve}} of the local community around each focal tree. It repeats the process for each provided \emph{r}.
}


\value{
  A list, with length = \emph{length(r)}. Within each element of the list, a matrix with dimensions \emph{np} X \emph{nsp}, where \emph{np} is the number of points of the focal point pattern and \emph{nsp} is the number of unique species (i.e., unique marks) in the target point pattern. Cell values in the matrix represent number of neighbours of the \emph{j} species for individual \emph{i} at the considered neibourhood size.
}
\seealso{
  \code{\link[spatstat]{marktable}} for a similar function.
}

\author{ Marcelino de la Cruz \email{marcelino.delacruz@urjc.es}}

\examples{

# Compute number of trees of different species  within circles of several  radii r 
# around each individual maple tree in lansing woods

data(lansing)

maple<- unmark(lansing[lansing$marks=="maple"])
maple
r= c(0.05, 0.08, 0.1)

counts <- mitable(ppp1=maple, ppp2=lansing, r=r) 
dim(counts[[1]])
head(counts[[1]])


}

\keyword{ spatial }

