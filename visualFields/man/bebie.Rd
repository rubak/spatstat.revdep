\name{bebie}
\alias{bebie}
\title{Bebie curve}
\description{Plots Bebie rank TD curve}
\usage{bebie( tdr, type = "conventional", diff = TRUE, percentiles = TRUE,
       correction = TRUE, txtfont = "serif", pointsize = 12, cex = 1 )}
\arguments{
  \item{tdr}{a \code{vf} object with the rank TD curve}
  \item{type}{whether to use a \code{conventional} way to plot the rank TD curve or \code{ghrank} type where the \code{vf} object passed is the reconstructed within-normal TD rank curve. Default is \code{conventional}}
  \item{diff}{add absolute td values or differences from mean normal. Default is \code{TRUE}}
  \item{percentiles}{add percentile lines. Default is \code{TRUE}}
  \item{correction}{add ranked TD line after correcting for general height difference. Default is \code{TRUE}}
  \item{txtfont}{font of the text with visual-sensitivity values. Default is \code{serif}}
  \item{pointsize}{size of the text with visual-sensitivity values. Default is \code{12}}
  \item{cex}{a numerical value giving the amount by which plotting text and symbols should be magnified relative to the default, which is \code{1}}
}

\value{returns a graph with the Bebie rank TD curve (also known as the cumulative defect curve [1]) if the option \code{diff} is set to \code{FALSE}. If the option \code{diff} is set to true, then it returns the difference between the subject's rank curve and the mean normal curve. To use this option (and the options with \code{percentile = TRUE}), the tables \code{nvtdrank}, \code{perctdrank}, and \code{perctdrankadj} must exist in the nv-object, for the perimeter \code{tperimetry}, the pattern of locations \code{talgorithm}, and thepresentation algorithm \code{tpattern}. See \code{\link{nvsapdefault}} to see the structure the tables must have}

\references{
[1] H. Bebie, J. Flammer, and T. Bebie. \emph{The cumulative defect curve: separation of local and diffuse components of visual field damage}. Graefe's Archive Clinical Experimental Ophthalmology, 227, 1989
}

\examples{
  # conventional "Bebie" rank TD curve
  td <- tdval( vf91016right[15,] )
  tdr <- tdrank( td )
  bebie( tdr )
  
  # "ghrank" type of "Bebie" rank TD curve
  td  <- tdval( vf91016right[15,] )
  tdr <- tdrank( td )
  ghr <- ghranktd( td )
  bebie( ghr, type = "ghrank" )
}

\author{Ivan Marin-Franch}