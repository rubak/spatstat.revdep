\name{xmlvfval}
\alias{xmlvfval}
\title{extracts visual-field sensitivity values}
\description{extracts visual-field sensitivity values}
\usage{xmlvfval( xmllines, patternMap, extractionType = c( "average" ) )}
\arguments{
  \item{xmllines}{lines from loaded XML files}
  \item{patternMap}{pattern of stimulus locations. Default is \code{saplocmap$p24d2}}
  \item{extractionType}{what type of extraction we want: all re-tested sensitivities \code{"all"} or just the mean \code{"average"} over re-tested values? Default is \code{"average"}}
}

\seealso{ \code{\link{loadvfxml}}, \code{\link{xmldevval}} }

\author{Ivan Marin-Franch}