\name{quadscheme.lt}
\alias{quadscheme.lt}
\title{Create line transect quadrature for spatstat}
\description{Creates a quadrature for spatstat from a study area, observations, and lines
}
\usage{
quadscheme.lt(study.area, observations, lines, width = NULL,
                      epsvu = c(1, 0.01), show.warnings=FALSE)
}
\arguments{
  \item{study.area}{owin class giving the boundaries of the study area}
  \item{observations}{
  \preformatted{- data frame of observations with the following structure
             label   - label linking it to a unique line
             x       - x coordinate
             y       - y coordinate
             distance- perpendicular distance from center line
             ...     - any number of covariates
           }}
  \item{lines}{
  \preformatted{- data frame of lines with the following structure
            label - unique label
            x0    - x coordinate of beginning of line
            y0    - y coordinate of beginning of line
            x1    - x coordinate of end of line
            y1    - y coordinate of end of line
            width - optional full width of transect around line
            angle - angle of rotation to get to vertical
            ...   - any number of covariates
           }}
  \item{width}{ if no width field is given in lines then it must be specified here
                as a constant width for all lines}
  \item{epsvu}{pixel dimensions epsvu[1] in v direction (height) and epsvu[2] in
                u direction (width) (these are used once line is rotated vertically.}
  \item{show.warnings}{if TRUE, warnings from quadrature construction will be
                       shown.}
}
\value{
  \item{Q}{quadscheme as defined in spatstat}
  \item{transects}{list of transect polygon dataframes}
  \item{lines.psp}{line segment process}
   }
\author{ Devin Johnson; Jeff Laake}
\seealso{\code{\link{LTDataFrame}}}

