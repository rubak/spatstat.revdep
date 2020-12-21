\name{create.points.by.offset}
\alias{create.points.by.offset}
\title{ Create point dataframe offset from line}
\description{
For a set of observations with x,y locations on the line and a perpendicular
distance, create a new observation dataframe with true x,y point locations.
}
\usage{
create.points.by.offset(lines, observations)
}
\arguments{
  \item{lines}{
   \preformatted{- data frame of lines with the following structure
             label - unique label
             x0    - x coordinate of beginning of line
             y0    - y coordinate of beginning of line
             x1    - x coordinate of end of line
             y1    - y coordinate of end of line
             width - optional full width of each transect
                      line is in center of transect
             ...   - any number of covariates
       }}
  \item{observations}{
   \preformatted{- data frame of observations with the following structure
             label   - label linking it to a unique line
             x       - x coordinate
             y       - y coordinate
             distance- perpendicular distance;
                       positive=right side; negative=left side
               ...   - any number of covariates
      }}
}
\value{
  observations dataframe with true x,y locations
}
\author{ Jeff Laake}
\seealso{\code{\link{offset.points}}}

