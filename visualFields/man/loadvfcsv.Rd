\name{loadvfcsv}
\alias{loadvfcsv}
\title{load visual fields from a CSV file}
\description{loads visual fields from a CSV file}
\usage{loadvfcsv( filename, patternMap )}
\arguments{
  \item{filename}{filename}
  \item{patternMap}{pattern of stimulus locations}
}
\details{The columns in the CSV file must follow the format of \code{\link{vf}}. Make sure that all dates have the format \code{MM/DD/YYYY} and all times have the format \code{HH:MM:SS}. Excel tends to change the format, and any edits in Excel that are safe may cause problems when loading the data as \code{\link{vf}}-object}

\value{returns a \code{\link{vf}}-object with all the data in the CSV file \code{filename}}

\seealso{ \code{\link{loadvfxml}}, \code{\link{loadvfxmlbatch}} }

\author{Ivan Marin-Franch}