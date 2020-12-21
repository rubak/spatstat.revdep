\name{loadvfxml}
\alias{loadvfxml}
\title{loads visual fields from a XML file}
\description{loads visual fields from a XML file}
\usage{loadvfxml( filename, patternMap, typeData = "vf", typeSubject = "pwg",
           extractionType = c( "average" ), daysyear = NULL )}
\arguments{
  \item{filename}{filename}
  \item{patternMap}{pattern of stimulus locations. Default is \code{saplocmap$p24d2}}
  \item{typeData}{Type of data to load; visual field (\code{vf}), total deviations (\code{td}), pattern deviations (\code{pd}), global indices (\code{gi}), visual-field index (\code{vfi}), total-deviation p-values (\code{tdp}), pattern-deviation p-values (\code{pdp}), global indices probability maps (\code{gip}), visual-field-index probability map (\code{vfip}). Default is \code{vf}}
  \item{typeSubject}{Type of subject, control (\code{ctr}) or patient (\code{pwg}). Default is (\code{pwg})}
  \item{extractionType}{When typeData is (\code{vf}) what type of extraction we want: all re-tested sensitivities \code{"all"} or just the mean \code{"average"} over re-tested values? Default is \code{"average"}}
  \item{daysyear}{if \code{NULL}, calculates real age. If informed, then calculates year based on the number of days per year (e.g. 362.25)}
}
\details{The XML files format must be as from the extraction of the HFA device}

\value{returns a \code{\link{vf}}-object with one row containing the information for the subject loaded in the XML}

\seealso{ \code{\link{loadvfcsv}}, \code{\link{loadvfxmlbatch}} }

\author{Ivan Marin-Franch}