\name{loadvfxmlbatch}
\alias{loadvfxmlbatch}
\title{loads visual fields from a set of XML files}
\description{loads visual fields from a set of XML files}
\usage{loadvfxmlbatch( filename, patternMap, typeData = "vf" )}
\arguments{
  \item{filename}{filename in CSV format with the list of XML files to upload and whether the subject whose data is in the XML file is a patient with glaucoma \code{pwg} or a control \code{ctr} }
  \item{patternMap}{pattern of stimulus locations. Default is \code{saplocmap$p24d2}}
  \item{typeData}{Type of data to load; visual field (\code{vf}), total deviations (\code{td}), pattern deviations (\code{pd}), global indices (\code{gi}), visual-field index (\code{vfi}), total-deviation p-values (\code{tdp}), pattern-deviation p-values (\code{pdp}), global indices probability maps (\code{gip}), visual-field-index probability map (\code{vfip}). Default is \code{vf}}
}

\value{returns a \code{\link{vf}}-object with as many rows as XML files exists in \code{filename}}

\seealso{ \code{\link{loadvfcsv}}, \code{\link{loadvfxml}} }

\author{Ivan Marin-Franch}