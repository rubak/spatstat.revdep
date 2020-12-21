\name{vfobject}
\alias{vf}
\title{visualField objects}
\description{This is the main object of the \code{\link{visualFields}} package. It is essentially a dataframe, but with a fixed number of columns (with pre-determined names) for information about the subject and test data and a variable number of columns for the perimetry results. These can be the sensitivities, or total-deviation values, or pattern-deviation values obtained from static automated perimetry (SAP), frequency-doubling perimetry (FDP), or any other perimetry device. (The number of columns for tested locations is variable as is different for different testing patterns, 24-2, 30-2, 10-2, etc.) Mean deviation, pattern standard deviation, vfi, etc are stored too in a visualField-type object}
\details{The fixed columns of the visualField object with information about subject and test are:
  \tabular{ll}{
    \bold{id}         \tab subject identification number\cr
    \bold{tperimetry} \tab test perimetry. The type of perimetry analysis performed. Possible values include\cr
                      \tab \code{"sap"} and \code{"fdp"}. The value of this column, \code{tperimetry}, is used\cr
                      \tab in conjunction with the value in \code{talgorithm}, and \code{tpattern} to find the\cr
                      \tab corresponding normative values (see help on \code{nv}) to use for data analysis (e.g.\cr
                      \tab calculation of total-deviation and pattern-deviation values and probability maps). At\cr
                      \tab the moment, only normative values for SAP, 24-2, SITA standard, is distributed with\cr
                      \tab \code{\link{visualFields}}. Nevertheless, \code{\link{visualFields}} contains a number of functions\cr
                      \tab that can be used for the generation of normative values (see \code{getnv},\cr
                      \tab \code{ageLinearModel}, \code{sdnv}, \code{tdval}, \code{pdval}, \code{locperc},\cr
                      \tab \code{vfstats}, \code{vfindex}, \code{gloperc}, \code{vfiperc}, \code{setnv}).\cr
    \bold{talgorithm} \tab test algorithm. The algorithm used for the perimetric test. Posible values are\cr
                      \tab \code{sitas} and \code{zest}. At the moment, only normative values for SAP, 24-2, SITA\cr
                      \tab standard, is distributed with \code{\link{visualFields}}\cr
    \bold{tpattern}   \tab test pattern. The pattern of locations used for the perimetric test. Posible values\cr
                      \tab are \code{p24d2} or \code{p10d2}. At the moment, only normative values for SAP, 24-2,\cr
                      \tab SITA standard, is distributed with \code{\link{visualFields}}\cr
    \bold{tdate}      \tab test date\cr
    \bold{ttime}      \tab test time\cr
    \bold{stype}      \tab type of subject. Values can be \code{ctr} for controls, \code{pwg} for patients with\cr
                      \tab glaucoma, \code{sus} for suspect subjects. This is just for information to display in\cr
                      \tab the printouts\cr
    \bold{sage}       \tab subject age. Important for the calculation of total-deviation values and probabiliby\cr
                      \tab maps.\cr
    \bold{seye}       \tab eye tested\cr
    \bold{sbsx}       \tab estimated x-position of the blind spot in degrees of angle of vision\cr
    \bold{sbsy}       \tab estimated y-position of the blind spot in degrees of angle of vision\cr
    \bold{sfp}        \tab false positives\cr
    \bold{sfn}        \tab false negatives\cr
    \bold{sfl}        \tab fixation losses\cr
    \bold{sduration}  \tab total duration of the test\cr
    \bold{spause}     \tab total time of pause
  }
  The reminder of the columns can be different things. For threshold sensitivity values, and total-deviation and pattern-deviation values, and their corresponding probability maps, they are:
  \tabular{ll}{
    \bold{L1 .. L54 .. L68 .. L76} \tab location number. There are up to 54 locations for the 24-2, up to 68 for\cr
                                   \tab the 10-2, and 76 for the 30-2. Information about the position of the\cr
                                   \tab locations, the size of the stimulus, and the x and y coordinates in\cr
                                   \tab degrees of visual angles are specified in \code{saplocmap} (for SAP)\cr
                                   \tab \code{fdplocmap} (for FDP)
  }
  For statistical values of the visual-fields results (mean deviation, pattern standard deviation, and others) and their corresponding probability mapped value, they are:
  \tabular{ll}{
    \bold{msens} \tab mean sensitivity value; or the probability mapped value\cr
    \bold{ssens} \tab standard deviation of the sensitivity values; or the probability mapped value\cr
    \bold{mtdev} \tab mean deviation (mean value of the total-deviation values; or the probability mapped value)\cr
    \bold{stdev} \tab standard deviation of the total-deviation values; or the probability mapped value\cr
    \bold{mpdev} \tab mean value of the pattern-deviation values; or the probability mapped value)\cr
    \bold{stdev} \tab standard pattern deviation (standard deviation pattern-deviation values; or the\cr
                 \tab probability mapped value
  }
  For visual field index (VFI) value and the corresponding probability mapped value, they are:
  \tabular{ll}{
    \bold{mvfi} \tab visual field indes (VFI); or the probability mapped value\cr
    \bold{svfi} \tab standard deviation of the VFI at each location; or the probability mapped value\cr
  }

}

\examples{
# DO NOT EXECUTE
# one can load sensitivities using loadvfcsv or loadvfxml the data so
# vf <- loadvfcsv( filename = "foo.csv", , patternMap = saplocmap$p24d2 )
# calculate total deviation values using \code{\link{visualFields}} normative values for
# SAP SITAS 24-2 (and Goldman size III stimulus)
# td <- tdval( vf )
# calculate pattern deviation values using total deviation values SAP SITAS 24-2
# pd <- tdval( td )
# OR
# pd <- tdval( tdval( vf ) )
# calculate total deviation proabiliby maps
# tdp <- tdpmap( td )
# calculate pattern deviation proabiliby maps
# pdp <- pdpmap( pd )
}

\seealso{\code{\link{vfsettings}}}

\author{Ivan Marin-Franch}