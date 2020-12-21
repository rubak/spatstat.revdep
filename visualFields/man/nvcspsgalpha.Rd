\name{nvcspsgalpha}
\alias{nvcspsgalpha}
\docType{data}
\title{IU normative reference values for Contrast Sensitiviy Perimetro. Alpha version!!!}
\description{
Normative reference values or (normative values or nv for short) for Contrast Sensitivy Perimcery (CSP) from IU (William H Swanson and Victor E Malinovsky) databases.
}
\usage{data( nvcspsgalpha )}
\format{This normative-value object is defined by some tables, variables, and sub-structures, some of which are mandatory and some which are not. This nv-object has the following objects:
\describe{
\item{\code{pmapsettings}}{\bold{mandatory}. Information about which percentiles are used for location-wise analysis in this nv and their corresponding color coding for display. This is a table with four columns. First column specifies the percentiles (in percent) for \code{cutoffs} and the other three columns are the corresponding RGB values (defined from 0 to 1) specifying the color code to use for values below that percentile}
\item{\code{globalco}}{\bold{mandatory}. Percentiles to be used for the analysis of globan indices such us mean deviation (MD), pattern standard deviation (PSD) or visual-field index (VFI)}
\item{\code{sgrnfl_zest}}{\bold{mandatory}. The actual normative values for the presentaion algorithm ZEST at the array of test locations SGRNFL}
\item{\code{nvname}}{\bold{mandatory}. Name assigned to this normative-value object}
}
}

\keyword{datasets}

\seealso{\code{\link{nvsapdefault}}}

\author{Ivan Marin-Franch, William H Swanson}