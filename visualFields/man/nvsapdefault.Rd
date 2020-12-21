\name{nvsapdefault}
\alias{nvsapdefault}
\docType{data}
\title{SUNY-IU normative reference values for static automated perimetry}
\description{
Normative reference values or (normative values or nv for short) for static automated perimetry (SAP) obtained from a combination of control subjects from SUNY and IU databases.
}
\usage{data( nvsapdefault )}
\format{This normative-value object is defined by some tables, variables, and sub-structures, some of which are mandatory and some which are not. Importantly an age linear model must be included for test pattern (24-2) and presentation algorithm (e.g. SITA standard). Also important are the tables specifying the settings. Think of the normative-value object as a structure with several levels. The first level has three main substructures with settings:
\describe{
\item{\code{pmapsettings}}{\bold{mandatory}. Information about which percentiles are used for location-wise analysis in this nv and their corresponding color coding for display. This is a table with four columns. First column specifies the percentiles (in percent) for \code{cutoffs} and the other three columns are the corresponding RGB values (defined from 0 to 1) specifying the color code to use for values below that percentile}
\item{\code{globalco}}{\bold{mandatory}. Percentiles to be used for the analysis of globan indices such us mean deviation (MD), pattern standard deviation (PSD) or visual-field index (VFI)}
\item{\code{nvname}}{\bold{mandatory}. Name assigned to this normative-value object}
From here age linear models, percentile values, etc are defined for combination of test patterns and presentation algorithms. It is important that the construction of the name for each of the normative-value for the particular test and presentation is in agreement with the values of \code{tperimetry} and \code{talgorithm} in \code{\link{vf}}. This is very important because it is the way \code{\link{visualFields}} can identify automatically which normative values are to be used with data for a particular visual field. For instance, for the test pattern 24-2 (\code{p24d2}) and the stimulus-presentation algorithm SITA standard (\code{sitas}), the substructure in the normative-values object must be called  \code{p24d2_sitas}. This substructure itself is composed of other strutures and tables. These are
\item{\code{demographics}}{\bold{optional}. It has stats about the controls subjects used for the calculation of the normative values}
\item{\code{agelm}}{\bold{mandatory}. A table with two columns, intercept and slope, specifying a linear model modeling the (linear) decrease at each location of the sensitivities in decibels (dB) per year. These were calculated with the \code{\link{ageLinearModel}}. See locations 26 and 35 have NAs. That is because these are the locations that correspond anatomically with the blind spot (see \code{\link{vfsettings}}) and are hence  excluded from any analysis}
\item{\code{sds}}{\bold{mandatory}. A table with three columns with the standard deviations for threshold sensitivities, total-deviation values (TD), and pattern-deviation values (PD)}
\item{\code{TDpercloc}}{\bold{mandatory}. A table with as many columns as cutoff percentiles (rows) were defined in \code{pmapsettings}. Each column is the TD cutoff value for each percentile. It is used to get the probability map}
\item{\code{PDpercloc}}{\bold{mandatory}. A table with as many columns as cutoff percentiles (rows) were defined in \code{pmapsettings}. Each column is the PD cutoff value for each percentile. It is used to get the probability map}
\item{\code{percglo}}{\bold{mandatory}. A table with as many columns as cutoff percentiles (rows) were defined in \code{globalco} for global indices. Each row has a different statistical index. Most of them are not really used by convention. The statistical indices considered are mean and standard deviation of the threshold sensitivities (\code{msens} and \code{ssens}), mean and standard deviation of the TD values (\code{mtdev} and \code{stdev}), and mean and standard deviation of the PD values (\code{mpdev} and \code{spdev}). The two indices that are used by convention are \code{mtdev} and \code{spdev}, that is mean deviation (MD) and pattern standard deviation (PSD)}
\item{\code{percvfi}}{\bold{mandatory}. A table with as many columns as cutoff percentiles (rows) were defined in \code{globalco} for global indices. This table has data only for the VFI (\code{mvfi}) and the standard deviation of the VFI at each locaton (\code{svfi})}
\item{\code{nvtdrank}}{\bold{optional}, but necessary to run \code{\link{bebie}} with the option \code{diff = TRUE}. A table with two columns, mean normal Bebie TD rank curve and the standard deviation at each ranked location.}
\item{\code{perctdrank}}{\bold{optional}, but necessary to run \code{\link{bebie}} with the option \code{percentiles = TRUE}. A table with as many columns as cutoff percentiles (rows) were defined in \code{pmapsettings}. Each column is the TD cutoff value for each percentile}
\item{\code{perctdrankadj7}}{\bold{optional}, but necessary to run \code{\link{bebie}} with the option \code{diff = TRUE} and \code{percentiles = TRUE}. Same as \code{perctdrank} but the difference from mean normal Bebie TD rank curve analysis}
  }
}

\references{
[1] ADD REFERENCE FOR 24-2

[2] H. J. Wyatt, M. W. Dul, and W. H. Swanson. \emph{Variability of visual field measurements is correlated with the gradient of visual sensitivity}. Vision Research, 47, 2007.

[3] A. Shafi, W. H. Swanson, and M. W. Dul. \emph{Structure and Function in Patients with Glaucomatous Defects Near Fixation}. Optometry and Vision Science, 88, 2011.

}

\keyword{datasets}

\author{Ivan Marin-Franch, William H Swanson}