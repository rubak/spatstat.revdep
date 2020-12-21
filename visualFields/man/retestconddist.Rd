\name{retestconddist}
\alias{retestconddist}
\title{Conditional retest distribution}
\description{Computes the conditional retest distribution and the \code{( 1 - alpha / 2 )} conditional retest intervals}
\usage{retestconddist( vf, nbase = 1, nfollow = 1, alpha = 0.1, typequantile = 7 )}
\arguments{
  \item{vf}{Visual field data. It has to have as many visual fields as nbase + nfollow}
  \item{nbase}{Number of visual fields to be used as baseline}
  \item{nfollow}{Number of visual fields to be used as follow up.}
  \item{alpha}{Significance to derive the \code{( 1 - alpha / 2 )} conditional retest intervals. Default value is \code{0.1}}
  \item{typequantile}{An integer between 1 and 9 selecting one of the nine quantile algorithms detailed below to be used. Default value is \code{7}}
}

\value{conditional retest distribution and the \code{( 1 - alpha / 2 )} conditional retest intervals}

\examples{
  vfcondretest <- retestconddist( vf91016left[c(1:4),], nbase = 2, nfollow = 2 )
}

\references{
[1] A. Heijl, A. Lindgren, and G. Lindgren. Test-retest variability in glaucomatous visual fields. American Academy of Ophthalmology, 180, 1989.
}

\author{Ivan Marin-Franch}