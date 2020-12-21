\name{poplr_cstat}
\alias{poplr_cstat}
\title{permutation of pointwise linear regression (PoPLR): calculation of the combined statistic}
\description{For details see [1]}
\usage{poplr_cstat( pval, truncVal = 1 )}
\arguments{
  \item{pval}{p-values. Typically the ones obtained from \code{\link{poplr_pstat}}}
  \item{truncVal}{p-value cut-off for truncation. Default is \code{1}}
}

\seealso{\code{\link{poplr}}, \code{\link{poplr_pstat}}, \code{\link{vflayout_poplr}}, \code{\link{vfplot_plr}}}

\references{
[1] N. O'Leary, B. C. Chauhan, and P. H. Artes. \emph{Visual field progression in glaucoma: estimating the overall significance of deterioration with permutation analyses of pointwise linear regression (PoPLR)}. Investigative Ophthalmology and Visual Science, 53, 2012
}

\examples{
  res <- poplr( vf91016right )
}

\author{Ivan Marin-Franch}