\name{poplr_pstat}
\alias{poplr_pstat}
\title{Permutation of Pointwise Linear Regression (PoPLR): calculation of the p-value}
\description{For details see [1]}
\usage{poplr_pstat( vf, porder, sltest = NULL )}
\arguments{
  \item{vf}{visual-field data. It can be td or pd as well}
  \item{porder}{order of permunations of visual-fields sensitivities. Each row contains a permutation of tests from vf}
  \item{sltest}{values for the 1-tailed hypothesis test for each location. The reference values are not restricted, but they should be either zero (was there any progression?) or negative (was the progression greater than test value?). Default is \code{NULL}}
}
\value{the function retunts two different structures depending on whehter the analysis is linear regression (\code{type = "slr"} or Spearman rank correlation \code{type = "rank"}). For \code{slr} analysis, the sturcture consists of four matrices with data: \code{pval} (p-value at each permutation and location of the significance of the 1-tailed hypothesis test specified by \code{sl_test} for each location), \code{se} (standard error), \code{sl} (slope), and \code{int} (intercept). For \code{rank} analysis, the sturcture consists of two matrices: \code{pval} (p-value at each permutation and location of the significance of the 1-tailed hypothesis test specified by \code{sl_test} for each location) and \code{rho} (the Spearman rank correlation)}

\seealso{\code{\link{poplr}}, \code{\link{poplr_cstat}}, \code{\link{vflayout_poplr}}, \code{\link{vfplot_plr}}}

\references{
[1] N. O'Leary, B. C. Chauhan, and P. H. Artes. \emph{Visual field progression in glaucoma: estimating the overall significance of deterioration with permutation analyses of pointwise linear regression (PoPLR)}. Investigative Ophthalmology and Visual Science, 53, 2012
}

\examples{
  res <- poplr( vf91016right )
}

\author{Ivan Marin-Franch}