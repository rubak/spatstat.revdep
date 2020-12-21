\name{pval}
\alias{pval}

\title{
P Value for a  Monte Carlo Test of Polygon Overlapping
}
\description{
Computes the p-value and sign of deviation for the hypothesis that  the first value in a vector  is larger or smaller (i.e., a two-sided test)  than the expected value represented by all the other values in the vector, or alternatively, computes a one-sided test.
}
\usage{
pval(x, alternative=c("two.sided", "less", "greater"))
}
\arguments{
  \item{x}{A vector, usually with some observed statistic in the first position (the observed overlapp) followed by a sequence of the same statistic computed for several realizations of a null model against which we test our hypothesis (i.e., a sequence of simulated, i.e., rotated, overlapps).
}
 \item{alternative}{a character string specifying the alternative hypothesis, must be one of "two.sided" (default), "greater" or "less". You can specify just the initial letter.}
}
\value{
Vector of length one whose absolute value represents the p-value and whose sign indicates wether the first value in \code{x} is larger (positive) or smaller (negative) than the expected value.
}

\author{
Marcelino de la Cruz
}
\examples{
pval(c(0,1:99))
pval(c(100,1:99))
pval(c(100,1:199))
pval(c(200,1:199))
pval(c(0,1:199))
}
\keyword{ distribution }

