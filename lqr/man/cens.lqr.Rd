\name{cens.lqr}
\alias{cens.lqr}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Skew-t quantile regression for censored and missing data
}
\description{
It fits a linear quantile regression model where the error term is considered to follow an SKT skew-t distribution, that is, the one proposed by Wichitaksorn et.al. (2014). Additionally, the model is capable to deal with missing and interval-censored data at the same time. Degrees of freedom can be either estimated or supplied by the user. It offers estimates and full inference. It also provides envelopes plots and likelihood-based criteria for assessing the fit, as well as fitted and imputed values.
}
\usage{
cens.lqr(y,x,cc,LL,UL,p=0.5,nu=NULL,precision=1e-06,envelope=FALSE)
}

\arguments{
  \item{y}{
the response vector of dimension \eqn{n} where \eqn{n} is the total of observations. It may contain both missing and censored values represented by \code{NaN}s.
}
\item{x}{
design matrix for the fixed effects of dimension \eqn{N x d} where \eqn{d} represents the number of fixed effects including the intercept, if considered.
}
  \item{cc}{
vector of censoring/missing indicators. For each observation it takes 0 if non-censored/missing, 1 if censored/missing.
}
  \item{LL}{
the vector of lower limits of dimension \eqn{n}x\eqn{1}. See details section.
}
  \item{UL}{
the vector of upper limits of dimension \eqn{n}x\eqn{1}. See details section.
}
\item{p}{
An unique quantile of interest to fit the quantile regression.
}
\item{nu}{
It represents the degrees of freedom of the skew-t distribution. When is not provided, we use the MLE.
}
  \item{precision}{
The convergence maximum error permitted. By default is 10^-6.
}
  \item{envelope}{
if \code{TRUE}, it will show a confidence envelope for a curve based on bootstrap replicates.
it is \code{FALSE} by default.
}
}
\details{
Missing or censored values in the response can be represented imputed as \code{NaN}s, since the algorithm only uses the information provided in the lower and upper limits LL and UL. The indicator vector \code{cc} must take the value of 1 for these observations.

*Censored and missing data*

If all lower limits are \code{-Inf}, we will be dealing with left-censored data.
Besides, if all upper limits are \code{Inf}, this is the case of right-censored data. Interval-censoring is considered when both limits are finites. If some observation is missing, we have not information at all, so both limits must be infinites.

Combinations of all cases above are permitted, that is, we may have left-censored, right-censored, interval-censored and missing data at the same time.
}
\value{
\item{iter}{number of iterations.}
\item{criteria}{attained criteria value.}
\item{beta}{fixed effects estimates.}
\item{sigma}{scale parameter estimate for the error term.}
\item{nu}{Estimate of \code{nu} parameter detailed above.}
\item{SE}{Standard Error estimates.}
\item{table}{Table containing the inference for the fixed effects parameters.}
\item{loglik}{Log-likelihood value.}
\item{AIC}{Akaike information criterion.}
\item{BIC}{Bayesian information criterion.}
\item{HQ}{Hannan-Quinn information criterion.}
\item{fitted.values}{vector containing the fitted values.}
\item{imputed.values}{vector containing the imputed values for censored/missing observations.}
\item{residuals}{vector containing the residuals.}
}
\author{Christian E. Galarza <\email{chedgala@espol.edu.ec}>, {Marcelo Bourguignon} <\email{m.p.bourguignon@gmail.com}> and
Victor H. Lachos <\email{hlachos@ime.unicamp.br}>

Maintainer: Christian E. Galarza <\email{chedgala@espol.edu.ec}>
}


\examples{
\dontrun{
##Load the data
data(ais)
attach(ais)

##Setting
y<-BMI
x<-cbind(1,LBM,Sex)

cc = rep(0,length(y))
LL = UL = rep(NA,length(y))

#Generating a 5\% of interval-censored values
ind = sample(x = c(0,1),size = length(y),replace = TRUE,prob = c(0.95,0.05))
ind1 = (ind == 1)

cc[ind1] = 1
LL[ind1] = y[ind1] - 10
UL[ind1] = y[ind1] + 10
y[ind1] = NA #deleting data

#Fitting the model

# A median regression with unknown degrees of freedom
out = cens.lqr(y,x,cc,LL,UL,p=0.5,nu = NULL,precision = 1e-6,envelope = TRUE)

# A first quartile regression with 10 degrees of freedom
out = cens.lqr(y,x,cc,LL,UL,p=0.25,nu = 10,precision = 1e-6,envelope = TRUE)
}
}

\references{
Galarza, C., Lachos V.H. and Bourguignon M. (2020). A skew-t quantile regression for censored and missing data. STAT (under review).

Galarza, C., Lachos, V. H., Cabral, C. R. B., & Castro, C. L. (2017). Robust quantile regression using a generalized class of skewed distributions. Stat,
6(1), 113-130.

Wichitaksorn, N., Choy, S. T., & Gerlach, R. (2014). A generalized class of skew distributions and associated robust quantile regression models. Canadian Journal of Statistics, 42(4), 579-596.
}

\keyword{package}
\keyword{quantile}
\keyword{regression}
\keyword{skew}

\seealso{\code{\link{lqr}},\code{\link{best.lqr}},\code{\link{Log.lqr}},
\code{\link{Log.best.lqr}},\code{\link{dSKD}},\code{\link[qrLMM]{QRLMM}},
\code{\link[qrNLMM]{QRNLMM}}}
