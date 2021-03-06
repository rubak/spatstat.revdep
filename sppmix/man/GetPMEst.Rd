% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/get_post.R
\name{GetPMEst}
\alias{GetPMEst}
\title{Retrieve the Surface of Posterior Means}
\usage{
GetPMEst(fit, num_comp = 1, burnin = floor(fit$L/10))
}
\arguments{
\item{fit}{Object of class \code{damcmc_res} or \code{bdmcmc_res}.}

\item{num_comp}{Number of components requested. The posterior will be calculated only
based on the posterior realizations that have this many mixture components. If missing the realizations
corresponding to the MAP number of components are returned. This parameter
is ignored if \code{fit} is of class \code{damcmc_res}.}

\item{burnin}{Number of initial realizations to discard. By default, it is 1/10 of the total number of iterations.}
}
\value{
An object of class \code{intensity_surface}.
}
\description{
The function first calculates the posterior means
of the parameters of the components of the
mixture intensity, based on a DAMCMC or
BDMCMC fit. Then the surface of posterior
means is calculated using the posterior
means of the parameters. For a BDMCMC fit,
the number of components should be
specified, and all realizations with
that number of components are gathered
to calculate the posterior intensity surface.

For examples see

\url{http://faculty.missouri.edu/~micheasa/sppmix/sppmix_all_examples.html
#GetPMEst}
}
\examples{
\donttest{
fit <- est_mix_damcmc(pp = spatstat.data::redwood, m = 3)
post_intsurf <- GetPMEst(fit, burnin = 1000)
plot(post_intsurf)
fit <- est_mix_bdmcmc(pp = spatstat.data::redwood, m = 5)
post_intsurf <- GetPMEst(fit, num_comp = 4, burnin = 1000)
plot(post_intsurf)
post_fixed = FixLS_da(fit,approx=FALSE, plot_result = TRUE)
plot(GetPMEst(post_fixed))}

}
\seealso{
\code{\link{est_mix_damcmc}},\code{\link{est_mix_bdmcmc}}
}
\author{
Jiaxun Chen, Sakis Micheas, Yuchen Wang
}
