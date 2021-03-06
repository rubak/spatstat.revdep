% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/markedIPPP_mix.R
\name{est_MIPPP_cond_loc}
\alias{est_MIPPP_cond_loc}
\alias{plot.MIPPP_fit}
\alias{summary.MIPPP_fit}
\title{Fit a MIPPP conditionally on location}
\usage{
est_MIPPP_cond_loc(pp, r, hyper = 1, L = 10000, burnin = floor(L/10),
  m = 3, fit_groundIPPP = FALSE, truncate = FALSE, grayscale = FALSE,
  startgamma, discrete_mark = TRUE, LL = 150, open_new_window = FALSE,
  show_plots = TRUE)

\method{plot}{MIPPP_fit}(x, surf, open_new_window = FALSE,
  grayscale = FALSE,
  main = "Locations and ground intensity surface of a marked IPPP", ...)

\method{summary}{MIPPP_fit}(object, ...)
}
\arguments{
\item{pp}{Marked point pattern of class \code{ppp}.}

\item{r}{Radius used to define the neighborhood system. Any two locations within this distance are considered neighbors.}

\item{hyper}{Hyperparameter for the proposal distribution of gamma. This is currently the standard deviation for the random walk Metropolis-Hastings steps (one step for each gamma). Use a small value.}

\item{L}{Number of iterations for the MCMC; default is 10000.}

\item{burnin}{Number of initial realizations to discard. By default, it is 1/10 of the total number of iterations.}

\item{m}{The number of components to fit for the ground process when \code{fit_groundIPPP=TRUE}.}

\item{fit_groundIPPP}{Logical variable requesting to fit and return the DAMCMC results of the ground process.}

\item{truncate}{Logical variable indicating whether or not we
we only work with events within the window defined
in the point pattern \code{pp}.}

\item{grayscale}{Logical to request plots in grayscale.}

\item{startgamma}{Initial value for the gamma vector. If missing the zero vector is used.}

\item{discrete_mark}{Logical flag indicating whether the mark is a discrete (numerical value) or not.
For continuous marks set this to FALSE.}

\item{LL}{Length of the side of the square grid.}

\item{open_new_window}{Open a new window for a plot.}

\item{show_plots}{Logical variable requesting to produce the probability field plots for each mark.}

\item{x}{An object of class \code{MIPPP_fit} (the result of a request from \code{\link{est_MIPPP_cond_loc}}).}

\item{surf}{An object of type \code{intensity_surface} representing the
IPPP surface for the ground process (conditioning on location only). This can be the surface of posterior means
or the MAP from a \code{damcmc_res} object or the MAP number of components
surface from a \code{bdmcmc_res} object.}

\item{main}{Main title for the plot.}

\item{...}{Additional arguments for the S3 method.}

\item{object}{An object of class \code{MIPPP_fit} (the result of a request from \code{\link{est_MIPPP_cond_loc}}).}
}
\value{
An object of class \code{MIPPP_fit}, which is simply a list containing the following components:
\item{gen_gammas}{Posterior realizations of the gammas.}
\item{prob_fields}{Probability fields of marks.}
\item{discrete_mark}{Same logical flag as the input argument.}
\item{r}{Same as the input argument.}
\item{pp}{Same as the input argument.}
\item{ground_fit}{An object of type \code{damcmc_res} which contains the results of a DAMCMC fit to the ground process. If \code{fit_groundIPPP=FALSE} this is \code{NULL}.}
\item{condition_on_loc}{Logical variable indicating the type of conditioning used in order to produce this MIPPP fit. For this function it is set to TRUE.}
}
\description{
This function fits a Marked IPPP (MIPPP) on a marked
point pattern by modeling the (joint)
intensity surface of the locations and the marks
using an IPPP for the locations (independent
of the mark values) and for discrete marks a Gibbs model for the
mark distribution which is conditionally defined
on all the locations. NOTE: The estimation procedure for continuous
marks (random fields) will be implemented
in future versions of the \code{sppmix} package.

For examples see

\url{http://faculty.missouri.edu/~micheasa/sppmix/sppmix_all_examples.html
#est_MIPPP_cond_loc}
}
\details{
We assume that the joint distribution of a
marked point pattern \code{N=[s,m(s)]} with \code{n}
events is of the form:

\code{p(N)=lambda^n*exp(-lambda)/(n!)*f(all s|theta1)*g(all m|theta2(s),all s)}

where \code{s} denotes a location and \code{m=m(s)}
a mark value at that location, lambda a parameter
with the interpretation as the average number of points
over the window of observation, and \code{f}, \code{g} are proper densities.

The location (or ground process) \code{f(all s|theta1)}
can be fit using any method for unmarked
point patterns (for us it is modeled using
an IPPP with a mixture of normals
intensity surface). The function fits
the parameters of the second part of this model by default. However,
setting \code{fit_groundIPPP=TRUE} will fit a mixture intensity
surface for the ground process and return it for future processing.
Alternatively, simply retrieve the marked point process returned
and fit the ground process using \code{\link{est_mix_damcmc}}
or  \code{\link{est_mix_bdmcmc}} (the marks will be ignored and only
the locations will be used).

Since \code{s} is observed over some
window and the marks are conditioned on
knowing the locations, then \code{g} is a
random field for each value of \code{m}.

The neighborhood system is controlled by
\code{r} and is crucial in this case. Small values
tend to produce probability fields with concentrated
masses about observed events of the process,
whereas, large neighborhoods allow us to borrow
strength across locations and result in much smoother
probability fields. This parameter is currently NOT estimated
by the \code{\link{sppmix}} package, but will
be implemented in future releases.

See Micheas (2014) for more details on special cases (2 marks) of
these Marked IPPP models via conditioning arguments.
}
\examples{
\donttest{
# Create a marked point pattern
x <- runif(100)
y <- runif(100)
#mark distribution is discrete uniform
m <- sample(1:2, 100, replace=TRUE)
m <- factor(m, levels=1:2)
pp <- spatstat.geom::ppp(x, y, c(0,1), c(0,1), marks=m)
# estimate the probability fields for each mark; since we have a discrete
# uniform for the mark distribution we should see probabilities about .5
# for both marks, as well as, the gamma credible sets should include 0,
# meaning that the marks are independent of location (probability .5 for
# each of the two mark values)
mpp_est <- est_MIPPP_cond_loc(pp, 0.1, hyper=0.2)
GetStats(mpp_est$gen_gammas[,1])$CredibleSet
GetStats(mpp_est$gen_gammas[,2])$CredibleSet
mpp_est <-est_MIPPP_cond_loc(pp, 0.3, hyper=0.2)
GetStats(mpp_est$gen_gammas[,1])$CredibleSet
GetStats(mpp_est$gen_gammas[,2])$CredibleSet
mpp_est <- est_MIPPP_cond_loc(pp, 0.5, hyper=0.2)
GetStats(mpp_est$gen_gammas[,1])$CredibleSet
GetStats(mpp_est$gen_gammas[,2])$CredibleSet
#Visualize the Tornado data about MO. We request to fit both the mark
#and ground process.
#plot the states, the tornado locations and the marks (strength of a tornado)
ret=PlotUSAStates(states=c('Iowa','Arkansas', 'Missouri','Illinois','Indiana',
 'Kentucky','Tennessee', 'Kansas','Nebraska','Texas','Oklahoma', 'Mississippi',
 'Alabama','Louisiana'), showcentroids=FALSE,shownames=TRUE, plotlevels = FALSE,
 main= "Tornadoes about MO, 2011")
#check out the mark values and their frequency
table(Tornadoes2011MO$marks)
#plot each point with a different shape according to its marks
ret$p+ggplot2::geom_point(data=as.data.frame( Tornadoes2011MO),ggplot2::aes(x=x,
 y=y,shape= as.factor(marks)))+ggplot2::guides(shape = ggplot2::guide_legend(
 title="Tornado power", ncol=2,byrow=TRUE))
#plot each point with a different color according to its marks
ret$p+ggplot2::geom_point(data=as.data.frame( Tornadoes2011MO),ggplot2::aes(x=x,
 y=y,color= as.factor(marks)))+ggplot2::guides(color = ggplot2::guide_legend(
 title="Tornado power", ncol=2,byrow=TRUE))
#plot each point with a different circle size according to its marks
ret$p+ggplot2::geom_point(data=as.data.frame( Tornadoes2011MO),ggplot2::aes(x=x,
 y=y,size=marks),shape=21)+ ggplot2::scale_size_continuous(breaks=sort(unique(
 Tornadoes2011MO$marks))) + ggplot2::guides(size =ggplot2::guide_legend(title=
 "Tornado power", ncol=2,byrow=TRUE))
# the marks must start from 1, recode the original
Tornadoes2011MO1=Tornadoes2011MO
Tornadoes2011MO1$marks=Tornadoes2011MO1$marks+1
mpp_est=est_MIPPP_cond_loc(Tornadoes2011MO1,r=1.5,hyper=0.01,
 startgamma = c(.1,.2,.3,.4,.5,.6),fit_groundIPPP=TRUE)
#Now generate an MIPPP with 3 marks
newMPP=rMIPPP_cond_loc(gammas=c(.1,.2,.5))
summary(newMPP)
plot(newMPP$surf,main="True IPPP intensity surface for the locations")
true_gammas=newMPP$gammas
genMPP=newMPP$genMPP
newMPP$r
mpp_est=est_MIPPP_cond_loc(genMPP,newMPP$r, hyper=0.2)
GetStats(mpp_est$gen_gammas[,1])$Mean
GetStats(mpp_est$gen_gammas[,2])$Mean
GetStats(mpp_est$gen_gammas[,3])$Mean
GetStats(mpp_est$gen_gammas[,1])$CredibleSet
GetStats(mpp_est$gen_gammas[,2])$CredibleSet
GetStats(mpp_est$gen_gammas[,3])$CredibleSet
summary(mpp_est)
plot(mpp_est)
plot(mpp_est,newMPP$surf)}

}
\references{
Hierarchical Bayesian Modeling of Marked Non-Homogeneous Poisson Processes with finite mixtures and inclusion of covariate information. Micheas, A.C. (2014). Journal of Applied Statistics, 41, 12, 2596-2615, DOI: 10.1080/02664763.2014.922167.
}
\seealso{
\code{\link{GetStats}},
\code{\link{rMIPPP_cond_loc}}
}
\author{
Sakis Micheas, Jiaxun Chen
}
