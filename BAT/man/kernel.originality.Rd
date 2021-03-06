% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/BAT.R
\name{kernel.originality}
\alias{kernel.originality}
\title{Functional originality of observations (species or individuals) in a n-dimensional hypervolume representing a given species or community.}
\usage{
kernel.originality(
  comm,
  trait,
  method = "gaussian",
  abund = TRUE,
  frac = 0.1,
  relative = FALSE,
  ...
)
}
\arguments{
\item{comm}{A 'Hypervolume' object constructed with the hypervolume R package or a sites x species matrix, with incidence or abundance data about the species in the community. Note that the use of 'HypervolumeList' object is not implemented for this function yet.}

\item{trait}{A matrix of traits for each species in comm (a species for each row and traits as columns). Must be provided only if 'comm' is a sites x species matrix.}

\item{method}{Method for constructing the 'Hypervolume' object. One of "box" (box kernel density estimation), "gaussian" (Gaussian kernel density estimation), or "svm" (one-class support vector machine). See respective functions of the hypervolume R package for details. Must be provided only if 'comm' is a sites x species matrix. Default is 'gaussian'.}

\item{abund}{A boolean (T/F) indicating whether abundance data should be used (TRUE) or converted to incidence (FALSE) before analysis. If not specified, default is TRUE. Ignored if a Hypervolume is provided as input data.}

\item{frac}{A value between 0.01 and 1, indicating the fraction of random points to be used in the estimation of originality. Default is 0.1.}

\item{relative}{A boolean (T/F) indicating whether originality should be relative to the most original species.}

\item{...}{further arguments to be passed for other methods in hypervolume package.}
}
\value{
A matrix with the originality values of each species or individual in each site.
}
\description{
Average dissimilarity between a species or individual and a sample of random points within the boundaries of the n-dimensional hypervolume.
}
\details{
A measure of the originality (sensu Pavoine et al., 2005) of each observation (species or individuals) used to construct the n-dimensional hypervolume. In a probabilistic hypervolume, originality is calculated as the average distance between each observation to a sample of stochastic points within the boundaries of the n-dimensional hypervolume (Mammola & Cardoso, 2020).
Originality is a measure of functional rarity (sensu Violle et al., 2017; Carmona et al., 2017) that allows to map the contribution of each observation to the divergence components of FD (Mammola & Cardoso, 2020).
The number of sample points to be used in the estimation of the originality is controlled by the frac parameter. Increase frac for less deviation in the estimation, but mind that computation time also increases. For large sample sizes, computation time can be very high (use method = 'box' for a quicker estimation).
If abundance data of species in the community are provided as input data (abund = TRUE), each species trait is weighted by replicating it by the abundance in the estimation of the hypervolume.
}
\examples{
comm <- rbind(c(0,3,2,1), c(1,5,6,2), c(0,0,2,1))
rownames(comm) <- c("Community_1","Community_2","Community_3")
colnames(comm) <- c("Sp_1","Sp_2","Sp_3","Sp_4")

trait <- cbind(c(2.2,4.4,6.1,8.3),c(0.5,1,0.5,0.4),c(0.7,1.2,0.5,0.4))
rownames(trait) <- c("Sp_1","Sp_2","Sp_3","Sp_4")
colnames(trait) <- c("Trait_1","Trait_2","Trait_3")

#Example with community and trait matrices as input data 
#kernel.originality(comm = comm, trait = trait, method='gaussian', abund = TRUE, frac = 0.01)

#Example with hypervolume as input data
#kernel.originality(comm = hypervolume_gaussian(trait))
}
\references{
Carmona, C.P., de Bello, F., Sasaki, T., Uchida, K. & Partel, M. (2017) Towards a common toolbox for rarity: A response to Violle et al. Trends in Ecology and Evolution, 32: 889-891.

Mammola, S. & Cardoso, P. (2020) Functional diversity metrics using kernel density n-dimensional hypervolumes. Methods in Ecology and Evolution, 11: 986-995.

Pavoine, S., Ollier, S. & Dufour, A.-B. (2005) Is the originality of a species measurable? Ecology Letters, 8: 579-586.

Violle, C., Thuiller, W., Mouquet, N., Munoz, F., Kraft, N.J.B., Cadotte, M.W., ... & Mouillot, D. (2017) Functional rarity: the ecology of outliers. Trends in Ecology and Evolution, 32: 356-367.
}
