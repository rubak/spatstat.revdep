\name{pisar}
\alias{pisar}
\alias{risar}
\alias{controldis}

\title{
Phylogenetic Individual Species Area Relationship
}
\description{
Compute Phylogenetic Individual Species Area Relationship function, i.e., PISAR(r) and its normalized version rISAR(r).
}
\usage{
pisar(mippp, mippp.sp = NULL, mimark = NULL, namesmark = NULL, d = NULL, r = NULL,
         buffer = 0, bfw = NULL)
risar(mippp, mippp.sp = NULL, mimark = NULL, namesmark = NULL, d = NULL, d0 = NULL,
        r = NULL, buffer = 0, bfw = NULL)
controldis(d, m, mimark)
}


\arguments{
  \item{mippp}{A multitype (a.k.a. multivariate) marked point pattern. An object with the \link[spatstat.geom]{ppp} format of \pkg{spatstat}.}
  \item{mippp.sp}{Univariate point pattern of the focal species. An object with the \link[spatstat.geom]{ppp} format of \pkg{spatstat}.}
  \item{mimark}{Character. Name of the focal species in the multitype \code{mippp}.}
  \item{namesmark}{Character. If the marks in \code{mippp} are within a \code{data.frame}, the name of the column with the species names.}
  \item{buffer}{One of "adapt", i.e., compute an adaptive buffer, or a number indicating the width of a fixed buffer area around the plot border.}
  \item{bfw}{An \link[spatstat.geom]{owin} object indicating the limits of the buffer area.}
  \item{r}{Vector of distances to compute PIDAR(r) functions.}
  \item{d}{A matrix expresing relationships (usually functional or phylogenetic) between species present in the multivariate point pattern.}
  \item{d0}{Another matrix expresing relationships (usually functional or phylogenetic) between species present in the multivariate point pattern.}
  \item{m}{A community (sites x species) data table.}
}
\details{
 The original definition of \emph{ISAR}(\emph{r}) (Wiegand et al. 2007) was reformulated as:

\deqn{ISAR_f (r) =\sum_{m=1}^{S} \delta_{fm} D_{fm} (r)}{PISAR[f](r) = sum[m] delta[f,m]D[f,m](r)}
(Wiegand and Moloney 2014; Wang et al. 2016), where \eqn{D_{fm}(r)}{D[f,m](r)} describe the probabilities that the nearest species \emph{m} neighbor of the typical
individual of the focal species \emph{f} is located within distance \emph{r}, and \eqn{\delta_{fm} }{delta[f,m]} yields a value of zero if \emph{f = m} and a value of one otherwise (note that in their original proposal ISAR was formulated as if the value assigned to \eqn{\delta_{fm} }{delta[f,m]} were 1 for all species pairs, including \emph{f = m}). Based in  this re-formulation, they defined the \emph{Phylogenetic Individual Species Area Relationship}, i.e., \emph{PISAR}(\emph{r}), as:


\deqn{PISAR_f (r) =\sum_{m=1}^{S} \delta_{fm}^{phy}D_{fm} (r)}{PISAR[f](r) = sum[m] delta[f,m]^phy D[f,m] }

where \eqn{\delta_{fm}^{phy}}{delta[f,m]^phy} is an index of phylogenetic (or functional) dissimilarity between species \emph{f}
and \emph{m}. \emph{PISAR}(\emph{r}) quantifies the expected phylogenetic (or functional) diversity of species within the
neighborhood with radius \emph{r} around the typical individual of the focal species \emph{f}.

 They also defined \emph{rISAR}(\emph{r}) as a  a function that is independent of local species richness within the neighborhood \emph{r}; for this, they divided the PISAR function by the ISAR function:

\deqn{rISAR_f (r) = \frac{\sum_{m=1}^{S} \delta_{fm}^{phy}D_{fm} (r)}{\sum_{m=1}^{S} \delta_{fm}D_{fm} (r)}  }{rISAR[f](r) = (sum[m] delta[f,m]^phy D[f,m]) /(sum[m] delta[f,m] D[f,m])}


If the placement of the focal species \emph{f} is unrelated with functional or phylogenetic relationships with their neighbors, the \eqn{rISAR_f(r)}{rISAR[f](r)} will approximate the mean pairwise functional (or phylogenetic) dissimilarity 
\eqn{\Delta_f^P = \sum_m \delta_{fm}^{phy}/(S-1)}{Delta{f}^P = sum{m} delta[f,m]^phy /(S-1)} between an individual of the
focal species \emph{f} and all other species in the plot.


The function \code{controldis} controls that the order of species in the phylogenetic distance matrix matches the order of species amomng the levels of species marks in the point pattern, and extracts the vector of distances from all species to the focal one (\code{mimark}).  


}
\value{
  \code{pisar} and \code{risar} return an object of class \code{fv}, see  \link[spatstat.geom]{fv.object}, which can be plotted directly using  \link[spatstat.core]{plot.fv}. 

Essentially a data frame containing a column named \emph{r} with the vector of values of the argument r at which the PISAR(r) or rISAR(r) function had been estimated and  another column, named \emph{"risar"} or \emph{"pidar"},which  contains an estimate of the selected function.

\code{controldis} returns either the vector of distances between the focal and the rest of species or a vector of 1's if there is not phyloegenetic distance provided.

}

\section{Simulation envelopes}{
To compute simulation envelopes for  \code{pisar} or \code{risar} functions,  use \code{\link[spatstat.core]{envelope}}. See the examples in this help page and in \code{\link[idar]{multifocalsimulator}} to know how to compute simulation envelopes from appropriate null models.
}


\references{
Wang, X., et al. (2016). Stochastic dilution effects weaken deterministic effects of niche-based processes in species rich forests. \emph{Ecology, in press}

Wiegand,T., Gunatilleke, C.V.S., Gunatilleke, I.A.U.N. and  Huth, A. (2007) How individual species structure diversity in tropical forests. \emph{PNAS} 104, 19029-19033.

Wiegand, T., and K.A. Moloney. (2014). \emph{A handbook of spatial point pattern analysis in ecology}. Chapman and Hall/CRC press, Boca Raton, FL 
}
\author{ Marcelino de la Cruz \email{marcelino.delacruz@urjc.es}}


 \section{NOTE}{
 When computing risar it is necessary to provide a phylogenetic or functional distance matrix to the argument \code{d}. By default,  argument \code{d0} will be set to a vector of 1's. It is however possible to provide a different matrix to \code{d0} and compute instead, e.g.,  a ratio of phylogenetic to functional diversity.
 
 }

\seealso{
See also \code{\link{isar}} for other individual diversity area functions.
}
\examples{

data(SF)
data(SFphylotree)

# Discard the size mark and keep the species mark in SF ppp:
sfsp<- ppp(SF$x, SF$y, window=SF$window, marks=SF$marks$species)

# compute phylogenetic distance among species
dphy <- cophenetic(SFphylotree)

# compute and plot PISAR function for sp_44
pisar_44 <- pisar(sfsp, mimark="sp_44", r=1:15,  d=dphy)
plot(pisar_44)

\dontrun{

# Compute rISAR and plot envelopes for an inhomogeneous Poisson model
#  of each species in San Francisco plot 
# BEWARE: THIS TAKES QUITE A FEW MINUTES !!!

# Split sfsp point pattern ppp by species
sfsp.sp<- split(sfsp)

#  Species with >= 10 individuals
sfsp10 <- sapply(sfsp.sp, function(x) x$n>=10)
#names of those species
nombressf<- names(sfsp10[sfsp10])

# parameters for the simulations, estimation of intensity, etc.
nsim<-199
nmin<-10
sigma <- 8
r<- seq(1,15, by=0.5)

# list to store results
risar.sf<- list()
# start computation
for( sp in nombressf){
   print(sp)
   # estimate intensity of the focal species
   lambda<- density(unmark(sfsp[sfsp$marks==sp]), sigma=sigma, positive=TRUE)      
   # obtain simulated patterns were all species ecept the focal remain fixed 
   #  and the focal varies according to an inhomomgeneous Poiisson process 
   simulados<- multifocalsimulator(sfsp, mimark=sp,
              simulate=expression(rpoispp(lambda)), nsim=nsim,nmin=nmin)
   # compute risar
   risar.sf[[sp]] <- envelope(sfsp, risar, mimark=sp, d=dphy, r=r,
        simulate=simulados,nsim=nsim, savefuns=T, buffer=0)
}

# plot the results
 dev.new(height=7, width=16)
 par(mfrow=c(3,9))
 for(i in 1:27) plot(risar.sf[[i]], legend=F, main=nombressf[i])

}


}

\keyword{ spatial }
