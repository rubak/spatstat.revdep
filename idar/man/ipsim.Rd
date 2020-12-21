\encoding{UTF-8}
\name{ipsim}
\alias{ipsim}
\alias{ipsimlist}
\alias{simulador2}
\alias{multifocalsimulator}

\title{Simulate Multivariate Point Patterns}
\description{
  These functions simulate  multivariate point patterns from a variety of null models, in the way required to test IDAR(r) functions.
}

\usage{
  ipsim(pp, mimark,sigma=0, lambda=NULL, namesmark=NULL)
  ipsimlist(pp, mimark, listsim)
  simulador2(mimark, milambda, nsim=99)
  multifocalsimulator(pp, mimark, simulate,nsim=99,nmin=NULL)
}

\arguments{
   \item{pp}{A multitype (a.k.a. multivariate) marked point pattern. An object with the \link[spatstat.geom]{ppp} format of \pkg{spatstat}}
   \item{mimark}{Character. Name of the focal species in the multitype \code{pp}.}
   \item{namesmark}{Character. If the marks in \code{pp} are within a \code{data.frame}, the name of the column with the species names}
   \item{sigma}{Sigma for the Gaussian kernel to estimate the intensity of the point pattern to simulate}
   \item{listsim}{List with simulated point paterns from \code{simulador2}}
   \item{lambda}{intensity surface, e.g., an image from density.ppp or predict.ppm}
   \item{milambda}{intensity surface, e.g., an image from density.ppp or predict.ppm}
   \item{nsim}{number of simulations}
   \item{nmin}{expected minimum number of points in each simulated pattern}
    \item{simulate}{either a list of pre-computed univariate point patterns or an expression in the R language indicating how to simulate the patterns.}
 }

\details{
This functions produce simulated point patterns appropriate to to compute envelopes of IDAR(r) functions. The usual tests of IDAR(r) functions require that the multivariate (i.e. multispecies) pattern remains fixed, except for the focal species, that is simulated acording to, e.g., an (inhomogeneous) Poisson process.

 \code{ipsim} returns the multivariate \code{pp} pattern with all species in the same locations except the "focal" one (i.e., the one indicated by the argument \code{mimark}) that is simulated using \link[spatstat.core]{rpoispp}. If  an intensity surface (argument \code{lambda}) is provided, the focal species will be simulated from this surface. If no \code{lambda} is provided but the argument \code{sigma} is >0, an intensity surface will be estimated with a Gaussian kernel with the sigma provided (using \link[spatstat.core]{density.ppp}) and the simulation will be made form this surface. If no lambda is provided and \code{sigma}=0, a homogeneos Poisson process will be simulated.

 \code{simulador2} generates a list (length = \emph{nsim}) of marked (with mark = \emph{mimark}) univariate point patterns from an intensity surface using \link[spatstat.core]{rpoispp}.

 \code{ipsimlist} uses the results of \code{simulador2} and the multivariate \code{pp} pattern to generate a list of multivariate point patterns with all species in the same locations except the "focal" one (i.e., the one indicated by the argument \code{mimark}) that has the locations simulated with \code{simulador2}.
 
  \code{multifocalsimulator} is more flexible and allows the simulation of whichever null model of the focal species that could be described by an R expression.

}


\value{
  \code{ipsim} produces a multivariate point pattern (with the \link[spatstat.geom]{ppp} format of \pkg{spatstat}); \code{ipsimlist}  and \code{multifocalsimulator} produce a list of multivariate point patterns; \code{simulador2} produces a list of univariate marked patterns.
}

\author{ Marcelino de la Cruz \email{marcelino.delacruz@urjc.es}}
\examples{


# Build a multivariate point pattern where maple is simulated according to a Poisson process
# and where all the other species are keep fixed in their original coordinates.
# (The warnings are because in the original lansing point pattern there is a dupplicated
#  point)

data(lansing) 
ipsim(pp=lansing, mimark="maple")

# Build a multivariate point pattern where maple is simulated according to an Inhomogeneous
# Poisson process from an intensity surface estimated "on the fly" with a Gaussian kernel with
# sd = "sigma", and where all the other species are keep fixed in their original coordinates.

ipsim(pp=lansing, mimark="maple", sigma=0.1)

# Build a multivariate point pattern where maple is simulated according to an Inhomogeneous
# Poisson process from a predefined  intensity surface "lambda" and where all the other
# species are keep fixed in their original coordinates. "Lambda" is an im object resulting
# from density.ppp(), from predict.ppm() or converted from any other rasterized image. 

maple.lambda<- density.ppp(lansing[lansing$marks=="maple"])
ipsim(pp=lansing, mimark="maple", lambda=maple.lambda)

# Build a list of 19 multivariate point pattern where maple is simulated according to an
# Inhomogeneous Poisson process from a predefined  intensity surface "lambda" and where all
# the other species are keep fixed in their original coordinates. "Lambda" is an im object
# resulting  from density.ppp(), # from predict.ppm() or converted from any other rasterized
# image. 

  # Estimate the intensity of maple
   maple.lambda<- density(unmark(lansing[lansing$marks=="maple"]))

  # first simulate the individual maple patterns
   maple.sim<- simulador2(mimark="maple", milambda=maple.lambda, nsim=19)

  # Then, mix the simulated maple patterns with the rest of the multivariate pattern
  # (which remains "fixed")
   multi.maple.sim<- ipsimlist(pp=lansing, mimark="maple", listsim=maple.sim)


## Use of multifocalsimulator() ##

# The same but in a single step with multifocalsimulator(): Build a list of 19 multivariate
# point pattern where maple is simulated according to an Inhomogeneous Poisson process
# from a predefined  intensity surface "lambda" and where all the other species are keep 
# fixed in their original coordinates. "Lambda" is an im object resulting from density.ppp(),
# from predict.ppm() or converted from any other rasterized image. 

   # Estimate the intensity of maple
   maple.lambda<- density(unmark(lansing[lansing$marks=="maple"]))

   # get 99 simulated multivariate point patterns where only maple varies,
   # according to an inhomogeneous Poisson process
   multi.maple.sim <- multifocalsimulator(lansing, "maple", nsim=99,
                                         simulate=expression(rpoispp(maple.lambda)))

\dontrun{
   # Use the simulated multivariate patterns to compute  envelopes for the ISAR against
   # a null model of IPP for maple 
   isar.maple.env<- envelope(lansing, fun=isar, mimark="maple", nsim=99, savefuns=TRUE,
                                         r=seq(0.01, 0.25, le=100), simulate=multi.maple.sim)
   plot( isar.maple.env)



# Use multifocalsimulator()  to compute  envelopes for the ISAR against a null model of
# Poisson cluster for maple.

   # First, adjust a Poisson Cluster process to maple
   require(ecespa)
   maple.pc<- ipc.estK(unmark(lansing[lansing$marks=="maple"]))

  # generate list of simulated multivariate partterns (all other species fixed and maple
  # simulated according to the adjusted PC process):
   maple.pc.sim<-multifocalsimulator(pp=lansing, mimark="maple", nsim=99,nmin=NULL,
                                                 simulate=expression(rIPCP(maple.pc)))

  # compute envelopes
   isar.maple.pc.env<- envelope(lansing, fun=isar, mimark="maple", nsim=99, savefuns=TRUE, 
                                            simulate=maple.pc.sim, r=seq(0.01, 0.25, le=100))
   plot( isar.maple.pc.env)


# Compute  envelopes for the IFDAR against a null model of Poisson cluster for sp_44 in
# San Francisco forest.
   data(SF)
   data(SFtraits)

   # first, get the original point pattern but with marks only for the species (i.e., dsicard 
   # the data.frame of marks and keep only the vector of species names)
   SFsp<- unmark(SF)
   marks(SFsp)<- SF$marks$species

   # second, adjust a PCP to sp_44
   sp_44.pc<- ipc.estK(unmark(SFsp[SFsp$marks=="sp_44"]))

  # third, generate multivariate simulated patterns with only sp_44 varying according
  # to the adjusted PCP
   sp_44.pc.sim<-multifocalsimulator(pp=SFsp, mimark="sp_44",nsim=99,nmin=NULL,
                                                 simulate=expression(rIPCP(sp_44.pc)))

   # finally, compute envelopes
   ifdar.sp_44.pc.env<- envelope(SFsp, fun=ifdar, mimark="sp_44", traits=SFtraits, r=1:30,
                        correct.trait.na=TRUE, nsim=99, simulate=sp_44.pc.sim,  savefuns=TRUE)

   plot( ifdar.sp_44.pc.env)

}
}
\keyword{ spatial }
