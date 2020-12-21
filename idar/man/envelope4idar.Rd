\encoding{Latin1}
\name{envelope4idar}
\alias{envelope4idar}
\alias{idar2}
\alias{raoDmod}

\title{Compute Simulation Envelopes for IDAR(r) Functions}
\description{
Compute simulation envelopes for IDAR(r) functions.
}

\usage{
  envelope4idar(mippp, mippp.sp.sim, mippp.sp, mimark=NULL,
    namesmark=NULL, r, idar="isar", buffer="adapt", bfw=NULL,
    nsim=NULL, nrank=1, tree = NULL, traits = NULL, 
    cross.idar=FALSE, savefuns=TRUE, correct.phylo="exclude",
    correct.trait.na=FALSE, correct.trait="mean" )
idar2(mippp.sp, mippp, mimark, idar="isar", buffer,bfw, r, 
    cross.idar=FALSE,  tree = NULL, traits = NULL)
    raoDmod(comm, phy = NULL) 
  
}

\arguments{
  \item{mippp}{A multitype (a.k.a. multivariate) marked point pattern. An object with the \link[spatstat]{ppp} format of \pkg{spatstat}.}
  \item{mippp.sp}{Univariate point pattern of the focal species. An object with the \link[spatstat]{ppp} format of \pkg{spatstat}.}
  \item{mippp.sp.sim}{A list with simulations of the focal species point pattern created with \code{\link{simulador2}} }
  \item{mimark}{Character. Name of the focal species in the multitype \code{mippp}.}
  \item{namesmark}{Character. If the marks in \code{mippp} are within a \code{data.frame}, the name of the column with the species names}
  \item{r}{Vector of distances to compute IDAR(r) functions}  
  \item{idar}{Character. The name of the idar fucntion to be computed. Either "isar", "ipsvar", "ipsrar", "ipsear", "ipscar", "icwmar", "icwmar.O",  "iraodar"or "imntdar"}
  \item{buffer}{Character or numeric. Either "adapt" (i.e., compute an adaptive buffer), or a number indicating de width of a fixed buffer area around the plot border}
  \item{bfw}{An \link[spatstat]{owin} object indicating the limits of the buffer area.}
  \item{nsim}{The number of simulations.}
  \item{nrank}{Integer. Rank of the envelope value amongst the nsim simulated values. A rank of 1 means that the minimum and maximum simulated values will be used.}
  \item{tree}{A phylogenetic tree in \code{phylo} format (\pkg{ape}) or a phylogenetic covariance matrix}
  \item{traits}{A \link[base]{data.frame} of traits, or a distance matrix among species (in \code{\link[stats]{dist}} or \code{\link[base]{matrix}} format) computed on a data.frame of traits.}
  \item{cross.idar}{Logical. If \code{TRUE}, the focal pattern be excluded form the community being measured. }
   \item{savefuns}{Logical flag indicating whether to save all the simulated function values.}
   \item{correct.phylo}{Character. Either \code{"mean"} meaning \emph{"include missing species in the tree with a constant mean phylogenetic covariance"} or \code{"exclude"}, meaning \emph{"exclude missing species in the tree from the analysis"} }
   \item{correct.trait.na}{Logical flag indicating whether \code{NA} values in the matrix of traits should be "corrected": NA values will be assigned the mean trait value.}
   \item{correct.trait}{Character. Either \code{"mean"} or \code{"exclude"}. Species missing in the data.frame of traits will be assigned mean trait values or  will be excluded from the analysis, respectively. }
   \item{comm}{A data.frame with community data}
\item{phy}{a phylogenetic tree in \code{phylo} format (\pkg{ape}) or more probably, a distance matrix.}
}

\details{
In 2007, Wiegand et al. developed the concept of Individual Species-Area Relationship. Basically, this consist in computing species accumulation curves by samping areas with varying radius \eqn{r} around the individual trees of a \emph{focal} species. Here we provide a version of \emph{ISAR(r)} (\code{idar="isar"}), but we also extend this concept to other diversity-area relationships and provide functions to compute individual phylogenetic diversity-area and individual functional diversity-area relationships. The individual phylogenetic functions  are based in Helmus et al. (2007) measures, i.e., phylogenetic species variability (\code{idar="ipsvar"}), phylogenetic species richness (\code{idar="ipsrar"}), phylogenetic species evenness (\code{idar="ipsear"}), and phylogenetic species clustering (\code{idar="ipscar"}). Also, an individual version of the mean nearest taxon distance of Webb et al. (2002) (\code{idar="imntd"}).  The individual functional-diversity function (\code{idar="ifdar"}) is based in the functional dispersion measure (FDis) of \enc{Laliberté and Legendre (2010)}{Laliberte and Legendre (2010)}. Other available functions are based in a weighted community mean of traits (\code{idar="icwmar"}),  on  a verison of Rao quadratic entropy (\code{idar="iraodar"}), and on the same functions but computed in "rings" defined by the neighbour \code{r} values (\code{idar="icwmar.O"} and \code{idar="iraodar.O"}), which would eliminate the "accumulative" or "memory" effect on \emph{ICWMAR(r)} and \emph{IRAODAR(r)}. 

Although recent litterature (e.g., Wiegand and Moloney 2014) suggest that buffer correction is not necessary for this type of statistics, and by default all functions are estimated without buffer (e.g., \code{buffer=0}), several  edge correction coould be employed. For example,  an adaptative buffer correction could be used (\code{buffer="adapt"}), i.e., for each radius \eqn{r}, only individuals of the focal species that are placed at a distance \eqn{>=r} from the border of the plot are considered in the computation of the different meassures. It is also possible to set a fixed buffer width  (e.g., \code{buffer=30}), which will accelerate te computations but will discard many useful trees. It is also possible to provide also a fixed window (in the argument \code{bfw}) to indicate the limits of the buffer area. This could be useful to computing the IDAR(r) functions in different subsets of the original plot (e.g., in different "habitats").

In general, computing envelopes with \code{envelope4idar} is a little faster than using \link[spatstat]{envelope} and the individual functions (e.g., \link{ipsvar}, \link{ifdar}, etc). In addition, \code{envelope4idar} has the possibility of computing "crossed" individual functions, i.e., using a focal species that is not part of the community whose diversities are being measured. This allows evaluating the diversity of e.g., young trees around older trees, etc. This is accomplished by seting the argument \code{cross.idar=TRUE}.

While  \code{envelope4idar} manages data and results, \code{idar2} actually computes the individual functions (both observed and simulated) . In general, \code{idar2} would not be called directly by the user.

\code{raoDmod} is a modification of the function \link[picante]{raoD} of \pkg{picante} to accept distance matrices instead of phylogenetoc trees. It would not be called directly by the user.

}

\section{Warning}{
       The transcription of species names in the multivariate \code{mippp}, in the row names of  the data.frame of traits (or in the names or dimnames of  the distance matrices) should be identical. The same applies to the tiplabels of the phylogenetic tree.
     }


\value{
  An object of class "\code{fv}", see  \link[spatstat]{fv.object}, and \link[spatstat]{envelope},which can be plotted directly using  \link[spatstat]{plot.fv}. 
Essentially a data frame containing columns 
   \describe{
      \item{r}{the vector of values of the argument r at which the \code{idar(r)} function has been estimated}
 
      \item{obs}{values of the summary function for the data point pattern }
 
      \item{lo}{lower envelope of simulations}
 
      \item{hi}{upper envelope of simulations }

       \item{nmean}{estimated theoretical value of the summary function under the considered null model, computed by averaging simulated values}

   }


}

\seealso{
       \code{\link[picante]{psd}} for a description of the phylogenetic measures of Helmus et al. (2007).

        \code{\link[FD]{fdisp}} for a description of the functional dispersion measure (FDis) of \enc{Laliberté and Legendre (2010)}{Laliberte and Legendre (2010)}.
}


\references{

   Helmus M.R., Bland T.J., Williams C.K. and Ives A.R. (2007) Phylogenetic measures of biodiversity. \emph{American Naturalist}, 169, E68-E83.

   \enc{Laliberté, E. and  Legendre, P. (2010) A distance-based framework for measuring functional diversity from multiple traits. \emph{Ecology} 91, 299-305.}{Laliberte, E. and  Legendre, P. (2010) A distance-based framework for measuring functional diversity from multiple traits. \emph{Ecology} 91, 299-305.} 

   Wiegand,T., Gunatilleke, C.V.S., Gunatilleke, I.A.U.N. and  Huth, A. (2007) How individual species structure diversity in tropical forests. \emph{PNAS} 104, 19029-19033.

   Webb, C., Ackerly, D., McPeek, M., and Donoghue M.  2002. Phylogenies and community ecology. \emph{Annual Review of Ecology and Systematics} 33:475-505

}

\author{ Marcelino de la Cruz \email{marcelino.delacruz@urjc.es}}


\examples{



# compute envelope for isar around sp_44
data(SF)
data(SFtraits)
data(SFphylotree)

  
# 1) Create a list with simulations of the focal species point pattern  with simulador2()
#     or by hand. # Beware that each simulated ppp should be marked (with the mark
#      of the focal species)

#  Example for simulations of an inomogenous PP
# Adjust an IPP:
    
    sp_44.ipp<- density.ppp(unmark(SF[SF$marks$species=="sp_44"]))
   
# simulate 19 realizations of the adjusted IPP
# (BEWARE: in real tests you should use 199 or higher)
    sp_44.ipp.sim <- vector(mode="list", le=19)
    sp_44.ipp.sim <- lapply(sp_44.ipp.sim, function(x) x=rpoispp(sp_44.ipp))
   # mark each simulated pattern
    sp_44.ipp.sim.m <- lapply(sp_44.ipp.sim, function(x) 
                                       {marks(x) = factor(rep("sp_44", x$n)); return(x)})

# 2) compute.envelopes
     # ISAR
  isar.sp_44.ipp.env<- envelope4idar(mippp=SF, mippp.sp.sim= sp_44.ipp.sim.m,
                mimark="sp_44", namesmark="species",r=1:30, buffer=0) 
\dontrun{
     # IFDAR
 ifdar.sp_44.ipp.env<- envelope4idar(mippp=SF, mippp.sp.sim= sp_44.ipp.sim.m,
                mimark="sp_44",  namesmark="species", r=1:30, idar="ifdar", buffer=0, 
		 traits=SFtraits, correct.trait.na=TRUE)

    #IPSVAR
   ipsvar.sp_44.ipp.env<- envelope4idar(mippp=SF, mippp.sp.sim= sp_44.ipp.sim.m,
                 mimark="sp_44", namesmark="species", r=1:30, idar="ipsvar", buffer=0, 
		  tree=SFphylotree)

####################################
# Computing CROSS_IDAR
#--------------------------------------------------

# You need a focal pattern that its not part of the multivariate pattern which is "measured"
# to estimate diversity.
# For example, let's measure diversity of small trees  around large trees of the focal species.

# First, obtain the pattern of large and small trees
# the pattern of small trees will be the "measured" one, i.e., the argument "mippp"

SFlarge<- SF[SF$marks$dbh>=10]
SFsmall <- SF[SF$marks$dbh<10]

# pattern of the focal species (this will be the argument "mippp.sp")
sp_44.large<- SFlarge[SFlarge$marks$species=="sp_44"]

# list of simulated patterns of the focal species (e.g., from an IPP)
sp_44.large.ipp<- density.ppp(sp_44.large)
sp_44.large.ipp.sim <- vector(mode="list", le=99)
sp_44.large.ipp.sim <- lapply(sp_44.large.ipp.sim, function(x) x=rpoispp(sp_44.large.ipp))

# COMPUTE envelopes for cross-ISAR(r)
 isar.sp_44.large.cross.ipp.env<- envelope4idar(mippp=SFsmall, r=1:30, buffer=0,
			mippp.sp.sim= sp_44.large.ipp.sim, mippp.sp=sp_44.large,
                        namesmark="species", cross.idar =TRUE) 

# COMPUTE envelopes for cross-IFDAR(r)
ifdar.sp_44.large.cross.ipp.env<- envelope4idar(mippp=SFsmall, r=1:30, idar="ifdar", 
                       buffer=0, mippp.sp.sim= sp_44.large.ipp.sim,  mippp.sp=sp_44.large,
                       namesmark="species", traits=SFtraits, correct.trait.na=TRUE,
                       cross.idar =TRUE) 

# COMPUTE envelopes for cross-IPSVAR(r)
ipsvar.sp_44.large.cross.ipp.env<- envelope4idar(mippp=SFsmall, r=1:30, idar="ipsvar", 
                      buffer=0, mippp.sp.sim= sp_44.large.ipp.sim,  mippp.sp=sp_44.large,
                      namesmark="species",  tree=SFphylotree, cross.idar =TRUE) 


####################################
# Comparing the performance of envelope() and envelope4idar()
# 
#--------------------------------------------------
require(ecespa) # for the ipc.estK() function
data(SF)
SFsp<- unmark(SF)
marks(SFsp)<- SF$marks$species
sp_44.ppp<-unmark(SF[SF$marks$species=="sp_44"])
sp_44.pc<- ipc.estK(sp_44.ppp)


# use multifocalsimulator() to use the same simulations with both functions.
sp_44.pc.sim.mf0<-multifocalsimulator(pp=SFsp, mimark="sp_44",
                           simulate=expression(rIPCP(sp_44.pc)), nsim=99, nmin=sp_44.ppp$n)

# envelopes with function envelope()
gc()
t0<- Sys.time()
  ifdar.sp_44.pc.env.e2<- envelope(SFsp, fun=ifdar, mimark="sp_44", traits=SFtraits,
                            correct.trait.na=TRUE, nsim=99, simulate=sp_44.pc.sim.mf0, r=1:30,
                            savefuns=TRUE, buffer=0)
Sys.time()-t0

# envelopes with function envelope4idar()

# Here you should input simulated patterns only for the focal species so, first,
# extract it fom the list of simulated multivariate ppp

sp_44.pc.sim.mf00<- lapply(sp_44.pc.sim.mf0, function(x) {x=x[x$marks=="sp_44"];return(x)})

gc()
t0<- Sys.time()
ifdar.sp_44.pc.env2<- envelope4idar(mippp=SF, mippp.sp.sim= sp_44.pc.sim.mf00, mimark="sp_44",
				namesmark="species", r=1:30, idar="ifdar", buffer=0, 
                                nsim=99, traits=SFtraits, correct.trait.na=TRUE)
Sys.time()-t0


}
}

\keyword{ spatial }
