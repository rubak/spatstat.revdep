\encoding{Latin1}
\name{isar}
\alias{isar}
\alias{ipscar}
\alias{ipsear}
\alias{ipsvar}
\alias{ipsrar}
\alias{ifdar}
\title{Individual Diversity Area Relationships}
\description{
  Estimate different Individual Diversity-Area Relationships from a multivariate point pattern.
}

\usage{
  isar(mippp, mippp.sp=NULL, mimark=NULL,  namesmark=NULL, r=NULL,
   buffer=0, bfw=NULL)
  ipscar(mippp, mippp.sp=NULL, mimark=NULL,  namesmark=NULL,
   tree=NULL, r=NULL, buffer=0, bfw=NULL, correct.phylo="mean")
  ipsear(mippp, mippp.sp=NULL, mimark=NULL,  namesmark=NULL,
   tree=NULL, r=NULL, buffer=0, bfw=NULL, correct.phylo="mean")
  ipsvar(mippp, mippp.sp=NULL, mimark=NULL,  namesmark=NULL,
   tree=NULL, r=NULL, buffer=0, bfw=NULL, correct.phylo="mean")
  ipsrar(mippp, mippp.sp=NULL, mimark=NULL,  namesmark=NULL,
   tree=NULL, r=NULL, buffer=0, bfw=NULL, correct.phylo="mean")
  ifdar(mippp, mippp.sp=NULL, mimark=NULL,  namesmark=NULL,
   traits=NULL, r=NULL, buffer=0, bfw=NULL, correct.trait.na=FALSE,
     correct.trait="mean")
}

\arguments{
  \item{mippp}{A multitype (a.k.a. multivariate) marked point pattern. An object with the \link[spatstat.geom]{ppp} format of \pkg{spatstat}.}
  \item{mippp.sp}{Univariate point pattern of the focal species. An object with the \link[spatstat.geom]{ppp} format of \pkg{spatstat}.}
  \item{mimark}{Character. Name of the focal species in the multitype \code{mippp}.}
  \item{namesmark}{Character. If the marks in \code{mippp} are within a \code{data.frame}, the name of the column with the species names}
  \item{buffer}{One of "adapt", i.e., compute an adaptive buffer, or a number indicating de width of a fixed buffer area around the plot border}
  \item{bfw}{An \link[spatstat.geom]{owin} object indicating the limits of the buffer area.}
  \item{r}{Vector of distances to compute IDAR(r) functions}
  \item{tree}{A phylogenetic tree in \code{phylo} format (\pkg{ape}) or a phylogenetic covariance matrix}
  \item{traits}{A \link[base]{data.frame} of traits, or a distance matrix among species (in \code{\link[stats]{dist}} or \code{\link[base]{matrix}} format) computed on a data.frame of traits.}
  \item{correct.phylo}{Character. Either \code{"mean"} meaning \emph{"include missing species in the tree with a constant mean phylogenetic covariance"} or \code{"exclude"}, meaning \emph{"exclude missing species in the tree from the analysis"} }
   \item{correct.trait.na}{Logical flag indicating whether \code{NA} values in the matrix of traits should be "corrected": NA values will be assigned the mean trait value.}
   \item{correct.trait}{Character. Either \code{"mean"} or \code{"exclude"}. Species missing in the data.frame of traits will be assigned mean trait values or  will be excluded from the analysis, respectively. }

}

\details{
In 2007, Wiegand et al. developed the concept of Individual Species-Area Relationship. Basically, this consist in computing species accumulation curves by samping areas with varying radius \eqn{r} around the individual trees of a \emph{focal} species. Here we extend this concept to other diversity-area relationships and provide functions to compute individual phylogenetic diversity-area and individual functional diversity-area relationships. The individual phylogenetic functions are based in Helmus et al. (2007) measures, i.e., phylogenetic species variability (\code{ipsvar}), phylogenetic species richness (\code{ipsrar}), phylogenetic species evenness (\code{ipsear}), and phylogenetic species clustering (\code{ipscar}). The individual functional-diversity function (\code{ifdar}) is based in the functional dispersion measure (FDis) of \enc{Laliberté and Legendre (2010)}{Laliberte and Legendre (2010)}.

Although recent litterature (e.g., Wiegand and Moloney 2014) suggest that buffer correction is not necessary for this type of statistics, and by default all functionsare estimated without buffer (e.g., \code{buffer=0}), several  edge correction coould be employed. For example,  an adaptative buffer correction could be used (\code{buffer="adapt"}), i.e., for each radius \eqn{r}, only individuals of the focal species that are placed at a distance \eqn{>=r} from the border of the plot are considered in the computation of the different meassures. It is also possible to set a fixed buffer width  (e.g., \code{buffer=30}), which will accelerate te computations but will discard many useful trees. It is also possible to provide also a fixed window (in the argument \code{bfw}) to indicate the limits of the buffer area. This could be useful to computing the IDAR(r) functions in different subsets of the original plot (e.g., in different "habitats"). 


}

\section{Warning}{
       The transcription of species names in the multivariate \code{mippp}, in the row names of  the \code{data.frame} of traits (or in the names or \code{dimnames} of  the distance matrices) should be identical. The same applies to the \code{tiplabels} of the phylogenetic tree.
}


\value{
  An object of class "\code{fv}", see  \link[spatstat.core]{fv.object}, which can be plotted directly using  \link[spatstat.core]{plot.fv}. 

Essentially a data frame containing a column named \emph{r} with the vector of values of the argument r at which the IDAR(r) function has been estimated and  aonther column, named \emph{"isar", "ipsvar", "ipsrar", "ipsear", "ipscar"} or \emph{"ifdar"}, according to the selected idar argment. This column contains an estimate of the selected IDAR(R) function.

}

\section{Simulation envelopes}{
To compute simulation envelopes for the IDAR(r) functions,  use \code{\link[spatstat.core]{envelope}}. See the examples in this help page and in \code{\link{ipsim}} to know how to compute simulation envelopes from appropriate null models.

To compute envelopes for "crossed" IDAR(r) functions or to accelerate the computation of "single" IDAR(r) functions, use \code{\link{envelope4idar}}.
}


\seealso{
       \code{\link[picante]{psd}} for a description of the phylogenetic measures of Helmus et al. (2007).

        \code{\link[FD]{fdisp}} for a description of the functional dispersion measure (FDis) of \enc{Laliberté and Legendre (2010)}{Laliberte and Legendre (2010)}.
}


\references{

   Helmus M.R., Bland T.J., Williams C.K. and Ives A.R. (2007) Phylogenetic measures of biodiversity. \emph{American Naturalist}, 169, E68-E83.

   \enc{Laliberté, E. and  Legendre, P. (2010) A distance-based framework for measuring functional diversity from multiple traits. \emph{Ecology} 91, 299-305.}{Laliberte, E. and  Legendre, P. (2010) A distance-based framework for measuring functional diversity from multiple traits. \emph{Ecology} 91, 299-305.} 

   Wiegand,T., Gunatilleke, C.V.S., Gunatilleke, I.A.U.N. and  Huth, A. (2007) How individual species structure diversity in tropical forests. \emph{PNAS} 104, 19029-19033.

}

\author{ Marcelino de la Cruz \email{marcelino.delacruz@urjc.es}}

\examples{

# ISAR
   # Point pattern with a data.frame of marks
   data(SF)
   isar.sp_44 <- isar(mippp = SF, mimark="sp_44", namesmark="species", r=1:40) 
   plot(isar.sp_44)
 
   # Point pattern with just a vector of marks
    data(lansing)
    isar.blackoak <- isar(mippp = lansing, mimark="blackoak",  r=seq(0.01, 0.25, le=100))
    plot(isar.blackoak)

   # Examples of the use of different buffers
     # No  buffer at all (by deffault, buffer = 0)
     isar.sp_44.0 <- isar(mippp = SF, mimark="sp_44", namesmark="species", r=1:18)    

     # Adaptive buffer (for each r, use only points within a r distance form the border)
     isar.sp_44.a <- isar(mippp = SF, mimark="sp_44", namesmark="species", r=1:18,
                                buffer="adapt") 

     # Predefined window, for example with a buffer of 7 m within plot limits
     mibfw<- erosion(SF$win, r=7)
      isar.sp_44.w <- isar(mippp = SF, mimark="sp_44", namesmark="species", r=1:18, bfw=mibfw) 


   ######################
   ### Phylogenetic functions ###
   ######################
   
data(SFphylotree)

# IPSCAR
    ipscar.sp_44 <- ipscar(mippp = SF, mimark="sp_44", namesmark="species", r=1:40,
                                     tree=SFphylotree)
    plot(ipscar.sp_44)
 
# IPSEAR
    ipsear.sp_44 <- ipsear(mippp = SF, mimark="sp_44", namesmark="species", r=1:40,
                                  tree=SFphylotree)
    plot(ipsear.sp_44)

# IPSVAR
    ipsvar.sp_44 <- ipsvar(mippp = SF, mimark="sp_44", namesmark="species", r=1:40,
                                   tree=SFphylotree)
    plot(ipsvar.sp_44)

# IPSRAR
    ipsrar.sp_44 <- ipsrar(mippp = SF, mimark="sp_44", namesmark="species", r=1:40,
                                  tree=SFphylotree)
    plot(ipsrar.sp_44)
 
  #####################
  ###  Functional functions ###
  #####################
 
data(SFtraits)

# IFDAR
   # this will cause an error becuse some species have NA's in the vector of trait values
\dontrun{
  # ifdar.sp_44 <- ifdar(mippp = SF, mimark="sp_44", namesmark="species", traits=SFtraits,
  #                           r=1:40, correct.trait="exclude")
}
   # "correct" NA's in trait values  by assigning tospecies without traits the average of the trait 
   # for  all the other species
   ifdar.sp_44 <- ifdar(mippp = SF, mimark="sp_44", namesmark="species",  traits=SFtraits,
			      r=1:40, correct.trait.na=TRUE)

    #"correct"  the existence of NA's in trait values  by excluding species without traits from the
    # analysis
    ifdar.sp_44 <- ifdar(mippp = SF, mimark="sp_44", namesmark="species", traits=SFtraits,
                               r=1:40, correct.trait.na=TRUE,correct.trait="exclude")

     plot(ifdar.sp_44)

# For examples of envelopes for these functions  see the help page of ipsim() or envelope4idar()

}


\keyword{ spatial }
