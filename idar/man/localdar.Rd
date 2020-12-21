\name{localdar}
\alias{localdar}
\alias{fdismap}
\alias{raoDmap}
\title{
Map Local Diversity Area Relationships
}
\description{
Estimates and maps local diversity-area relationships.
}
\usage{
localdar(mippp, mippp.sp = NULL, nx = NULL, ny = NULL, mimark = NULL, idar = "isar", 
            buffer = 0, bfw = NULL, r, cross.idar = FALSE, tree = NULL, traits = NULL, 
            namesmark = NULL, correct.trait.na = TRUE, correct.trait = "mean",
            correct.phylo="mean")
fdismap(comm, traits)
raoDmap(comm, phy = NULL)
}

\arguments{
 \item{mippp}{A multitype (a.k.a. multivariate) marked point pattern. An object with the \link[spatstat]{ppp} format of \pkg{spatstat}.}
 \item{mippp.sp}{Univariate point pattern of the focal species. An object with the \link[spatstat]{ppp} format of \pkg{spatstat}.}
 \item{nx}{Number of points of the grid along the x axis.}
 \item{ny}{Number of points of the grid along the y axis.}
 \item{mimark}{Character. Name of the focal species in the multitype \code{mippp}.}
 \item{idar}{Character. The name of the idar function to be computed. Either "isar", "ipsvar", "ipsrar", "ipsear", "ipscar", "icwmar", "icwmar.O",  "iraodar"or "imntdar"}
 \item{buffer}{One of "adapt", i.e., compute an adaptive buffer, or a number indicating the width of a fixed buffer area around the plot border}
 \item{bfw}{An \link[spatstat]{owin} object indicating the limits of the buffer area.}
 \item{r}{Vector of distances to compute IDAR(r) functions}
 \item{cross.idar}{Logical. If \code{TRUE}, the focal pattern will be excluded form the community being measured. }
  \item{tree}{A phylogenetic tree in \code{phylo} format (\pkg{ape}) or a phylogenetic covariance matrix}
  \item{traits}{A \link[base]{data.frame} of traits, or a distance matrix among species (in \code{\link[stats]{dist}} or \code{\link[base]{matrix}} format) computed on a data.frame of traits.}
  \item{namesmark}{Character. If the marks in \code{mippp} are within a \code{data.frame}, the name of the column with the species names}
  \item{correct.trait.na}{Logical flag indicating whether \code{NA} values in the matrix of traits should be "corrected": NA values will be assigned the mean trait value.}
 \item{correct.trait}{Character. Either \code{"mean"} or \code{"exclude"}. Species missing in the data.frame of traits will be assigned mean trait values or  will be excluded from the analysis, respectively. }
  \item{correct.phylo}{Character. Either \code{"mean"} meaning \emph{"include missing species in the tree with a constant mean phylogenetic covariance"} or \code{"exclude"}, meaning \emph{"exclude missing species in the tree from the analysis"} }
  \item{comm}{A comunity data table (sites x species).}
 \item{phy}{A comunity data table (sites x species).}
  
}


\details{
\code{localdar} estimates any of the individual diversity area indices (isar, ipsvar,ipsrar, ipsear, ipscar, icwmar, icwmar.O,  iraodar or imntdar) at specific locations, such as the locations of trees of a "focal" point partern or at some grid of points (i.e., "mapping" local diversity area relationships). If no predefined focal point  pattern is provided (by the arguments \code{mippp.sp} or \code{mimark}), \code{localdar} will, by default, estimate the selected idar function in a 30 x 30 point-grid and return a map (the size of the grid can be modified by the arguments \code{nx} and \code{ny}).
\code{fdismap} and \code{raoDmap} are internal functions used by \code{localdar} to get the individual componentes of iraoD and FDis, instead of the averages obtained by \code{\link{envelope4idar}}.

}
\value{
If a focal point patternn has been provided by the arguments \code{mippp.sp} or \code{mimark}, \code{localdar} will return a list of marked point patterns (as long as the vector \code{r}) with the marks showing the estimation of the selected diversity index for the local community defined by a circle of radius \emph{r} around each of the points of the focal pattern. If no focal point pattern is provided, it will return a list of maps (as long as the vector \code{r})  each with the format \link[spatstat]{im} of \pkg{spatstat}, each pixel showing the estimation of the selected diversity index for the local community defined by a circle of radius \emph{r} around the pixel center.
}

\author{ Marcelino de la Cruz \email{marcelino.delacruz@urjc.es}}

\seealso{ \code{\link{envelope4idar}}}

\examples{

# Map of local species area-relationship [ISAR(r)] in San Francisco plot at r=5 and r= 6 m
  data(SF)
  data(SFphylotree)
  data(SFtraits)

  isarSF <- localdar(SF, r=5:6, namesmark="species")
  
  isarSF
  
  plot(isarSF[[1]])
  
  
# Map of local species area-relationship [ISAR(r)] in San Francisco plot at r=5 and r= 6 m,
# with an adaptive buffer

isarSF <- localdar(SF, r=5:6, namesmark="species",buffer="adapt")



# Map of local species area-relationship [ISAR(r)] in lansing woods at different radii,
# with an fixed  buffer (only for rectangular windows).

data(lansing)
lansing.bfw<- owin(c(0.2,0.8),c(0.2,0.8))
lansing.lsar.bf<-localdar(lansing, r=seq(0.05,0.2,by=0.05), bfw=lansing.bfw)


# Map of local species area-relationship [ISAR(r)] in San Francisco plot at r=5 and r= 6 m,
# with a buffer of 7 m within plot limits

  mibfw<- erosion(SF$win, r=7)

 isarSF <- localdar(SF, r=5:6, namesmark="species", bfw=mibfw)



# Estimate local species area-relationship [ISAR(r)] in the local communities 
#  in circles with radius  r=5 and r= 6 m around the individuals of sp_44

  sp_44_sar <- localdar(SF, r=5:6, namesmark="species", mimark="sp_44")
  sp_44_sar
  plot(sp_44_sar[[1]])
  sp_44_sar[[1]]$marks
  


# Estimate local species area-relationship [ISAR(r)] in the local communities 
#  in circles with radius  r=5 and r= 6 m around the individuals of sp_44
# EXCLUDING the focal species from species counts.

  sp_44_sar <- localdar(SF, r=5:6, namesmark="species", mimark="sp_44", cross.idar=TRUE)
  sp_44_sar
  plot(sp_44_sar[[1]])
  sp_44_sar[[1]]$marks
 
  
  


# Map and estimation of local Phylogenetic Species Variety
local_psvar<-localdar(SF, r=5:6, idar="ipsvar", tree=SFphylotree, namesmark="species")
sp44_psvar <- localdar(SF, r=5:6, idar="ipsvar", tree=SFphylotree, namesmark="species",
                                    mimark="sp_44", buffer="adapt")

# Map and estimation of local Functional Dispersion
local_fdar <- localdar(SF, nx=50, ny=25, r=5:6,  idar="ifdar", traits=SFtraits,
                                 namesmark="species", correct.trait.na=TRUE)
sp44_fdar <- localdar(SF, nx=50, ny=25, r=5:6,  idar="ifdar", traits=SFtraits,
                                  namesmark="species", mimark="sp_44", correct.trait.na=TRUE)


# Map of a local community weighted mean of wood density
# first, put the wood desnsity data as a named vector
 wood.density.vec<-unlist(SFtraits[,  "wood.density", drop=FALSE])
    names(wood.density.vec)<- rownames(SFtraits)
    
local_cwd <- localdar(SF, r=5:6,  idar="icwmar", traits=wood.density.vec,
                                  namesmark="species", correct.trait.na=TRUE)
local_O.cwm <-localdar(SF, r=5:6, idar="icwmar.O", traits=wood.density.vec, 
                                     namesmark="species", correct.trait.na=TRUE)

# Map of Rao's phylogenetic diversity
local_rao<- localdar(SF,  r=5:6, idar="iraodar",  tree=SFphylotree, namesmark="species")
local_O.rao <- localdar(SF, r=5:6, idar="iraodar.O",  tree=SFphylotree, namesmark="species")

# Map of local mean nearest taxon distance 
local_mntd <- localdar(SF, r=5:6, idar="imntdar",  tree=SFphylotree, namesmark="species")


}

\keyword{ spatial }

