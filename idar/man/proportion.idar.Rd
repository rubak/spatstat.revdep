\encoding{UTF-8}
\name{proportion.idar}
\alias{proportion.idar}
\alias{plot.pidar}


\title{Determine the Proportion of Accumulator and Repeller Species}
\description{
  This function determines the proportion of species in a community that "accumulate" or "repell" diversity, in the sense of Wiegand et al. (2007).
}

\usage{
  proportion.idar(envlist, alfa = 0.05)
     \method{plot}{pidar}(x,cols=c(1,2,3),type=c("l","o","o"), pch=c(NA,19,19),
    lty=c(1,1,1),  legend=TRUE, p.legend="topleft",...)
}

\arguments{
  \item{envlist}{A list whose elements are the result of applying \link[spatstat.core]{envelope} to a set of point patterns.}
  \item{alfa}{alpha value to calculate deviations from null model}
  \item{x}{The result of \code{proportion.idar}}
  \item{cols}{A vector (length =3) with the color names or codes for the curves of each of the 3 categories: "neutral", "accumulator", "repeller".}
  \item{type}{A vector (length =3), with the type of plot for each category. See \link[graphics]{plot.default}.}
  \item{pch}{A vector (length =3), with the point type  for each category. See \link[graphics]{points}}
  \item{lty}{A vector (length =3), with the line type  for each category. See \link[graphics]{par}}
   \item{legend}{Logical. Should a legend be added to the plot?}
  \item{p.legend}{Position of the legend. It can be a vector wit the xycoordinates of the upperleft corner of the legend box or a  keyword accepted by \code{\link[graphics]{legend}}}
  \item{...}{Additional graphical parameters passed both to funtions \code{plot} and \code{lines}}
}

\details{
This function determines the proportion of accumulator, repeller and neutral species at each scale \emph{r}, following the approach of Wiegand et al (2007). A species is classified as an accumulator at scale \emph{r} if there are less than \eqn{(nsim+1) *alpha/2} simulated values greater than the observed \emph{idar(r)}. On the contrary, a species is classified as repeller at scale if there are less than \eqn{(nsim+1) *alpha/2} simulated values smaller than the observed \emph{idar(r)}.  The percentage is computed over the total of species in\code{ envlist}. It is necessary that the objects in envlist (i.e. the envelope objects) had been computed with the argument "\code{simfuns=TRUE}".

}
\value{
\code{proportion.idar} produces an object of class "\code{pidar}", basically a list with components

   \describe{
        \item{percentage}{A \code{data.frame} with 3 columns ("p.accumulators", "p.repellers", "p.neutrals") indicating the percentage in each category for each \emph{r} (rows)}
        \item{nsp}{The total number of species for which the percentage has been computed.}
        \item{nsim}{The number of simulations}
        \item{alfa}{The alpha value employed to calculate the deviations from null model and the assignation to each of the categories}
        \item{r}{The vector of r values at which the idar(r) functions have been estimated}
         \item{behaviour}{A data.frame with the behaviour ("A"= accumulator, "R" = repeller) of each species at each radius r}
    }
   
   \code{plot.pidar} plots the result.

}
\references{
   Wiegand,T., Gunatilleke, C.V.S., Gunatilleke, I.A.U.N. and  Huth, A. (2007) How individual species structure diversity in tropical forests. \emph{PNAS} 104, 19029-19033.

}

\author{ Marcelino de la Cruz \email{marcelino.delacruz@urjc.es}}

\examples{

# Compute percentage of accumulator, neutral an repeller species (for ISAR) in Lansing woods
# In this example, against a null model of IPP.
data(lansing)

# Compute rISAR and plot envelopes for an inhomogeneous Poisson model
# of each species in San Francisco plot
# Split sfsp point pattern ppp by species
lansing.sp<- split(lansing)

# Set parameters for the simulations, estimation of intensity surface for IPP, etc.
# BEWARE: THIS is for R-TESTING ALONE. REAL TESTS SHOULD EMPLOY 199 simulations or higher
nsim<-8
r<- seq(0.01,0.25, by=0.01)
# Create list to store the results
isar.lansing<- list()
# start computation

for( i in 1: length(lansing.sp)){
   print(i)
   # estimate intensity of the focal species
   lambda<- density(lansing.sp[[i]])
   # obtain simulated patterns were all species except the focal remain fixed
   # and the focal varies according to an inhomomgeneous Poiisson process
   simulados<- multifocalsimulator(lansing, mimark=levels(lansing$marks)[i],
                            simulate=expression(rpoispp(lambda)), nsim=nsim)
  # compute isar and envelopes for each species
  # It is COMPULSORY that the argument "savefuns" it is set to TRUE
   isar.lansing[[i]] <- envelope(lansing, isar, mimark=levels(lansing$marks)[i],  r=r,
         simulate=simulados,nsim=nsim, savefuns=TRUE, buffer=0)
}


prop.isar.lansing <- proportion.idar(isar.lansing)
head(prop.isar.lansing $percentage)
head(prop.isar.lansing $behaviour)
plot(prop.isar.lansing, p.legend=c(0.15,60))

\dontrun{

# Compute percentage of accumulator, neutral an repeller species (for ISAR) in San Francisco
# forest.  In this example, against a null model of IPP.
data(SF)

# Discard the size mark and keep the species mark in SF ppp:
sfsp<- ppp(SF$x, SF$y, window=SF$window, marks=SF$marks$species)

# Compute ISAR and plot envelopes for an inhomogeneous Poisson model
# of each species in San Francisco plot
# Split sfsp point pattern ppp by species
sfsp.sp<- split(sfsp)
# Select species with >= 10 individuals (to get some statistical power)
sfsp10 <- sapply(sfsp.sp, function(x) x$n>=10)
#Get names of those species
nombressf<- names(sfsp10[sfsp10])
# Set parameters for the simulations, estimation of intensity surface for IPP, etc.
nsim<-199
nmin<-10
sigma <- 8
r<- seq(1,15, by=0.5)
# Create list to store the results
isar.sf<- list()
# start computation
# BEWARE: THIS TAKES QUITE A FEW MINUTES!!!
for( sp in nombressf){
   print(sp)
   # estimate intensity of the focal species
   lambda<- density(unmark(sfsp[sfsp$marks==sp]), sigma=sigma, positive=TRUE)
   # obtain simulated patterns were all species except the focal remain fixed
   # and the focal varies according to an inhomomgeneous Poiisson process
   simulados<- multifocalsimulator(sfsp, mimark=sp,
                            simulate=expression(rpoispp(lambda)), nsim=nsim,nmin=nmin)
  # compute isar and envelopes for each species
  # It is COMPULSORY that the argument "savefuns" it is set to TRUE
   isar.sf[[sp]] <- envelope(sfsp, isar, mimark=sp,  r=r,
         simulate=simulados,nsim=nsim, savefuns=TRUE, buffer=0)
}


prop.isar.SF <- proportion.idar(isar.sf)
head(prop.isar.SF$percentage)
head(prop.isar.SF$behaviour)
plot(prop.isar.SF)

}
}

\keyword{ spatial }
