\encoding{latin1}
\name{ecespa}
\alias{ecespa-package}
\alias{ecespa}
\docType{package}
\title{
Functions for spatial point pattern analysis in ecology

}
\description{
This is a summary of the features of \pkg{ecespa} by way of examples of its main functions.
}

\author{
Marcelino de la Cruz Rot, with contributions of Philip M. Dixon and Jose M. Blanco-Moreno and heavily borrowing  Baddeley's & Turner's \pkg{spatstat} code.

}

\references{
De la Cruz, M. 2006. \enc{Introducción al análisis de datos mapeados o algunas de las (muchas) cosas  que puedo hacer si tengo coordenadas}{Introduccion al analisis  de datos mapeados o algunas de las (muchas) cosas  que puedo hacer si tengo coordenadas}. \emph{Ecosistemas}  15 (3): 19-39. 

De la Cruz, M. 2008. \enc{Métodos para analizar datos puntuales. 
  En: \emph{Introducción al Análisis Espacial de Datos en Ecología y Ciencias Ambientales: Métodos y Aplicaciones } 
  (eds. Maestre, F. T., Escudero, A. y Bonet, A.), pp 76-127. Asociación Española de Ecología Terrestre, Universidad Rey Juan Carlos
  y Caja de Ahorros del Mediterráneo, Madrid.}{Metodos para analizar datos puntuales. 
  En: \emph{Introduccion al Analisis Espacial de Datos en Ecologia y Ciencias Ambientales: Metodos y Aplicaciones } 
  (eds. Maestre, F. T., Escudero, A. y Bonet, A.), pp 76-127. Asociacion Espanola de Ecologia Terrestre, Universidad Rey Juan Carlos
  y Caja de Ahorros del Mediterraneo, Madrid.}
  
De la Cruz, M., Romao, R.L.,  Escudero, A. and Maestre, F.T. 2008. Where do seedlings go? A spatio-temporal analysis of
 early mortality in a semiarid specialist. \emph{Ecography}, 31(6): 720-730. \url{http://dx.doi.org/10.1111/j.0906-7590.2008.05299.x}.

Diggle, P. J. 2003. \emph{Statistical analysis of spatial point patterns}. Arnold, London. 

Dixon, P.M. 2002. Nearest-neighbor contingency table analysis of spatial segregation for several species.
  \emph{Ecoscience}, 9 (2): 142-151. \url{https://doi.org/10.1080/11956860.2002.11682700}. 
  
Dixon, P. M. 2002. Ripley's K function. In \emph{The encyclopedia of environmetrics} 
(eds. El-Shaarawi, A.H. & Piergorsch, W.W.), pp. 1976-1803. John Wiley & Sons Ltd, NY.
  
Escudero, A., Romao, R.L., De la Cruz, M. & Maestre, F. 2005. Spatial pattern and neighbour effects on 
\emph{Helianthemum squamatum} seedlings in a Mediterranean gypsum community. \emph{Journal of  Vegetation  Science},  16: 383-390. \url{https://doi.org/10.1111/j.1654-1103.2005.tb02377.x}.

Getis, A. and Franklin, J. 1987. Second-order neighbourhood analysis of mapped point patterns. \emph{Ecology} 68: 473-477. \url{https://doi.org/10.2307/1938452}.

Loosmore, N.B. and Ford, E.D. (2006) Statistical inference using the G or K point pattern spatial statistics. \emph{Ecology} 87, 1925-1931. \url{ https://doi.org/10.1890/0012-9658(2006)87[1925:SIUTGO]2.0.CO;2}.

Lotwick,H.W. & Silverman, B. W. 1982. Methods for analysing spatial processes of several types of points. 
 \emph{Journal of the Royal Statistical Society} B,  44: 406-413. \url{https://doi.org/10.1111/j.2517-6161.1982.tb01221.x}.
  
Olano, J.M., Laskurain, N.A., Escudero, A. and De la Cruz, M. 2009. Why and where adult trees die in a secondary temperate forest? 
The role of neighbourhood. \emph{Annals of Forest Science}, 66: 105. \url{http://dx.doi.org/10.1051/forest:2008074}.
  
Penttinen, A. 2006. Statistics for Marked Point Patterns. In \emph{The Yearbook of the Finnish Statistical Society}, pp. 70-91. 

 \enc{Rey-Benayas, J.M., de la Montaña, E., Pérez-Camacho, L., de la Cruz, M., Moreno, D., Parejo, J.L. and Suárez-Seoane, S.}{Rey-Benayas, J.M., de la Montana, E., Perez-Camacho, L., de la Cruz, M., Moreno, D., Parejo, J.L. and Suarez-Seoane, S.} 2010. 
Inter-annual dynamics and spatial congruence of a nocturnal bird assemblage inhabiting a Mediterranean agricultural mosaic. \emph{Ardeola},57(2): 303-320.

Syrjala, S. E. 1996. A statistical test for a difference between the spatial distributions of two populations. \emph{Ecology} 77: 75-80. \url{ https://doi.org/10.2307/2265656}.

Waagepetersen, R. P. 2007. An estimating function approach to inference for inhomogeneous Neymann-Scott processes. \emph{Biometrics} 63: 252-258. \url{https://doi.org/10.1111/j.1541-0420.2006.00667.x}.

}


\keyword{ package }
\keyword{ spatial }

\examples{
\dontrun{

#############################################
### Transfom easily data from a data.frame into the ppp format 
### of spatstat:

data(fig1)

plot(fig1) #typical xyplot

fig1.ppp <- haz.ppp (fig1)

fig1.ppp

plot(fig1.ppp) # point pattern plot of spatstat



#############################################
###  Summarize the joint pattern of points and marks at different scales
###  with the normalized mark-weighted K-function (Penttinen, 2006). 
###  Compare this function in two consecutive cohorts of Helianthemum
###  squamatum seedlings:

 ## Figure 3.10 of De la Cruz (2008):
  
  data(seedlings1)
  
  data(seedlings2)
  
  s1km <- Kmm(seedlings1, r=1:100)
  
  s2km <- Kmm(seedlings2, r=1:100)
  
  plot(s1km, ylime=c(0.6,1.2), lwd=2, maine="", xlabe="r(cm)")

  plot(s2km,  lwd=2, lty=2, add=T )

  abline(h=1, lwd=2, lty=3)
  
  legend(x=60, y=1.2, legend=c("Hs_C1", "Hs_C2", "H0"),
	 lty=c(1, 2, 3), lwd=c(3, 2, 2), bty="n")
 
## A pointwise test of normalized Kmm == 1 for seedlings1:

   s1km.test <- Kmm(seedlings1, r=1:100, nsim=99)

   plot(s1km.test,  xlabe="r(cm)")
   
   
   

#############################################
###  Explore the local relationships between marks and locations (e.g. size 
###  of one cohort of H. squamatum seedlings). Map the marked point pattern 
###  to a random field for visual inspection, with the normalized mark-sum
###  measure (Penttinen, 2006).

data(seedlings1)
   
 seed.m <- marksum(seedlings1, R=25)

 plot(seed.m, what="marksum", sigma = 5)  # raw mark-sum measure; sigma is bandwith for smoothing

 plot(seed.m, what="pointsum", sigma = 5) # point sum measure
   
 plot(seed.m,  what="normalized", dimyx=200, contour=TRUE, sigma = 5) # normalized  mark-sum measure

# the same with added grid

 plot(seed.m,  what="normalized", dimyx=200, contour=TRUE, sigma = 5, grid=TRUE) 



#############################################
###  Test against the null model of "independent labelling",
###  i.e. test asociation/repulsion between  a "fixed" pattern (e.g. adult
###  H. squamatum plants) and a "variable" pattern (e.g. of surviving and 
###  dead seedlings), with 2.5\% and 97.5\% envelopes of 999 random 
###  labellings (De la Cruz & al. 2008).

data(Helianthemum)


cosa <- K012(Helianthemum, fijo="adultHS", i="deadpl", j="survpl",
             r=seq(0,200,le=201), nsim=999, nrank=25, correction="isotropic")

plot(cosa$k01, sqrt(./pi)-r~r,  col=c(3, 1, 3), lty=c(3, 1, 3), las=1,
         ylab=expression(L[12]), xlim=c(0, 200), 
	 main="adult HS vs. dead seedlings", legend=FALSE)

plot(cosa$k02, sqrt(./pi)-r~r, col=c(3, 1, 3), lty=c(3, 1, 3), las=1, 
         ylab=expression(L[12]), xlim=c(0, 200),
 	 main="adult HS vs. surviving seedlings", legend=FALSE)



#############################################
###  Test differences of agregation and segregation between two patterns, 
###  e.g. surviving and dying H. squamatum seedlings (De la Cruz & al. 2008). 

data(Helianthemum)

cosa12 <- K1K2(Helianthemum, j="deadpl", i="survpl", r=seq(0,200,le=201),
		 nsim=999, nrank=1, correction="isotropic")

plot(cosa12$k1k2, lty=c(2, 1, 2), col=c(2, 1, 2), xlim=c(0, 200),
         main= "survival- death",ylab=expression(K[1]-K[2]), legend=FALSE)

plot(cosa12$k1k12, lty=c(2, 1, 2), col=c(2, 1, 2), xlim=c(0, 200),
	 main="segregation of surviving seedlings",ylab=expression(K[1]-K[12]), legend=FALSE)

plot(cosa12$k2k12, lty=c(2, 1, 2), col=c(2, 1, 2), xlim=c(0, 200),
         main= "segregation of dying seedlings",ylab=expression(K[2]-K[12]), legend=FALSE)


#############################################
###  Test 'univariate' and 'bivariate' point patterns 
###  against non-Poisson (in-)homogeneous models 
###  (De la Cruz and Escudero, submited).

 data(urkiola)

   #####################
   ## univariate example

   # get univariate pp
   I.ppp <- split.ppp(urkiola)$birch

   # estimate inhomogeneous intensity function
   I.lam <- predict (ppm(I.ppp, ~polynom(x,y,2)), type="trend", ngrid=200)

   # Compute and plot envelopes to Kinhom, simulating from an Inhomogeneous
   #  Poisson Process:
   
   I2.env <- envelope( I.ppp,Kinhom, lambda=I.lam, correction="trans", 
                              nsim=99, simulate=expression(rpoispp(I.lam)))
   plot(I2.env, sqrt(./pi)-r~r, xlab="r (metres)", ylab= "L (r)", col=c(1,3,2,2),legend=FALSE) 

   # It seems that there is short scale clustering; let's fit an Inhomogeneous 
   # Poisson Cluster Process: 

   I.ki <- ipc.estK(mippp=I.ppp, lambda=I.lam, correction="trans")

   # Compute and plot envelopes to Kinhom, simulating from the fitted IPCP:

   Ipc.env <- Ki(I.ki, correction="trans", nsim=99, ngrid=200)

   plot (Ipc.env, xlab="r (metres)", ylab= "L (r)")
 
   #####################
   ## bivariate example: test independence between birch and quercus in Urkiola

   J.ppp <- split.ppp(urkiola)$oak
   
   # We want to simulate oak from a homogeneous Poisson model:
   J.ppm <- ppm(J.ppp, trend=~1, interaction=Poisson() )
   
   IJ.env <- Kci (mod1=I.ki, mod2=J.ppm, nsim=99)
   
   plot(IJ.env, type=12)
   
   plot(IJ.env, type=21)




#############################################
###  Simulate envelopes from the fitted values of a logistic model,
###  as in Olano et al. (2009)
   
   
   data(quercusvm)

   # read fitted values from logistic model:
   
   
   probquercus <-c(0.99955463, 0.96563477, 0.97577094, 0.97327199, 0.92437309,
   0.84023396, 0.94926682, 0.89687281, 0.99377915, 0.74157478, 0.95491518,
   0.72366493, 0.66771787, 0.77330148, 0.67569082, 0.9874892, 0.7918891, 
   0.73246803, 0.81614635, 0.66446411, 0.80077908, 0.98290508, 0.54641754,
   0.53546689, 0.73273626, 0.7347013, 0.65559655, 0.89481468, 0.63946334,
   0.62101995, 0.78996371, 0.93179582, 0.80160346, 0.82204428, 0.90050059,
   0.83810669, 0.92153079, 0.47872421, 0.24697004, 0.50680935, 0.6297911, 
   0.46374812, 0.65672284, 0.87951682, 0.35818237, 0.50932432, 0.92293014,
   0.48580241, 0.49692053, 0.52290553, 0.7317549, 0.32445982, 0.30300865,
   0.73599359, 0.6206056, 0.85777043, 0.65758613, 0.50100406, 0.31340849, 
   0.22289286, 0.40002879, 0.29567678, 0.56917817, 0.56866864, 0.27718552,
   0.4910667, 0.47394411, 0.40543788, 0.29571349, 0.30436276, 0.47859015,
   0.31754526, 0.42131675, 0.37468782, 0.73271225, 0.26786274, 0.59506388, 
   0.54801851, 0.38983575, 0.64896835, 0.37282031, 0.67624306, 0.29429766,
   0.29197755, 0.2247629, 0.40697843, 0.17022391, 0.26528042, 0.24373722,
   0.26936163, 0.13052254, 0.19958585, 0.18659692, 0.36686678, 0.47263005,
   0.39557661, 0.68048997, 0.74878567, 0.88352322, 0.93851375)
   
  

   ################################ 
   ## Envelopes for an homogeneous point pattern:
   
   cosap <- Kinhom.log(A=quercusvm, lifemark="0", nsim=99, prob=probquercus)

   plot(cosap)

   
   ################################ 
   ## Envelopes for an inhomogeneous point pattern:
   
   ## First, fit an inhomogeneous Poisson model to alive trees :
   
   quercusalive <- unmark(quercusvm[quercusvm$marks == 0])

    mod2 <- ppm(quercusalive, ~polynom(x,y,2))

    ## Now use mod2 to estimate lambda for K.inhom:
    
    cosapm <- Kinhom.log(A=quercusvm, lifemark="0", prob=probquercus,
                                   nsim=99,  mod=mod2)
    plot(cosapm)




#############################################
###  Test segregation based on the counts in the contingency table
###  of nearest neighbors in a multitype point pattern (Dixon, 2002)

data(swamp)

dixon2002(swamp,nsim=99)




#############################################
###  Fit the Poisson cluster point process to a point pattern with 
###  the method of minimum contrast (Diggle 2003).

data(gypsophylous)

# Estimate K function ("Kobs").

gyps.env <- envelope(gypsophylous, Kest, correction="iso", nsim=99)

plot(gyps.env, sqrt(./pi)-r~r, ylab="L(r)", legend=FALSE)

# Fit Poisson Cluster Process. The limits of integration 
# rmin and rmax are setup to 0 and 60, respectively. 

cosa.pc <- pc.estK(Kobs = gyps.env$obs[gyps.env$r<=60],
		           r = gyps.env$r[gyps.env$r<=60])

# Add fitted Kclust function to the plot.

lines(gyps.env$r,sqrt(Kclust(gyps.env$r, cosa.pc$sigma2,cosa.pc$rho)/pi)-gyps.env$r,
       lty=2, lwd=3, col="purple")

# A kind of pointwise test of the gypsophylous pattern been a realisation
# of the fitted model, simulating with sim.poissonc and using function J (Jest).

gyps.env.sim <- envelope(gypsophylous, Jest, nsim=99,
                    simulate=expression(sim.poissonc(gypsophylous,
		    sigma=sqrt(cosa.pc$sigma2), rho=cosa.pc$rho)))

 plot(gyps.env.sim,  main="",legendpos="bottomleft")






#############################################
###  Compute Syrjala's test for the difference between the spatial 
###  distributions of two populations, as in Rey-Benayas et al. 
### (submited)


 
   data(syr1); data(syr2); data(syr3)
   
   plot(syrjala.test(syr1, syr2, nsim=999)) 
   
   plot(syrjala.test(syr1, syr3, nsim=999)) 
   



}
}


