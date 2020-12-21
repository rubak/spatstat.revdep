#Plots_off()
oldoptions <- options(warn=-1)
oldpar <- par()
open_new_plot=FALSE
cat("NOTE: if you're running rstudio make sure\nthe plotting window is large otherwise you will\nget a margin error.")

plotstring("\nWe discuss the BDMCMC Bayesian approach\nfor fitting a Poisson point process with an\nintensity surface that is a mixture of normal\ncomponents. Assume that there are edge\neffects and entertain a random number of\nmixture components (i.e.,\nm is a parameter of the model).")

plotstring("We generate a mixture with\nm=4 components and build the\nPoisson intensity surface. The\noriginal window used is [-2,2]x[-2,2]")
truemix4=rnormmix(m = 4, sig0 = .1, df = 5,xlim= c(-2,2), ylim = c(-2,2))
plot(truemix4,xlim= c(-2,2), ylim = c(-2,2),whichplots=0,open_new_window=open_new_plot)+add_title("True mixture of normals density")
trueintsurfmix4=to_int_surf(truemix4,lambda = 150,win =spatstat::owin( c(-2,2),c(-2,2)))
bigwin=spatstat::owin(c(-2,2),c(-2,2))

plotstring("Then based on a point pattern from this model,\nour goal is to estimate the parameters of the\nmixture model (including the number of\ncomponents) and lambda, in order to\nrecover the true Poisson intensity surface.")

#analyzing a point pattern when we know the true Poisson surface
ppmix4 <- rsppmix(intsurf = trueintsurfmix4)# draw points
plotstring("\nThe truncate parameter is TRUE by default\nin the rsppmix function, in order\nto generate a point pattern on the\nspecified window from the Poisson process\nwith the created intensity surface.")

plotstring("The window of observation has a major\neffect in calculations in this case,\nespecially if we have a lot of events\nnear the boundary. Parameter truncate\nis set to TRUE in all the functions\ncalled below. First we display the pattern.")
plotmix_2d(trueintsurfmix4,ppmix4,colors = TRUE,open_new_window=open_new_plot,win=trueintsurfmix4$window)+add_title("True Poisson intensity surface along with the point pattern, W=[-2,2]x[-2,2]",lambda =trueintsurfmix4$lambda,m=trueintsurfmix4$m,n=ppmix4$n)

plotstring("\nThe function est_mix_bdmcmc can be\nused to fit the BDMCMC with a random\nnumber of mixture components. We\ncall the routine with the truncate\nparameter set to TRUE (edge effects)\nand entertain mixtures with number of\ncomponents in the range of 1 to 10.")

plotstring("\nDone with calculations. Let us check\nout the distribution of the number\nof components and the trace plot.\nAlso note the frequency table for the\nnumbers of components that were visited\nby the BDMCMC chain.")

#fit a Poisson with mixture intensity surface
BDMCMCfit=est_mix_bdmcmc(pp = ppmix4,
        m = 10,L=20000,truncate = TRUE)

plot_CompDist(BDMCMCfit,open_new_window=open_new_plot)

plotstring("We anticipate that the best estimator\nof the true Poisson intensity surface will\nbe the Bayesian model average (BMA) of all\nrealizations (after burnin). Use function GetBMA\nin order to obtain and plot the BMA surface.\nCalculation of the BMA can be slow...")

plotstring("\nDone with calculations. Let us check\nout the BMA surface.")
BDMCMCfit_BMA=GetBMA(BDMCMCfit)

plot_density(as.data.frame(BDMCMCfit_BMA))+ggplot2::ggtitle("Bayesian model average intensity surface")

plot_density(as.data.frame(BDMCMCfit_BMA),TRUE)+ggplot2::ggtitle("Contours of the Bayesian model average intensity surface")

plotstring("\nThe BMA intensity surface does not\nallow us to achieve mixture deconvolution,\ni.e., obtain specific values for the\nparameters of the mixture model,\nincluding the exact number of components,\nas well as, the corresponding probabilities,\nmeans and covariance matrices.")

plotstring("\nThis can be crucial in some applications,\ne.g., when the component means represent\nthe locations of monitoring stations, we\nneed to explicitly estimate the values of the\nmixture component means.")

plotstring("\nAs a result, researchers can make\ndecisions about where to allocate\nresources more efficiently.")

BDtab=GetBDTable(BDMCMCfit,FALSE)#retrieve frequency table and MAP estimate for number of components

MAPm=BDtab$MAPcomp

plotstring("\nThe maximum a posteriori (MAP) estimate\nfor the number of components, is a natural\nchoice that allows us to achieve mixture\ndeconvolution. Use functions GetBDTable and\nGetBDCompfit to retrieve the MAP and\nthe corresponding realizations from\nthe BDMCMC fit.")

BDMCMCfit=drop_realization(BDMCMCfit,.1*BDMCMCfit$L)

BDMCMCfit=drop_realization(BDMCMCfit,
                           (BDMCMCfit$Badgen==1))

#retrieve all BDMCMC realizations corresponding to mixture with MAP components
BDMCMCfitMAPcomp=GetBDCompfit(BDMCMCfit,MAPm,burnin=0)
BDMCMCfitMAPcomp
BDMCMCfitMAPcompgens=BDMCMCfitMAPcomp$BDgens

plotmix_2d(BDMCMCfitMAPcomp$BDsurf,ppmix4,colors = TRUE,open_new_window=open_new_plot,win=bigwin)+add_title("MAP Poisson intensity surface along with the point pattern",lambda =BDMCMCfitMAPcomp$BDsurf$lambda,m=BDMCMCfitMAPcomp$BDsurf$m,n=ppmix4$n,L=BDMCMCfitMAPcomp$BDgens$L)

plotstring("\nWe'll work more with the intensity surface\ncorresponding to the MAP # of components.\nPlot the chains for the mixture parameters\nps, mus, and sigmas, in their own plots\nfor all components. This is useful in graphically\nassessing if the label switching\nproblem is present in the MAP intensity surface.")

#plot the chains
plot_chains(BDMCMCfitMAPcompgens,open_new_window=open_new_plot,separate = FALSE)

#check for label switching
plotstring("\nNow we use the function check_labels()\nto check if the label switching\nproblem is present.")
labelswitch=check_labels(BDMCMCfitMAPcompgens)

if(labelswitch)
{
  plotstring("\nLabel switching is present. Use\nfunction FixLS_da() to\nfix the problem. First use an\nidentifiability constraint to\npermute the labels.")
  post_fixedBDMCMCfitIC = FixLS_da(BDMCMCfitMAPcompgens,burnin=0)
  plot_chains(post_fixedBDMCMCfitIC,open_new_window=open_new_plot,separate = FALSE)

  plotstring("Membership indicator variables\nbased on an identifiability\nconstraint (IC).")
  print(plot_ind(post_fixedBDMCMCfitIC,open_new_window=open_new_plot)+add_title("Posterior means of the membership indicators (IC permuted labels)", m = post_fixedBDMCMCfitIC$m, n = post_fixedBDMCMCfitIC$data$n))

  #plot the surface and the point pattern
  permSurfaceofPostMeansIC=GetPMEst(post_fixedBDMCMCfitIC,burnin=0)

  plotstring("The estimated Poisson surface of\nposterior means based on the\nidentifiability constraint follows.")
  print(plotmix_2d(permSurfaceofPostMeansIC,ppmix4,colors = TRUE,open_new_window=open_new_plot,win=bigwin)+add_title("Poisson surface of posterior means (IC)",lambda =permSurfaceofPostMeansIC$lambda,m=permSurfaceofPostMeansIC$m,n=ppmix4$n,L=post_fixedBDMCMCfitIC$L))

  print(plotmix_2d(permSurfaceofPostMeansIC,ppmix4,colors = TRUE,contour = TRUE,open_new_window=open_new_plot,win=bigwin)+add_title("Contours of the Poisson surface of posterior means (IC)",lambda =permSurfaceofPostMeansIC$lambda,m=permSurfaceofPostMeansIC$m,n=ppmix4$n,L=post_fixedBDMCMCfitIC$L))

  plotstring("Now use a decision theoretic approach\nbased on minimization of the SEL\n(Squared Error Loss), is order to find\nthe best permutation, and apply it\nto all the posterior realizations.")

  plotstring("\nCalculations can take a long time.\nWe expect the chains to look better now.")
  post_fixedBDMCMCfitSEL = FixLS_da(BDMCMCfitMAPcompgens,approx=FALSE,burnin=0)
  plot_chains(post_fixedBDMCMCfitSEL,open_new_window=open_new_plot,separate = FALSE)

  plotstring("Membership indicator variables\nbased on the best permutation.")
  print(plot_ind(post_fixedBDMCMCfitSEL,open_new_window=open_new_plot)+add_title("Posterior means of the membership indicators (Permuted labels)", m = post_fixedBDMCMCfitSEL$m, n = post_fixedBDMCMCfitSEL$data$n))

  #plot the surface and the point pattern
  permSurfaceofPostMeansSEL=GetPMEst(post_fixedBDMCMCfitSEL,burnin=0)

  plotstring("The estimated Poisson surface of\nposterior means based on the\nbest permutation follows.")
  print(plotmix_2d(permSurfaceofPostMeansSEL,ppmix4,colors = TRUE,open_new_window=open_new_plot,win=bigwin)+add_title("Poisson surface of posterior means (Permuted labels)",lambda =permSurfaceofPostMeansSEL$lambda,m=permSurfaceofPostMeansSEL$m,n=ppmix4$n,L=post_fixedBDMCMCfitSEL$L))

  print(plotmix_2d(permSurfaceofPostMeansSEL,ppmix4,colors = TRUE,contour = TRUE,open_new_window=open_new_plot,win=bigwin)+add_title("Contours of the Poisson surface of posterior means (Permuted labels)",lambda =permSurfaceofPostMeansSEL$lambda,m=permSurfaceofPostMeansSEL$m,n=ppmix4$n,L=post_fixedBDMCMCfitSEL$L))

}

if(!labelswitch)
{
  plotstring("\nNo label switching detected. The\nsurface of posterior means\nis a valid estimator of the true\nPoisson intensity surface.")

  SurfaceofPostMeans=GetPMEst(BDMCMCfitMAPcompgens,burnin=0)

  print(plot(ppmix4,trueintsurfmix4$mus,open_new_window=open_new_plot)+add_title("Point pattern along with the true component means",lambda = trueintsurfmix4$lambda,m=trueintsurfmix4$m,n=ppmix4$n))

  print(plot(ppmix4,SurfaceofPostMeans$mus,open_new_window=open_new_plot)+add_title("Point pattern along with the posterior means of the mixture components",lambda = SurfaceofPostMeans$lambda,m=SurfaceofPostMeans$m,n=ppmix4$n,L=BDMCMCfitMAPcompgens$L))

  #plot the surface and the point pattern
  print(plotmix_2d(SurfaceofPostMeans,ppmix4,colors = TRUE,open_new_window=open_new_plot,win=bigwin)+add_title("Poisson surface of posterior means",lambda =SurfaceofPostMeans$lambda,m=SurfaceofPostMeans$m,n=ppmix4$n,L=.9*BDMCMCfit$L))

  print(plotmix_2d(SurfaceofPostMeans,ppmix4,colors = TRUE,contour = TRUE,open_new_window=open_new_plot,win=bigwin)+add_title("Contours of the Poisson surface of posterior means",lambda =SurfaceofPostMeans$lambda,m=SurfaceofPostMeans$m,n=ppmix4$n,L=.9*BDMCMCfit$L))

  plotstring("We also plot the posterior means\nof the membership indicator variables\n(allocation indicators of a point to a\ncomponent). These are the posterior\nprobabilities of a point belonging\nto a specific component.")

  print(plot_ind(BDMCMCfitMAPcompgens,open_new_window=open_new_plot)+add_title("Posterior means of the membership indicators (Unpermuted labels)", m = BDMCMCfitMAPcompgens$m, n = BDMCMCfitMAPcompgens$data$n))
}

# Plot the average of the surfaces of the posterior realizations
plotstring("\nThe average of the surfaces for\neach posterior realization (for the MAP) does\nnot suffer from the label switching\nproblem. Computation of the average of\nthe surfaces can be slow.")

plotstring("\nDone with calculations. Parameter LL\naffects the grid size of the function\nplot_avgsurf. Note that the surface is also\nreturned by this function.\nWe expect a much better fit.\nLet's check the average of the surfaces.")
avgsurf=plot_avgsurf(BDMCMCfitMAPcompgens,win=bigwin, LL = 100,showplot =FALSE)

p<-plot_density(as.data.frame(avgsurf))+ggplot2::ggtitle("Average surface of the MAP posterior realization surfaces\nWindow=[-4,4]x[-4,4], x denotes a true component mean")
#show the point pattern points
pp_df <- data.frame(ppmix4$x,ppmix4$y)
names(pp_df) <- c("x", "y")
p<-p + ggplot2::geom_point(data = pp_df,size=0.8)
#show the true means
mean_df <- data.frame(do.call(rbind, trueintsurfmix4$mus))
names(mean_df) <- c("x", "y")
p + ggplot2::geom_point(data = mean_df, color = "red",shape = "x", size = 5)

p<-plot_density(as.data.frame(avgsurf),contour =TRUE)+ggplot2::ggtitle("Contours of the average surface of the MAP posterior realization surfaces\nWindow=[-4,4]x[-4,4], x denotes a true component mean")
#show the point pattern points
pp_df <- data.frame(ppmix4$x,ppmix4$y)
names(pp_df) <- c("x", "y")
p<-p + ggplot2::geom_point(data = pp_df,size=0.8)
#show the true means
mean_df <- data.frame(do.call(rbind, trueintsurfmix4$mus))
names(mean_df) <- c("x", "y")
p + ggplot2::geom_point(data = mean_df, color = "red",shape = "x", size = 5)

plotstring("Finally, we produce 3d plots\nof the intensity surfaces we discussed.")
#devAskNewPage(ask = FALSE)

plot(truemix4,whichplots=1,xlim=bigwin$xrange,ylim=bigwin$yrange,title1="True normal mixture density with 4 components")
plot(trueintsurfmix4,main="True normal mixture intensity surface with 4 components",win=bigwin)

burnin=.1*BDMCMCfit$L
title1 = paste("Bayesian model average of",floor(BDMCMCfit$L-burnin),"posterior realizations")
plotmix_3d(BDMCMCfit_BMA,title1=title1)

plot(BDMCMCfitMAPcomp$BDsurf,main=paste("MAP Mixture intensity surface with",MAPm,"components"),win=bigwin)


plotmix_3d(avgsurf, title1 = paste("Average surface of",BDMCMCfitMAPcompgens$L,"posterior surfaces (for MAP m)"))

if(!labelswitch)
{
  print(plot(SurfaceofPostMeans, main = "Poisson surface of posterior means (MAP components)",win=bigwin))
}

if(labelswitch)
{
  permSurfaceofPostMeansIC=GetPMEst(post_fixedBDMCMCfitIC,burnin=0)
  plot(permSurfaceofPostMeansIC,win=bigwin, main = "Posterior surfacePoisson surface of posterior means (IC)")

  permSurfaceofPostMeansSEL=GetPMEst(post_fixedBDMCMCfitSEL,burnin=0)
  plot(permSurfaceofPostMeansSEL,win=bigwin, main = "Poisson surface of posterior means (Permuted)")
}

plotstring("See the vignettes and help pages of\nthese functions for more details.\nThanks")

suppressWarnings( par(oldpar))
options(oldoptions)
