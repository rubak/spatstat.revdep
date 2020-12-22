#Plots_off()
oldoptions <- options(warn=-1)
oldpar <- par()
#devAskNewPage(ask = TRUE)
cat("NOTE: if you're running rstudio make sure\nthe plotting window is large otherwise you will\nget a margin error.")

open_new_plot=FALSE
#openwin_sppmix(TRUE)

#create a normal mixture
#object "normmix"
plotstring("There are 3 types of plots in sppmix:\n1) generic 2d plots\n2) 2d plots using the ggplot package\n3) 3d plots using the rgl package")

#1)Generic plots
plotstring("We illustrate generic plots first.\nWe display a point pattern\nrepresenting earthquake locations\nin CA during the year 2014.\nOnly those events with magnitudes\n>=3 in the Richter scale are displayed.")

plot2dPP(CAQuakes2014.RichterOver3.0,
         open_new_window=open_new_plot,title1="Earthquake locations in CA, 2012")

plotstring("Next we display a point pattern\nrepresenting tornado locations\nabout MO during the year 2011.")

plot2dPP(Tornadoes2011MO,open_new_window=open_new_plot,
         title1="Tornadoes about MO, USA, 2011")

#fit a Poisson with mixture intensity surface
plotstring("\nFor this tutorial we will work\nwith a known intensity and generate a\n point pattern from it. In order to\n show the full variety of plots, we\n fit the DAMCMC to this pattern.\nSee the model fitting tutorial\nfor more details on DAMCMC.")
truemix5=rnormmix(m = 5, sig0 = .1, df = 5,xlim= c(-3,3), ylim = c(-3,3))
trueintsurfmix5=to_int_surf(truemix5,lambda = 150,win =spatstat.geom::owin( c(-3,3),c(-3,3)))
ppmix5 <- rsppmix(intsurf = trueintsurfmix5)# draw points

#retrieve the surface of posterior means
plotstring("Use plot2dPP to display the point\npattern along with the true component\nmeans and the posterior means\nfrom the DAMCMC model fit.\nLabel switching is possible. See the\nhelp page for FixLS_da() for details.")
DAMCMCfit=est_mix_damcmc(pp = ppmix5, m = 5,L=10000)
SurfaceofPostMeans=GetPMEst(DAMCMCfit)

plot2dPP(ppmix5,trueintsurfmix5$mus,open_new_window=open_new_plot,
         title1="Point pattern along with the true component means")
plotstring("If we observe clustered posterior means\nthis is an indication of label switching.\nWe show how to handle this situation\nbelow and in the other demos and\ntutorials (vignettes) involving model fitting.")
plot2dPP(ppmix5,SurfaceofPostMeans$mus,open_new_window=open_new_plot,
         title1="Point pattern along with the component posterior means")

#2) ggplots
plotstring("Now we discuss 2d plots using\nthe ggplot package. Let's reproduce\nthe Earthquake and Tornadoes examples\nusing function 'PlotUSAStates'.")

plotstring("Fit an IPPP with intensity surface modeled\nby a mixture with 8 normal components.\nWe work with the California Earthquake data.")
fitDA=est_mix_damcmc(CAQuakes2014.RichterOver3.0, m=8, L = 20000)

plotstring("Now retrieve the surface of Maximum a\nPosteriori (MAP) estimates of the mixture\nparameter. Note that the resulting surface\nis not affected by label switching.")
MAPsurf=GetMAPEst(fitDA)

plotstring("Plot the states and the earthquake\nlocations along with the fitted MAP\nIPPP intensity surface.")
ret=PlotUSAStates(states=c('California','Nevada','Arizona'),
                  open_new_window=open_new_plot,                  showcentroids=FALSE,shownames=TRUE,
                  main="Earthquakes in CA, 2014",
                  pp=CAQuakes2014.RichterOver3.0,
                  surf=MAPsurf,boundarycolor="white",
                  namescolor="white")

plotstring("Repeat the process for the Tornado data.\nWe visualize the Tornado data about MO,\nUSA, by plotting the states and the\ntornado locations")
fitDA=est_mix_damcmc(Tornadoes2011MO, m=5, L = 20000)

MAPsurf=GetMAPEst(fitDA)

ret=PlotUSAStates(states=c('Iowa','Arkansas','Missouri','Illinois','Indiana','Kentucky','Tennessee','Kansas','Nebraska','Texas','Oklahoma','Mississippi','Alabama','Louisiana'),
                  open_new_window=open_new_plot,showcentroids=FALSE,shownames=TRUE,
                  plotlevels = FALSE,
                  main="Tornadoes about MO, 2011",
                  pp=Tornadoes2011MO,
                  surf=MAPsurf,boundarycolor="white",
                  namescolor="white")


plotstring("Next we illustrate standard plots of the true density\nand intensity surfaces, along with\nthe points of the point pattern.\nWe also present the contour plots.")

plot(truemix5,whichplots=0,xlim=c(-3,3),ylim=c(-3,3),open_new_window=open_new_plot,
     title1="Plot of the true mixture density with 5 bivariate normal components")+add_title("Plot of the true mixture density with 5 bivariate normal components")
plot(truemix5,whichplots=0,xlim=c(-3,3),ylim=c(-3,3),contour = TRUE,open_new_window=open_new_plot)+add_title("Contour plot of the true mixture of normals density")

plotmix_2d(trueintsurfmix5,open_new_window=open_new_plot)+add_title("True Poisson intensity surface",lambda =trueintsurfmix5$lambda,m=trueintsurfmix5$m)
plotmix_2d(trueintsurfmix5,contour = TRUE,open_new_window=open_new_plot)+add_title("Contour plot of the true intensity surface",lambda =trueintsurfmix5$lambda,m=trueintsurfmix5$m)

plotmix_2d(trueintsurfmix5,ppmix5,colors = TRUE,open_new_window=open_new_plot)+add_title("True Poisson intensity surface along with the point pattern",lambda =trueintsurfmix5$lambda,m=trueintsurfmix5$m)
plotmix_2d(trueintsurfmix5,ppmix5,colors = TRUE,contour = TRUE,open_new_window=open_new_plot)+add_title("Contour plot of the true intensity surface along with the point pattern",lambda =trueintsurfmix5$lambda,m=trueintsurfmix5$m)

plotstring("We illustrate plots from the\nDAMCMC model fit. The first surface\ncorresponds to the posterior means\nand it may suffer from label switching.\nStart with plots of the point\n pattern along with the component means.")
#plot.sppmix
plot(ppmix5,trueintsurfmix5$mus,open_new_window=open_new_plot)+add_title("Point pattern along with the true component means",lambda = trueintsurfmix5$lambda,m=trueintsurfmix5$m,n=ppmix5$n)
plot(ppmix5,SurfaceofPostMeans$mus,open_new_window=open_new_plot)+add_title("Point pattern along with the posterior means of the mixture components",lambda = SurfaceofPostMeans$lambda,m=SurfaceofPostMeans$m,n=ppmix5$n,L=.9*DAMCMCfit$L)

#plot the surface and the point pattern
plotmix_2d(SurfaceofPostMeans,ppmix5,open_new_window=open_new_plot)+add_title("Poisson surface of posterior means",lambda =SurfaceofPostMeans$lambda,m=SurfaceofPostMeans$m)
plotmix_2d(SurfaceofPostMeans,contour = TRUE,ppmix5,open_new_window=open_new_plot)+add_title("Contours of the Poisson surface of posterior means",lambda =SurfaceofPostMeans$lambda,m=SurfaceofPostMeans$m)

# Plot the average of the surfaces of the posterior realizations
plotstring("\nThe second surface is the average\nof the surfaces for each posterior\nrealization and does not suffer\nfrom the label switching problem.\nComputation of the average of surfaces\ncan be slow depending on parameter\nchoices and number of posterior realizations.")

plotstring("\nDone with calculations. Parameter\nLL affects the grid size\nof the function plot_avgsurf.\nNote that the surface is also\nreturned by this function.\nWe expect a much better fit.\nLet's check the average of the surfaces.")
avgsurf=plot_avgsurf(DAMCMCfit, LL = 100,showplot =FALSE)


p<-plot_density(as.data.frame(avgsurf))+ggplot2::ggtitle("Average surface of the posterior realization surfaces\nWindow=[-3,3]x[-3,3], x denotes a true component mean")
#show the point pattern points
pp_df <- data.frame(ppmix5$x,ppmix5$y)
names(pp_df) <- c("x", "y")
p<-p + ggplot2::geom_point(data = pp_df,size=0.8)
#show the true means
mean_df <- data.frame(do.call(rbind, trueintsurfmix5$mus))
names(mean_df) <- c("x", "y")
p + ggplot2::geom_point(data = mean_df, color = "red",shape = "x", size = 5)

p<-plot_density(as.data.frame(avgsurf),contour =TRUE)+ggplot2::ggtitle("Contours of the average surface of the posterior realization surfaces\nWindow=[-3,3]x[-3,3], x denotes a true component mean")
#show the point pattern points
pp_df <- data.frame(ppmix5$x,ppmix5$y)
names(pp_df) <- c("x", "y")
p<-p + ggplot2::geom_point(data = pp_df,size=0.8)
#show the true means
mean_df <- data.frame(do.call(rbind, trueintsurfmix5$mus))
names(mean_df) <- c("x", "y")
p + ggplot2::geom_point(data = mean_df, color = "red",shape = "x", size = 5)

plotstring("Plot the chains for the\nmixture parameters ps, mus, and\nsigmas. We can do them in separate\nplots or in the same plot\nwhich is useful in graphically\nassessing if the label switching\nproblem is present.")
#plot the chains
plot_chains(DAMCMCfit,open_new_window=open_new_plot)
plot_chains(DAMCMCfit,open_new_window=open_new_plot,separate = FALSE)

#check for label switching
plotstring("\nNow use function check_labels()\nto see if the label switching\nproblem is present.")
labelswitch=check_labels(DAMCMCfit)
if(labelswitch)
{
  plotstring("\nLabel switching is present. Use\nfunction FixLS_da() to\nfix the problem. First use an\nidentifiability constraint to\npermute the labels.")
  post_fixedDAMCMCfitIC = FixLS_da(DAMCMCfit)
  plot_chains(post_fixedDAMCMCfitIC,open_new_window=open_new_plot,separate = FALSE)

  plotstring("\nNow use a decision theoretic approach\nbased on minimization of the SEL\n(Squared Error Loss), in order to find\nthe best permutation, and apply it\nto all the posterior realizations.")

  plotstring("\nCalculations can take a long time. We expect\nthat the chains will look better now.")
  post_fixedDAMCMCfitSEL = FixLS_da(DAMCMCfit,approx=FALSE)

  plot_chains(post_fixedDAMCMCfitSEL,open_new_window=open_new_plot,separate = FALSE)
}
if(!labelswitch)
  plotstring("\nNo label switching detected. The\nsurface of posterior means\nis a valid estimator of the true\nPoisson intensity surface.")

plotstring("We also plot the posterior means\nof the membership indicator variables\n(allocation indicators of a point to a\ncomponent). These are the posterior\nprobabilities of a point belonging\nto a specific component.")
print(plot_ind(DAMCMCfit,open_new_window=open_new_plot)+add_title("Posterior means of the membership indicators (Unpermuted labels)", m = DAMCMCfit$m, n = DAMCMCfit$data$n))

if(labelswitch)
{
  plotstring("Membership indicator variables\nbased on an identifiability\nconstraint (IC).")
  print(plot_ind(post_fixedDAMCMCfitIC,open_new_window=open_new_plot)+add_title("Posterior means of the membership indicators (IC permuted labels)", m = post_fixedDAMCMCfitIC$m, n = post_fixedDAMCMCfitIC$data$n))

  plotstring("Membership indicator variables\nbased on the best permutation.")
  print(plot_ind(post_fixedDAMCMCfitSEL,open_new_window=open_new_plot)+add_title("Posterior means of the membership indicators (Permuted labels)", m = post_fixedDAMCMCfitSEL$m, n = post_fixedDAMCMCfitSEL$data$n))

}

plotstring("Finally, we produce 3d plots\nof the intensity surfaces we discussed.")

plot(truemix5,whichplots=1,xlim=c(-3,3),ylim=c(-3,3),title1="True normal mixture density with 5 components")
plot(trueintsurfmix5,main="True normal mixture intensity surface with 5 components")

plotmix_3d(avgsurf, title1 = paste("Average of",.9*DAMCMCfit$L,"posterior realizations"))

plot(SurfaceofPostMeans, main = "Poisson surface of posterior means")

if(labelswitch)
{
  permSurfaceofPostMeansIC=GetPMEst(post_fixedDAMCMCfitIC)
  plot(permSurfaceofPostMeansIC, main = "Poisson surface of posterior means (IC)")

  permSurfaceofPostMeansSEL=GetPMEst(post_fixedDAMCMCfitSEL)
  plot(permSurfaceofPostMeansSEL, main = "Poisson surface of posterior means (Permuted)")
}

plotstring("Note that this analysis can be\nrepeated with truncate=TRUE in every function\ncalled, in order to handle edge effects.\nSee the model fitting demos\nfor more details.\nThanks")

suppressWarnings( par(oldpar))
options(oldoptions)
#on.exit(devAskNewPage(ask = FALSE))
#par(ask=FALSE)
