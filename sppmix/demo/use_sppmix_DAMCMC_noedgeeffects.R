#Plots_off()
oldoptions <- options(warn=-1)
oldpar <- par()
open_new_plot=FALSE
cat("NOTE: if you're running rstudio make sure\nthe plotting window is large otherwise you will\nget a margin error.")

plotstring("We discuss the DAMCMC Bayesian approach\nfor fitting a Poisson point process\nwith an intensity surface that is a\nmixture of normal components.\nAssume that there are no edge effects\nand a fixed number of mixture components.")

plotstring("We generate a mixture with\n4 components and build the\nPoisson intensity surface. The\noriginal window used is [-2,2]x[-2,2]")
truemix4=rnormmix(m = 4, sig0 = .1, df = 5,xlim= c(-2,2), ylim = c(-2,2))
plot(truemix4,xlim= c(-2,2), ylim = c(-2,2),whichplots=0,open_new_window=open_new_plot)+add_title("True mixture of normals density")
trueintsurfmix4=to_int_surf(truemix4,lambda = 150,win =spatstat.geom::owin( c(-2,2),c(-2,2)))

plotstring("Then based on a point pattern\nfrom this model, our goal is\nto estimate the parameters of the mixture\nmodel and lambda, in order to recover\nthe true Poisson intensity surface.")

#analyzing a point pattern when we know the true Poisson surface
ppmix4 <- rsppmix(intsurf = trueintsurfmix4,truncate = FALSE)# draw points
plotstring("We set the truncate parameter to FALSE\nin the rsppmix function, in order\nto generate a point pattern on the\nwhole plane from the Poisson process\nwith the created intensity surface.")

plotstring("The window of observation has no\neffect in calculations in this case.\nParameter truncate is set to FALSE in\nall functions called. We display the\npattern over several windows.")
plotmix_2d(trueintsurfmix4,ppmix4,colors = TRUE,open_new_window=open_new_plot,win=spatstat.geom::owin(c(-1,1),c(-1,1)))+add_title("True Poisson intensity surface along with the point pattern, W=[-1,1]x[-1,1]",lambda =trueintsurfmix4$lambda,m=trueintsurfmix4$m,n=ppmix4$n)
plotmix_2d(trueintsurfmix4,ppmix4,colors = TRUE,open_new_window=open_new_plot,win=trueintsurfmix4$window)+add_title("True Poisson intensity surface along with the point pattern, W=[-2,2]x[-2,2]",lambda =trueintsurfmix4$lambda,m=trueintsurfmix4$m,n=ppmix4$n)
plotmix_2d(trueintsurfmix4,ppmix4,colors = TRUE,open_new_window=open_new_plot,win=spatstat.geom::owin(c(-3,3),c(-3,3)))+add_title("True Poisson intensity surface along with the point pattern, W=[-3,3]x[-3,3]",lambda =trueintsurfmix4$lambda,m=trueintsurfmix4$m,n=ppmix4$n)
plotmix_2d(trueintsurfmix4,ppmix4,colors = TRUE,open_new_window=open_new_plot,win=spatstat.geom::owin(c(-4,4),c(-4,4)))+add_title("True Poisson intensity surface along with the point pattern, W=[-4,4]x[-4,4]",lambda =trueintsurfmix4$lambda,m=trueintsurfmix4$m,n=ppmix4$n)
bigwin=spatstat.geom::owin(c(-4,4),c(-4,4))

plotstring("\nThe function est_mix_damcmc can be\nused to fit the DAMCMC for a fixed\nnumber of mixture components. We\ncall the routine with the truncate\nparameter set to FALSE (no edge effects)\nand entertain a mixture with m=4 components.")

plotstring("\nDone with calculations. Now we check for\nlabel switching. If present we permute the\nlabels and use the surface of posterior\nmeans to estimate the true intensity. If there\nis no label switching we do not\nhave to permute the posterior realizations\nfrom the DAMCMC fit.")

#fit a Poisson with mixture intensity surface
DAMCMCfit=est_mix_damcmc(pp = ppmix4, m = 4,L=10000)

plotstring("Plot the chains for the\nmixture parameters ps, mus, and\nsigmas. We can do them in the same plot\nwhich is useful in graphically\nassessing if the label switching\nproblem is present.")
#plot the chains
plot_chains(DAMCMCfit,open_new_window=open_new_plot,separate = FALSE)

#check for label switching
plotstring("\nNow we use the function check_labels()\nto check if the label switching\nproblem is present.")
labelswitch=check_labels(DAMCMCfit)

if(labelswitch)
{
  plotstring("\nLabel switching is present. Use\nfunction FixLS_da() to\nfix the problem. First use an\nidentifiability constraint to\npermute the labels.")
  post_fixedDAMCMCfitIC = FixLS_da(DAMCMCfit)
  plot_chains(post_fixedDAMCMCfitIC,open_new_window=open_new_plot,separate = FALSE)

  plotstring("Membership indicator variables\nbased on an identifiability\nconstraint (IC).")
  print(plot_ind(post_fixedDAMCMCfitIC,open_new_window=open_new_plot)+add_title("Posterior means of the membership indicators (IC permuted labels)", m = post_fixedDAMCMCfitIC$m, n = post_fixedDAMCMCfitIC$data$n))

  #plot the surface and the point pattern
  permSurfaceofPostMeansIC=GetPMEst(post_fixedDAMCMCfitIC)

  plotstring("The estimated Poisson surface of\nposterior means based on the\nidentifiability constraint follows.")
  print(plotmix_2d(permSurfaceofPostMeansIC,ppmix4,colors = TRUE,open_new_window=open_new_plot,win=bigwin)+add_title("Poisson surface of posterior means (IC)",lambda =permSurfaceofPostMeansIC$lambda,m=permSurfaceofPostMeansIC$m,n=ppmix4$n,L=.9*DAMCMCfit$L))

  print(plotmix_2d(permSurfaceofPostMeansIC,ppmix4,colors = TRUE,contour = TRUE,open_new_window=open_new_plot,win=bigwin)+add_title("Contours of the Poisson surface of posterior means (IC)",lambda =permSurfaceofPostMeansIC$lambda,m=permSurfaceofPostMeansIC$m,n=ppmix4$n,L=.9*DAMCMCfit$L))

  plotstring("Now use a decision theoretic approach\nbased on minimization of the SEL\n(Squared Error Loss), is order to find\nthe best permutation, and apply it\nto all the posterior realizations.")

  plotstring("\nCalculations can take a long time.\nWe expect the chains to look better now.")
  post_fixedDAMCMCfitSEL = FixLS_da(DAMCMCfit,approx=FALSE)
  plot_chains(post_fixedDAMCMCfitSEL,open_new_window=open_new_plot,separate = FALSE)

  plotstring("Membership indicator variables\nbased on the best permutation.")
  print(plot_ind(post_fixedDAMCMCfitSEL,open_new_window=open_new_plot)+add_title("Posterior means of the membership indicators (Permuted labels)", m = post_fixedDAMCMCfitSEL$m, n = post_fixedDAMCMCfitSEL$data$n))

  #plot the surface and the point pattern
  permSurfaceofPostMeansSEL=GetPMEst(post_fixedDAMCMCfitSEL)

  plotstring("The estimated Poisson surface of\nposterior means based on the\nbest permutation follows.")
  print(plotmix_2d(permSurfaceofPostMeansSEL,ppmix4,colors = TRUE,open_new_window=open_new_plot,win=bigwin)+add_title("Poisson surface of posterior means (Permuted labels)",lambda =permSurfaceofPostMeansSEL$lambda,m=permSurfaceofPostMeansSEL$m,n=ppmix4$n,L=.9*DAMCMCfit$L))

  print(plotmix_2d(permSurfaceofPostMeansSEL,ppmix4,colors = TRUE,contour = TRUE,open_new_window=open_new_plot,win=bigwin)+add_title("Contours of the Poisson surface of posterior means (Permuted labels)",lambda =permSurfaceofPostMeansSEL$lambda,m=permSurfaceofPostMeansSEL$m,n=ppmix4$n,L=.9*DAMCMCfit$L))

}

if(!labelswitch)
{
  plotstring("\nNo label switching detected. The\nsurface of posterior means\nis a valid estimator of the true\nPoisson intensity surface.")

  SurfaceofPostMeans=GetPMEst(DAMCMCfit)

  print(plot(ppmix4,trueintsurfmix4$mus,open_new_window=open_new_plot)+add_title("Point pattern along with the true component means",lambda = trueintsurfmix4$lambda,m=trueintsurfmix4$m,n=ppmix4$n))

  print(plot(ppmix4,SurfaceofPostMeans$mus,open_new_window=open_new_plot)+add_title("Point pattern along with the posterior means of the mixture components",lambda = SurfaceofPostMeans$lambda,m=SurfaceofPostMeans$m,n=ppmix4$n,L=.9*DAMCMCfit$L))

  #plot the surface and the point pattern
  print(plotmix_2d(SurfaceofPostMeans,ppmix4,colors = TRUE,open_new_window=open_new_plot,win=bigwin)+add_title("Poisson surface of posterior means",lambda =SurfaceofPostMeans$lambda,m=SurfaceofPostMeans$m,n=ppmix4$n,L=.9*DAMCMCfit$L))

  print(plotmix_2d(SurfaceofPostMeans,ppmix4,colors = TRUE,contour = TRUE,open_new_window=open_new_plot,win=bigwin)+add_title("Contours of the Poisson surface of posterior means",lambda =SurfaceofPostMeans$lambda,m=SurfaceofPostMeans$m,n=ppmix4$n,L=.9*DAMCMCfit$L))

  plotstring("We also plot the posterior means\nof the membership indicator variables\n(allocation indicators of a point to a\ncomponent). These are the posterior\nprobabilities of a point belonging\nto a specific component.")

  print(plot_ind(DAMCMCfit,open_new_window=open_new_plot)+add_title("Posterior means of the membership indicators (Unpermuted labels)", m = DAMCMCfit$m, n = DAMCMCfit$data$n))
}

# Plot the average of the surfaces of the posterior realizations
plotstring("\nThe average of the surfaces for\neach posterior realization does\nnot suffer from the label switching\nproblem. Computation of the average of\nthe surfaces can be slow.")

plotstring("\nDone with calculations. Parameter LL\naffects the grid size of the function\nplot_avgsurf. Note that the surface is also\nreturned by this function.\nWe expect a much better fit.\nLet's check the average of the surfaces.")
avgsurf=plot_avgsurf(DAMCMCfit,win=bigwin, LL = 100,showplot =FALSE)

p<-plot_density(as.data.frame(avgsurf))+ggplot2::ggtitle("Average surface of the posterior realization surfaces\nWindow=[-4,4]x[-4,4], x denotes a true component mean")
#show the point pattern points
pp_df <- data.frame(ppmix4$x,ppmix4$y)
names(pp_df) <- c("x", "y")
p<-p + ggplot2::geom_point(data = pp_df,size=0.8)
#show the true means
mean_df <- data.frame(do.call(rbind, trueintsurfmix4$mus))
names(mean_df) <- c("x", "y")
p + ggplot2::geom_point(data = mean_df, color = "red",shape = "x", size = 5)

p<-plot_density(as.data.frame(avgsurf),contour =TRUE)+ggplot2::ggtitle("Contours of the average surface of the posterior realization surfaces\nWindow=[-4,4]x[-4,4], x denotes a true component mean")
#show the point pattern points
pp_df <- data.frame(ppmix4$x,ppmix4$y)
names(pp_df) <- c("x", "y")
p<-p + ggplot2::geom_point(data = pp_df,size=0.8)
#show the true means
mean_df <- data.frame(do.call(rbind, trueintsurfmix4$mus))
names(mean_df) <- c("x", "y")
p + ggplot2::geom_point(data = mean_df, color = "red",shape = "x", size = 5)

plotstring("Finally, we produce 3d plots\nof the intensity surfaces we discussed.")
devAskNewPage(ask = FALSE)

plot(truemix4,whichplots=1,xlim=bigwin$xrange,ylim=bigwin$yrange,title1="True normal mixture density with 4 components")
plot(trueintsurfmix4,main="True normal mixture intensity surface with 4 components",win=bigwin)

plotmix_3d(avgsurf, title1 = paste("Average surface of",.9*DAMCMCfit$L,"posterior surfaces"))

if(!labelswitch)
{
  print(plot(SurfaceofPostMeans, main = "Poisson surface of posterior means",win=bigwin))
}

if(labelswitch)
{
  permSurfaceofPostMeansIC=GetPMEst(post_fixedDAMCMCfitIC)
  plot(permSurfaceofPostMeansIC, main = "Poisson surface of posterior means (IC)",win=bigwin)

  permSurfaceofPostMeansSEL=GetPMEst(post_fixedDAMCMCfitSEL)
  plot(permSurfaceofPostMeansSEL, main = "Poisson surface of posterior means (Permuted)",win=bigwin)
}

plotstring("See the vignettes and help pages of\nthese functions for more details.\nNote that this analysis can be\nrepeated with truncate=TRUE in every function\ncalled, in order to handle edge effects.\nThanks")

suppressWarnings( par(oldpar))
options(oldoptions)
