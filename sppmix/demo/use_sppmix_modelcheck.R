#Plots_off()
oldoptions <- options(warn=-1)
oldpar <- par()
open_new_plot=FALSE

plotstring("\nWe briefly discuss model selection\nand checking. As we saw in demos 5 and 6\nthe BDMCMC can be used to obtain the\ndistribution of the number of components\nand the MAP number of components.")

plotstring("\nWe show how to perform a Monte Carlo\ngoodness of fit test as well as, compute\nthe AIC, BIC and ICLC model\nselection criteria. The function selectMix\ncomputes the criteria whereas the function\nmc_gof performs the MC gof test.")

plotstring("\nIn addition, function kstest2d provides the\nKolmogorov-Smirnov non-parametric gof test.\nWe will assume no edge effects for this demo.")

plotstring("We generate a mixture with\n4 components and build the\nPoisson intensity surface. The\noriginal window used is [-2,2]x[-2,2]")
truemix4=rnormmix(m = 4, sig0 = .1, df = 5,xlim= c(-2,2), ylim = c(-2,2))
trueintsurfmix4=to_int_surf(truemix4,lambda = 150,win =spatstat::owin( c(-2,2),c(-2,2)))
#analyzing a point pattern when we know the true Poisson surface
ppmix4 <- rsppmix(intsurf = trueintsurfmix4,truncate = FALSE)# draw points
plotstring("We set the truncate parameter to FALSE\nin the rsppmix function, in order\nto generate a point pattern on the\nwhole plane from the Poisson process\nwith the created intensity surface.")

plotstring("The window of observation has no\neffect in calculations in this case.\nParameter truncate is set to FALSE in\nall functions called. We display the\npattern over several windows.")
plotmix_2d(trueintsurfmix4,ppmix4,colors = TRUE,open_new_window=open_new_plot,win=spatstat::owin(c(-1,1),c(-1,1)))+add_title("True Poisson intensity surface along with the point pattern, W=[-1,1]x[-1,1]",lambda =trueintsurfmix4$lambda,m=trueintsurfmix4$m,n=ppmix4$n)
plotmix_2d(trueintsurfmix4,ppmix4,colors = TRUE,open_new_window=open_new_plot,win=trueintsurfmix4$window)+add_title("True Poisson intensity surface along with the point pattern, W=[-2,2]x[-2,2]",lambda =trueintsurfmix4$lambda,m=trueintsurfmix4$m,n=ppmix4$n)
plotmix_2d(trueintsurfmix4,ppmix4,colors = TRUE,open_new_window=open_new_plot,win=spatstat::owin(c(-3,3),c(-3,3)))+add_title("True Poisson intensity surface along with the point pattern, W=[-3,3]x[-3,3]",lambda =trueintsurfmix4$lambda,m=trueintsurfmix4$m,n=ppmix4$n)
plotmix_2d(trueintsurfmix4,ppmix4,colors = TRUE,open_new_window=open_new_plot,win=spatstat::owin(c(-4,4),c(-4,4)))+add_title("True Poisson intensity surface along with the point pattern, W=[-4,4]x[-4,4]",lambda =trueintsurfmix4$lambda,m=trueintsurfmix4$m,n=ppmix4$n)
bigwin=spatstat::owin(c(-4,4),c(-4,4))

plotstring("\nWe run the selectMix function to obtain\nmodel selection criteria. Note the output\nwith the AIC, BIC and ICLC values\nin the console. We compare mixture\nmodels with 1 to 5 components, and the\ncriteria are calculated by fitting the models\nwith the specific components.")

plotstring("\nCalculations can take a while. Note that the\nsmaller a criterion value is, the better the\nmodel fits. If all goes well, we should be\ngetting the model with 4 components as\nthe most appropriate, unless two\ncomponents are on top of each other, in which\ncase the number of components proposed\nwill be much less than the true number.")

ModelSel2=selectMix(ppmix4,1:5,truncate=FALSE)

plotstring("\nNow we fit several mixture models to\nestimate the intensity surface of the point\npattern, and then run the gof test for\neach one. We should be getting a good fit for\nm=4 components.")

plotstring("\nFirst fit 2 components and perform the gof.\nThen fit 3, 4 (the true number of components),\nand then 5. Each time we estimate\nthe surface using the surface of posterior\nmeans from a DAMCMC fit (permuted\nusing the best permutation under SEL).")

plotstring("\nCalculations for 2 components are performed.\nSee the console for the results.")

DAMCMCfit2=est_mix_damcmc(pp = ppmix4, m = 2,L=10000)
post_fixedDAMCMCfitIC2 = FixLS_da(DAMCMCfit2,approx = FALSE)
permSurfaceofPostMeansIC2=GetPMEst(post_fixedDAMCMCfitIC2)
print(mc_gof(ppmix4, permSurfaceofPostMeansIC2, 0.05))

plotstring("\nCalculations for 3 components are performed.\nSee the console for the results.")

DAMCMCfit3=est_mix_damcmc(pp = ppmix4, m = 3,L=10000)
post_fixedDAMCMCfitIC3 = FixLS_da(DAMCMCfit3,approx = FALSE)
permSurfaceofPostMeansIC3=GetPMEst(post_fixedDAMCMCfitIC3)
print(mc_gof(ppmix4, permSurfaceofPostMeansIC3, 0.05))

plotstring("\nCalculations for 4 components are performed.\nSee the console for the results.")

DAMCMCfit4=est_mix_damcmc(pp = ppmix4, m = 4,L=10000)
post_fixedDAMCMCfitIC4 = FixLS_da(DAMCMCfit4,approx = FALSE)
permSurfaceofPostMeansIC4=GetPMEst(post_fixedDAMCMCfitIC4)
print(mc_gof(ppmix4, permSurfaceofPostMeansIC4, 0.05))

plotstring("\nCalculations for 5 components are performed.\nSee the console for the results.")

DAMCMCfit5=est_mix_damcmc(pp = ppmix4, m = 5,L=10000)
post_fixedDAMCMCfitIC5 = FixLS_da(DAMCMCfit5,approx = FALSE)
permSurfaceofPostMeansIC5=GetPMEst(post_fixedDAMCMCfitIC5)
print(mc_gof(ppmix4, permSurfaceofPostMeansIC5, 0.05))

plotstring("Now consider the K-S test. We run\nfunction kstest2d and compare the pattern\nunder consideration against patterns from\nmodels with 2-5 components generated from\nthe models we just fit.")
pp2 <- rsppmix(intsurf = permSurfaceofPostMeansIC2,truncate = FALSE)# draw points
pp3 <- rsppmix(intsurf = permSurfaceofPostMeansIC3,truncate = FALSE)# draw points
pp4 <- rsppmix(intsurf = permSurfaceofPostMeansIC4,truncate = FALSE)# draw points
pp5 <- rsppmix(intsurf = permSurfaceofPostMeansIC5,truncate = FALSE)# draw points

plotstring("\nCalculations for 2 components are performed.\nSee the console for the results.")
kstest2d(pp2, ppmix4)
plotstring("\nCalculations for 3 components are performed.\nSee the console for the results.")
kstest2d(pp3, ppmix4)
plotstring("\nCalculations for 4 components are performed.\nSee the console for the results.")
kstest2d(pp4, ppmix4)
plotstring("\nCalculations for 5 components are performed.\nSee the console for the results.")
kstest2d(pp5, ppmix4)

plotstring("Now use function kstest2dsurf to\nperform the K-S test against a proposed\nsurface and report the results in the console.")
plotstring("\nCalculations for 2 components are performed.\nSee the console for the results.")
kstest2dsurf(ppmix4, permSurfaceofPostMeansIC2)
plotstring("\nCalculations for 3 components are performed.\nSee the console for the results.")
kstest2dsurf(ppmix4, permSurfaceofPostMeansIC3)
plotstring("\nCalculations for 4 components are performed.\nSee the console for the results.")
kstest2dsurf(ppmix4, permSurfaceofPostMeansIC4)
plotstring("\nCalculations for 5 components are performed.\nSee the console for the results.")
kstest2dsurf(ppmix4, permSurfaceofPostMeansIC5)

plotstring("See the vignettes and help pages of\nthese functions for more details.\nNote that this analysis can be\nrepeated with truncate=TRUE in every function\ncalled, in order to handle edge effects.\nThanks")

suppressWarnings( par(oldpar))
options(oldoptions)
