#Plots_off()
oldoptions <- options(warn=-1)
oldpar <- par()
#devAskNewPage(ask = TRUE)
cat("NOTE: if you're running rstudio make sure\nthe plotting window is large otherwise you will\nget a margin error.")

open_new_plot=FALSE
#openwin_sppmix(TRUE)

#create a normal mixture
#object "normmix"
plotstring("\nIn this demo we simulate and fit a Marked\nInhomogeneous Poisson Point Process (MIPPP)\nwhere we assume that the events consist\nof a location (ground process) component\nand a mark component, with the mark\ndistribution depending on locations.\nThis leads to mark distributions described\nby random fields over the window of observation.")

plotstring("We start with some simulations.\nWe use function 'rMIPPP_cond_loc'\nto simulate both the discrete and continuous\nmark case. Currently, Bayesian estimation\nof the model parameters is implemented\nin function 'est_MIPPP_cond_loc' ONLY\nfor the discrete mark case.")

plotstring("If the marks are discrete valued, we use\na Markov Random Field (e.g., Ising for two marks)\nand generate the mark values at the\nobserved points as well as, produce\nthe probability fields over the window.\nIn the continuous mark case, we generate\nGaussian random fields (or their transformations)\nin order to build the mark probability fields.")

newMPP=rMIPPP_cond_loc(
  open_new_window=open_new_plot,
  gammas=c(.1,.2,.3),r=.5,
  bigwin=spatstat.geom::owin(c(-5,5),c(-5,5)))

plotstring("The parameters 'gammas' act as weights\nin the discrete mark distribution.\nLarge values lead to small probabilities\nof observing the corresponding mark value.\nFunction 'rMIPPP_cond_loc' returns both the\nground process and the mark distribution,\nas well as, plots the probability\nfields of observing a mark.")

plotstring("The neighborhood system is controlled by a\nparameter 'r', with large values leading to\nsmoother fields (since we borrow strength\nacross space) whereas small values tend to\nconcentrate the mark masses about\nthe observed locations.")

summary(newMPP)

true_gammas=newMPP$gammas
genMPP=newMPP$genMPP
genMPP
newMPP$r

plotstring("Function summary gives us information\non the return value of the generated\nMIPPP from function 'rMIPPP_cond_loc'.")

plotstring("We fit the Bayesian model for the\nMIPPP and obtain posterior realizations\nof the model parameters, including\nthe 'gammas' and the component mixture\nparameters for the ground process.")

mpp_est <- est_MIPPP_cond_loc(genMPP,newMPP$r,
  open_new_window=open_new_plot,
  hyper=0.2,L=20000)

plotstring("The result of a call to 'est_MIPPP_cond_loc'\nis an 'MIPPP_fit' object to which we can\napply the generic functions 'summary'\nand 'plot'. We also use function 'GetStats'\nin order to perform a posterior analysis\non the 'gammas' parameters.")

plot(mpp_est,open_new_window=open_new_plot)

GetStats(mpp_est$gen_gammas[,1])$Mean
GetStats(mpp_est$gen_gammas[,2])$Mean
GetStats(mpp_est$gen_gammas[,3])$Mean

GetStats(mpp_est$gen_gammas[,1])$CredibleSet
GetStats(mpp_est$gen_gammas[,2])$CredibleSet
GetStats(mpp_est$gen_gammas[,3])$CredibleSet

plotstring("\nNow we generate a MIPPP with continuous\nmarks and mark distribution defined\nconditionally on location (random fields).\nWe call the function 'rMIPPP_cond_loc'\nwith parameter 'discrete_mark=FALSE'.")

plotstring("\nThe mark probability fields implemented\ninclude Gaussian random fields for real marks\nand Chi-square random fields for positive\nreal valued marks. The spatial covariance\nfunction for the GRFs is modeled using\nthe Matern family of covariance functions.\nNote that the GRFs are stationary (mean\ndoes not depend on location) and\nisotropic (no direction is preferred).")

plotstring("\nReal marks first (Gaussian random field)...")
newMPPcontGRF=rMIPPP_cond_loc(discrete_mark = FALSE,
  open_new_window=open_new_plot,
  bigwin=spatstat.geom::owin(c(-5,5),c(-5,5)))

plotstring("\nPositive real marks next (Chi-square random field)...")
newMPPcontChiSQ=rMIPPP_cond_loc(mark_distr_choice=1,
  open_new_window=open_new_plot,
  discrete_mark = FALSE,
  bigwin=spatstat.geom::owin(c(-5,5),c(-5,5)))

plotstring("\nTurn to an application. We visualize\naggregate income levels in MO by county\nusing data from the American Community\nSurvey (ACS). We plot in the original\nscale first. We use function 'PlotUSAStates'\nto produce these plots.")

plotstring("\nWe pass the marks vector\n'MOAggIncomeLevelsPerCounty' to the function\n'PlotUSAStates', which contains the\naggregate income values of Missourian\ncounties.")

ret=PlotUSAStates(showcounties=TRUE,
                  open_new_window = open_new_plot,
                  states=c('Missouri'),showcentroids=TRUE,
                  typecentroid=1,discretelevels=FALSE,shownames=TRUE,
                  plotlevels=TRUE,marks=MOAggIncomeLevelsPerCounty,
                  main="Aggregate Income in MO, 2014",
                  guidemain = "Income level",
                  namescolor="gray",boundarycolor="gray")

plotstring("\nNow we plot in the log scale.")
ret=PlotUSAStates(showcounties=TRUE,
                  open_new_window = open_new_plot,
                  states=c('Missouri'),showcentroids=TRUE,
                  typecentroid=1,discretelevels=FALSE,shownames=TRUE,
                  plotlevels=TRUE,marks=log(MOAggIncomeLevelsPerCounty),
                  main="Aggregate Income in MO, 2014",
                  guidemain = "Income level\n(log scale)",
                  namescolor="gray",boundarycolor="gray")

plotstring("\nThe 'marker points' for a random set\nare essentially the most south-western points\n(for convex sets, e.g., the county boundary\ndescribes a random set called the grain\nand the union of all the grains\nis a random set known as the Boolean model).")

plotstring("\nMarker points for a Boolean model\nare used to summarize the random set\nvia the created pattern which has\nnice properties. next we plot the marker\npoints, county boundaries and names.")
ret=PlotUSAStates(showcounties=TRUE,
                  open_new_window = open_new_plot,
                  states=c('Missouri'),showcentroids=TRUE,typecentroid = 1,
                  discretelevels=FALSE,shownames=TRUE,
                  plotlevels=FALSE,marks=log(MOAggIncomeLevelsPerCounty),
                  main="Marker points for Missouri counties")

plotstring("\nNow plot only the marker points\nand treat this as a marked IPPP.")
ret=PlotUSAStates(showcounties=TRUE,
                  open_new_window = open_new_plot,
                  states=c('Missouri'),showcentroids=TRUE,typecentroid = 1,
                  discretelevels=FALSE,shownames=FALSE,
                  plotlevels=FALSE,marks=log(MOAggIncomeLevelsPerCounty),
                  main="Marker points for Missouri counties",boundarycolor = NULL)

plotstring("Let us discretize log(income) to\n3 levels; low if <=20 (mark value 1),\naverage if >20 and <=23 (mark value 2),\nand high if >23 (mark value 3).")
newmarks=rep(0,length(MOAggIncomeLevelsPerCounty))

newmarks[log(MOAggIncomeLevelsPerCounty)<=20]=1

newmarks[log(MOAggIncomeLevelsPerCounty)>20
         & log(MOAggIncomeLevelsPerCounty)<=23]=2

newmarks[log(MOAggIncomeLevelsPerCounty)>23]=3

table(newmarks)

levels=c("low","average","high")

plotstring("We plot the aggregate income per\ncounty in the new discrete scale next.\nNote that: 1=low, 2=average, 3=high.")
ret=PlotUSAStates(showcounties=TRUE,
                  open_new_window = open_new_plot,
                  states=c('Missouri'),showcentroids=TRUE,
                  typecentroid=1,discretelevels=TRUE,
                  shownames=TRUE,plotlevels=TRUE,
                  main="Aggregate Income in MO, 2014",
                  marks=newmarks,levels=levels,
                  guidemain = "Income level",
                  namescolor="gray",boundarycolor="gray")

plotstring("Now fit a marked IPPP model on the\npattern of marked marker points.")
MPP=ret$PPPMarker

mpp_est <- est_MIPPP_cond_loc(MPP,r=1,
                              open_new_window=open_new_plot,
                              hyper=0.2,L=20000)

plotstring("The function 'plot_MPP_probs' summarizes\nthe mark probabilities of a marked point\npattern, by displaying for each location,\nthe probabilities of observing each of\nthe discrete marks at that location.\nThis function works for discrete marks only.")
plot_MPP_probs(mpp_est, open_new_window = open_new_plot)

plotstring("Note that this analysis can be\nrepeated with truncate=TRUE in every function\ncalled, in order to handle edge effects.\nSee the model vignettes for more details.\nThanks")

suppressWarnings( par(oldpar))
options(oldoptions)
#on.exit(devAskNewPage(ask = FALSE))
#par(ask=FALSE)
