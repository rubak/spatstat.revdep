#Plots_off()
oldoptions <- options(warn=-1)
#oldpar <- par(ask = interactive() && dev.interactive(orNone=TRUE))
oldpar <- par()

cat("NOTE: if you're running rstudio make sure\nthe plotting window is large otherwise you will\nget a margin error.")

open_new_plot=FALSE
#openwin_sppmix(TRUE)

#create a normal mixture
#object "normmix"
plotstring("First build a normmix object.\nThis is simply a mixture with\nnormal components.")
truemix <- normmix(ps = c(.3, .7),
                mus = list(c(0.2, 0.2), c(.8, .8)),
                sigmas = list(.01*diag(2), .01*diag(2)))

truemix

#use summary to get extra info on
#the normmix object
plotstring("Use summary to see details\non the normmix object and\na simple plot() to display it.\nWe can also change the observation\nwindow and plot the contours.")
summary(truemix)

plot(truemix,whichplots=0,open_new_window=open_new_plot)

plot(truemix,whichplots=0,xlim=c(0,1),ylim=c(0,1),open_new_window=open_new_plot)

plot(truemix,whichplots=0,xlim=c(-1,2),ylim=c(-1,2),open_new_window=open_new_plot)

plot(truemix,whichplots=0,xlim=c(0,1),ylim=c(0,1),contour = TRUE,open_new_window=open_new_plot)

plotstring("We can generate the parameters ps,\nmus and sigmas of the normal mixture\nover the window [-3,3]x[-3,3] and \nplot the density")

truemix3=rnormmix(m = 3, sig0 = .1, df = 5,xlim= c(-3,3), ylim = c(-3,3))
plot(truemix3,whichplots=0,xlim=c(-3,3),ylim=c(-3,3),open_new_window=open_new_plot)

truemix5=rnormmix(m = 5, sig0 = .1, df = 5,xlim= c(-3,3), ylim = c(-3,3))
plot(truemix5,whichplots=0,xlim=c(-3,3),ylim=c(-3,3),open_new_window=open_new_plot)

plotstring("Build an intensity_surface object\nbased on the last normmix object by\nsetting the parameter lambda=100,\nwhich denotes the average number\nof points over the window\n[-3,3]x[-3,3].")
intsurf=to_int_surf(truemix5,lambda = 100,
                    win =spatstat::owin( c(-3,3),c(-3,3)))

#plot(intsurf)

plotmix_2d(intsurf,open_new_window=open_new_plot)

plotmix_2d(intsurf,contour = TRUE,open_new_window=open_new_plot)+add_title("Contour plot of the intensity surface",lambda =intsurf$lambda,m=intsurf$m)
#using normmix we can create an
#intensity_surface object directly
plotstring("Now build an intensity_surface object.\nParameter lambda=200, denotes the\naverage number of points over the\nwindow [-1,1]x[-2,3]")
demo_truemix3comp <- normmix(ps=c(.2, .5,.3),
  mus=list(c(-0.3, -1.3), c(.1,.5),c(0.7, 1.7)),
  sigmas = list(.3*diag(2),.5*diag(2),
   .2*diag(2)), lambda = 200,
          win = spatstat::owin(c(-1,1),c(-2,3)))
demo_truemix3comp

#use summary to get extra info on
#the surface object
summary(demo_truemix3comp)
plotstring("Use summary to see details\non the intensity_surface object\nand a plot() call to plot it.")

plotmix_2d(demo_truemix3comp,open_new_window=open_new_plot)
plotmix_2d(demo_truemix3comp,contour = TRUE,open_new_window=open_new_plot)

plotstring("Generate a point pattern from a\nPoisson with the given intensity surface.\nWe can also change the lambda\nand the window on the fly and\nplot the intensity surface.")
pp1 <- rsppmix(intsurf = intsurf)# draw points
plot(pp1, mus = intsurf$mus,open_new_window=open_new_plot)#plot the mixture means as well
plotmix_2d(intsurf, pp1,colors = TRUE, open_new_window=open_new_plot)
plotmix_2d(intsurf, pp1,colors = TRUE,contour = TRUE,open_new_window=open_new_plot)+add_title("Contour plot of the intensity surface",lambda =intsurf$lambda,m=intsurf$m,n=pp1$n )
plotmix_2d(intsurf, pp1,colors = TRUE,win=spatstat::owin(c(-1,2),c(-1,2)),open_new_window=open_new_plot)
plotmix_2d(intsurf, pp1,colors = TRUE,win=spatstat::owin(c(-1,2),c(-1,2)),contour = TRUE,open_new_window=open_new_plot)+add_title("Contour plot of the intensity surface",lambda =intsurf$lambda,m=intsurf$m,n=pp1$n )

plotstring("Finally, we produce 3d plots\nof the intensity surfaces we saw.")

plot(demo_truemix3comp,main="True normal mixture intensity surface with 3 components")

plot(truemix5,whichplots=1,xlim=c(-3,3),ylim=c(-3,3),title1="True normal mixture density with 5 components")

plot(intsurf,main="True normal mixture intensity surface with 5 components")

plotstring("More details can be found in\nthe vignettes and help pages of these\nobjects and their related functions.\nThanks")

suppressWarnings( par(oldpar))
options(oldoptions)
#par(ask=FALSE)
