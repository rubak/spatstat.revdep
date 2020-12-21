FixLS_daRetBestPerm<- function(fit, burnin = floor(fit$L / 10),
                    xlab = "x",ylab = "y",
                    approx = TRUE,
                    plot_result = FALSE,
                    run_silent = FALSE)
{
  if(burnin>=fit$L)
  {
    if(!run_silent)
      cat("\nBad burnin value. Using the default.")
    burnin = floor(fit$L/10)
  }
  fit <- drop_realization(fit, burnin)
  win <- spatstat.geom::domain(fit$data)
  if(class(fit)=="bdmcmc_res")
  {
    tab=GetBDTable(fit,F)
    fit=GetBDCompfit(fit,tab$MAPcomp)$BDgens
    burnin = fit$L / 10
  }
  m <- fit$m
  if(m>5 && approx ==FALSE && !run_silent)
    cat("\nWARNING: m>5, this will take a long long time to complete...\n")
  xlims1 <- c(win$xrange)
  ylims1 <- c(win$yrange)
  L <- dim(fit$genps)[1]
  BestPermutation=1
  if (approx == TRUE)
  {
    if(!run_silent)
      cat("\nStarting the approximate relabeling algorithm.\n")
    iter=1
    lagnum=lagstart=min(.05*(fit$L-burnin),50)+iter
    while(check_labels(fit,lagnum=lagnum,showmessages=FALSE))
    {
      if(!run_silent)
        cat("\nIteration",iter,", lag=",lagnum,"\n")
      iter=iter+1
      lagnum=lagstart+iter*10
      permgens <- PostGenGetBestPermIdenConstraint_sppmix(fit)
      fit$allgens_List = permgens$allgens_List
      fit$genps = permgens$genps
      fit$genmus = permgens$genmus
      fit$gensigmas = permgens$gensigmas
      fit$genzs = permgens$genzs
      fit$ApproxCompMass=PermuteZs_sppmix(fit$ApproxCompMass ,permgens$best_perm)
      #      fit$genlamdas = permgens$genlamdas
      BestPermutation=permgens$best_perm
    }
  } else
  {
    if(!run_silent)
      cat("\nStarting the decision theoretic relabeling algorithm.\n")
    iter=1
    lagnum=lagstart=min(.05*(fit$L-burnin),50)+iter
    while(check_labels(fit,lagnum=lagnum,showmessages=FALSE))
    {
      if(!run_silent)
        cat("\nIteration",iter,", lag=",lagnum,"\n")
      iter=iter+1
      lagnum=lagstart+iter*10
      permgens <- PostGenGetBestPerm_sppmix(fit$allgens_List)
      fit$allgens_List = permgens$permuted_gens
      fit$genps = permgens$permuted_ps
      fit$genmus = permgens$permuted_mus
      fit$gensigmas = permgens$permuted_sigmas
      fit$genzs=PermuteZs_sppmix(fit$genzs ,permgens$best_perm)
      fit$ApproxCompMass=PermuteZs_sppmix(fit$ApproxCompMass ,permgens$best_perm)
      #      fit$genlamdas = fit$genlamdas
      BestPermutation=permgens$best_perm
    }
  }
  post_ps <- colMeans(fit$genps)
  mus <- apply(fit$genmus, 1:2, mean)
  mean_mat <- function(mats) Reduce("+", mats) / length(mats)
  sigmas <- apply(fit$gensigmas, 2, mean_mat)
  mean_lambda <- mean(fit$genlamdas)
  post_mus <- post_sigmas <- vector("list", m)
  for (i in 1:m) {
    post_mus[[i]] <- mus[i, ]
    post_sigmas[[i]] <- matrix(sigmas[, i], 2, 2)
  }
  post_normix = normmix(post_ps, post_mus, post_sigmas, lambda = mean_lambda,
                        win = win)
  if (plot_result == TRUE)
  {
    print(plot.sppmix(fit$data, post_normix$mus))
    print(plotmix_2d(post_normix,fit$data))
    if (approx == TRUE)
      plot.intensity_surface(post_normix,
                             main = "Surface of posterior means (Identifiability Constraint)")
    else
      plot.intensity_surface(post_normix,
                             main = "Surface of posterior means (Permuted Labels)")
  }
  return(list(fit=fit,BestPermutation=BestPermutation))
}

testit<-function(truncate=FALSE,open_new_plot=FALSE)
{
  ContinentalUSA_state_names=c(
    "California","Florida","Georgia","Idaho"
    ,"Illinois","Iowa","Kentucky","Louisiana"
    ,"Maryland","Michigan","Minnesota","Missouri"
    ,"New York"
    ,"Oregon"
    ,"Tennessee"
    ,"Texas"
    ,"Virginia"
    ,"Wisconsin"
    ,"Arizona"
    ,"Arkansas"
    ,"Colorado"
    ,"Indiana"
    ,"Connecticut"
    ,"Nebraska"
    ,"New Mexico"
    ,"North Carolina"
    ,"Ohio"
    ,"Maine"
    ,"Massachusetts"
    ,"Mississippi"
    ,"Montana"
    ,"Oklahoma"
    ,"South Carolina"
    ,"South Dakota"
    ,"Utah"
    ,"Washington"
    ,"West Virginia"
    ,"Wyoming"
    ,"Delaware"
    ,"Rhode Island"
    ,"Alabama"
    ,"North Dakota"
    ,"Pennsylvania"
    ,"Vermont"
    ,"Kansas"
    ,"Nevada"
    ,"New Hampshire"
    ,"New Jersey"  )
  truemix <- rnormmix(m = 5, sig0 = .1, df = 5,xlim=c(-2,2),ylim=c(-2,2))
  trueintsurf=to_int_surf(truemix,lambda = 100, win = spatstat.geom::owin(c(-2,2),c(-2,2)))
  plot(trueintsurf,main = "True Poisson intensity surface (mixture of normal components)")
  genPPP <- rsppmix(trueintsurf,truncate=FALSE)
  ModelSel1=selectMix(genPPP,1:7,runallperms=0,truncate=FALSE)

  if(0)
{
  truemix4=rnormmix(m = 4, sig0 = .1, df = 5,xlim= c(-2,2), ylim = c(-2,2))
#  plot(truemix4,xlim= c(-2,2), ylim = c(-2,2),whichplots=0,open_new_window=open_new_plot)+add_title("True mixture of normals density")
  trueintsurfmix4=to_int_surf(truemix4,lambda = 150,win =spatstat.geom::owin( c(-2,2),c(-2,2)))
  bigwin=spatstat.geom::owin(c(-4,4),c(-4,4))
  ppmix4 <- rsppmix(intsurf = trueintsurfmix4,truncate = truncate,win=bigwin)# draw points
  print(plotmix_2d(trueintsurfmix4,ppmix4,open_new_window=open_new_plot,win=spatstat.geom::owin(c(-4,4),c(-4,4)))+add_title("True Poisson intensity surface along with the point pattern, W=[-4,4]x[-4,4]",lambda =trueintsurfmix4$lambda,m=trueintsurfmix4$m,n=ppmix4$n))
  BDMCMCfit=est_mix_bdmcmc(pp = ppmix4, m = 10,
                           L=10000,truncate = truncate)
  BDtab=GetBDTable(BDMCMCfit)#retrieve frequency table and MAP estimate for number of components
  MAPm=BDtab$MAPcomp
  BDMCMCfitMAPcomp=GetBDCompfit(BDMCMCfit,MAPm)
  print(range(BDMCMCfitMAPcomp$BDgens$genmus))

  BDMCMCfit=drop_realization(BDMCMCfit)

  BDMCMCfit=drop_realization(BDMCMCfit,
        (BDMCMCfit$Badgen==1))

  plot_CompDist(BDMCMCfit,open_new_window=open_new_plot)
  BDtab=GetBDTable(BDMCMCfit)#retrieve frequency table and MAP estimate for number of components
  MAPm=BDtab$MAPcomp
  BDMCMCfitMAPcomp=GetBDCompfit(BDMCMCfit,MAPm)
#  BDMCMCfitMAPcomp
  BDMCMCfitMAPcompgens<-BDMCMCfitMAPcomp$BDgens
  print(range(BDMCMCfitMAPcompgens$genmus[,1,]))
  print(range(BDMCMCfitMAPcompgens$genmus[,2,]))

  Get_User_Input_sppmix("Continue?")

#  plotmix_2d(BDMCMCfitMAPcomp$BDsurf,ppmix4,open_new_window=open_new_plot,win=bigwin)+add_title("MAP Poisson intensity surface along with the point pattern",lambda =BDMCMCfitMAPcomp$BDsurf$lambda,m=BDMCMCfitMAPcomp$BDsurf$m,n=ppmix4$n,L=BDMCMCfitMAPcomp$BDgens$L)

#  plot_chains(BDMCMCfitMAPcompgens,open_new_window=open_new_plot,separate = FALSE)

  labelswitch=check_labels(BDMCMCfitMAPcompgens,burnin=0)

  post_fixedBDMCMCfitIC = FixLS_da(BDMCMCfitMAPcompgens,burnin=0)
#  plot_chains(post_fixedBDMCMCfitIC,open_new_window=open_new_plot,separate = FALSE)
#  print(plot_ind(post_fixedBDMCMCfitIC,burnin=0,open_new_window=open_new_plot)+add_title("Posterior means of the membership indicators (IC permuted labels)", m = post_fixedBDMCMCfitIC$m, n = post_fixedBDMCMCfitIC$data$n))
  permSurfaceofPostMeansIC=GetPMEst(post_fixedBDMCMCfitIC,burnin=0)
  print(plotmix_2d(permSurfaceofPostMeansIC,ppmix4,open_new_window=open_new_plot,win=bigwin)+add_title("Poisson surface of posterior means (IC)",lambda =permSurfaceofPostMeansIC$lambda,m=permSurfaceofPostMeansIC$m,n=ppmix4$n,L=post_fixedBDMCMCfitIC$L))

  post_fixedBDMCMCfitSEL = FixLS_da(BDMCMCfitMAPcompgens,approx=FALSE,burnin=0)
#  plot_chains(post_fixedBDMCMCfitSEL,open_new_window=open_new_plot,separate = FALSE)
#  print(plot_ind(post_fixedBDMCMCfitSEL,burnin=0,open_new_window=open_new_plot)+add_title("Posterior means of the membership indicators (best permutation)", m = post_fixedBDMCMCfitSEL$m, n = post_fixedBDMCMCfitSEL$data$n))
  permSurfaceofPostMeansSEL=GetPMEst(post_fixedBDMCMCfitSEL,burnin=0)
  print(plotmix_2d(permSurfaceofPostMeansSEL,ppmix4,open_new_window=open_new_plot,win=bigwin)+add_title("Poisson surface of posterior means (best permutation)",lambda =permSurfaceofPostMeansSEL$lambda,m=permSurfaceofPostMeansSEL$m,n=ppmix4$n,L=post_fixedBDMCMCfitSEL$L))
}

}

Get_User_Input_sppmix<- function(prompt_string="",modeYN=1)
{
  options(warn=-1)
  #modeYN=1, ask for yes or no
  #modeYN=0, ask for a value
  if(modeYN)
  {
    ret=0;#no, 1 is yes
    while(1)
    {
      ANSWER <- readline(paste(prompt_string,
                               "(Y)es or (N)o? or (Q)uit "))
      if (substr(ANSWER, 1, 1) == "n"
          || substr(ANSWER, 1, 1) == "N")
      {
        ret=0
        break
      }
      if (substr(ANSWER, 1, 1) == "y"
          || substr(ANSWER, 1, 1) == "Y")
      {
        ret=1
        break
      }
      if (substr(ANSWER, 1, 1) == "q"
          || substr(ANSWER, 1, 1) == "Q")
      {
        stop("Execution ended by the user")
      }
    }
  }
  else
  {
    ret=0;#default return value
    while(1)
    {
      #message("Enter ",prompt_string,":")
      val <- readline(paste("Enter",prompt_string,": "))
      #scan(what=double())
      #      check=is.na(as.numeric(val))
      if(is.na(as.numeric(val)))
        message("Enter a number, not letters")
      else
      {
        ret=as.numeric(val)
        break
      }
    }
  }
  options(warn=0)
  return(ret)
}

CreateDemoSurfs<- function()
{
#  R -e "rmarkdown::render('knitr_example.Rmd')"
  Plots_off()
  demo_truemix3comp <- normmix(ps=c(.2, .5,.3),
   mus=list(c(-0.3, -1.3), c(.1,.5),
            c(0.7, 1.7)),
   sigmas = list(.3*diag(2),.5*diag(2),
                 .2*diag(2)))
  demo_intsurf3comp=to_int_surf(
    demo_truemix3comp,lambda = 200,
    win = spatstat.geom::owin(c(-1,1),c(-2,3)))
  #generate a point pattern
  PP<-rsppmix(demo_intsurf3comp)
  #normmix plots
  plot(demo_truemix3comp,
       xlim=demo_intsurf3comp$window$xrange,
       ylim=demo_intsurf3comp$window$yrange)

  # postfit=est_mix_damcmc(pp1,m=2,L=1200)
  #
  # avgsurf=plot_avgsurf(postfit, LL = 50, burnin = 1000)
  #
  # plotmix_3d(avgsurf)
  #
  # plot_density(avgsurf)
  #intensity_surface plots
  # #2d plots
  # plotmix_2d(demo_intsurf)#elevation plot

  # plotmix_2d(demo_intsurf,contour=TRUE)#contour plot
  #
  # #get a point pattern
  # pp1 <- rsppmix(demo_intsurf)
  #
  # #plots with the point pattern
  #
  # plot(pp1)
  # plot(pp1,demo_intsurf$mus)
  #
  # plot2dPP(pp1)
  # plot2dPP(pp1,demo_intsurf$mus)
  #
  plot(demo_intsurf3comp,main = "True Poisson intensity surface (mixture of normal components)")


  plot(demo_intsurf3comp)
  plotmix_2d(demo_intsurf3comp)
#  devtools::use_data(demo_truemix3comp,overwrite = TRUE)
#  devtools::use_data(demo_intsurf3comp,overwrite = TRUE)
}

GetMixtureLimitsList<- function(mix)
{
  #easy way to determine limits
  #go through each mean x-y coord and
  #find the smallest mu-5sig and max mu+5sig
  xlim=c(10000000000,-10000000000)
  ylim=c(10000000000,-10000000000)
  for(j in 1:mix$m)
  {
    cxmin=as.vector(mix$mus[[j]])[1]
    -5*sqrt(as.vector(mix$sigmas[[j]])[1])
    cxmax=as.vector(mix$mus[[j]])[1]
    +5*sqrt(as.vector(mix$sigmas[[j]])[1])
    cymin=as.vector(mix$mus[[j]])[2]
    -5*sqrt(as.vector(mix$sigmas[[j]])[4])
    cymax=as.vector(mix$mus[[j]])[2]
    +5*sqrt(as.vector(mix$sigmas[[j]])[4])
    if(cxmin<xlim[1])
      xlim[1]=cxmin
    if(cxmax>xlim[2])
      xlim[2]=cxmax
    if(cymin<ylim[1])
      ylim[1]=cymin
    if(cymax>ylim[2])
      ylim[2]=cymax
  }
  return(list(xlim,ylim))
}

#' Opens a new graphics window
#'
#' @description
#' This function is independent of the OS present, and is useful
#' when working outside of RStudio. The latter GUI
#' places all plots under the plots tab, but if working
#' in the R GUI the plots will be overwritten if you don't open a new device.
#'
#' For examples see
#'
#' \url{http://faculty.missouri.edu/~micheasa/sppmix/sppmix_all_examples.html
#' #openwin_sppmix}
#'
#' @param check2open Logical: TRUE to open a newplot, FALSE do not open.
#' @details This function is used by almost all plotting functions of the \code{\link{sppmix}} package.
#' @author Sakis Micheas
#' @examples
#' \donttest{
#' openwin_sppmix(TRUE)}
#'
#' @export
openwin_sppmix <- function(check2open=FALSE){
  if(check2open)
  {
    dev.new()
#    if(.Platform$OS.type == "windows") {
#      windows(width=400, height=300,xpos=10, ypos=10)
#    }
#    if(.Platform$OS.type == "unix") {
#      X11(width = 8, height = 6, xpos = 0, ypos = 0)
#    }
  }
}

MakeNormMixFromMixtureList<- function(mix)
{
  #takes a mixture list and returns
  #a normmix object
  m=length(mix);
  ps=vector("double", m);
#  cat(m,"\n")
  mus=vector("list", m);
  sigmas=vector("list", m);
  for(i in 1:m)
  {
    ps[[i]]=as.numeric(mix[[i]]$p)
#    cat(mix[[i]]$p,"\n")
    mus[[i]]=as.vector(mix[[i]]$mu)
    sigmas[[i]]=as.matrix(mix[[i]]$sigma)
  }
  #  cat(sum(ps))
  norm_mix=normmix(ps, mus, sigmas,estimated=TRUE);
  return (norm_mix)
}

MakeMixtureListFromNormMix<- function(norm_mix)
{
  #takes a mixture list and returns
  #a normmix object
  stopifnot(is.normmix(norm_mix))
  m=length(norm_mix$ps);
  mixlist=vector("list", m)
  for(i in 1:m)
  {
    compi=list(p=norm_mix$ps[i],
      mu=norm_mix$mus[[i]],
      sigma=norm_mix$sigmas[[i]])
    mixlist[[i]]=compi
  }
  return (mixlist)
}

#' Run an sppmix package demo or vignette
#'
#' @description
#' The function starts the different tutorials of
#' the \code{sppmix} package.
#'
#' For examples see
#'
#' \url{http://faculty.missouri.edu/~micheasa/sppmix/sppmix_all_examples.html
#' #Demo_sppmix}
#'
#' @param whichdemo A number indicating which
#' demo to run. Do not pass any number in order
#' to see all the availabe choices.
#'
#' @author Sakis Micheas
#' @examples
#' \donttest{
#' Demo_sppmix()#shows all available demos and opens the vignettes page
#' Demo_sppmix(1)#demo on sppmix objects}
#'
#' @export
Demo_sppmix<-function(whichdemo=NULL)
{
  if(is.null(whichdemo))
  {
    cat("\nAvailable demos (choose a number):\n")
    cat("1) illustrates sppmix objects\n")
    cat("2) plotting in the sppmix package\n")
    cat("3) IPPP model fitting: fixed number of components and no edge effects present.\n")
    cat("4) IPPP model fitting: fixed number of components with edge effects present.\n")
    cat("5) IPPP model fitting: random number of components and no edge effects present.\n")
    cat("6) IPPP model fitting: random number of components with edge effects present.\n")
    cat("7) model selection and checking in sppmix.\n")
    cat("8) Marked IPPP model fitting.\n")
    browseVignettes(package = 'sppmix')
  }
  else
  {
    if(whichdemo==1)# || whichdemo=="objects")
    {
      print(demo(topic = 'use_sppmix_objects', package = 'sppmix'))
    } else
    if(whichdemo==2)# || whichdemo=="plots")
    {
      print(demo(topic = 'use_sppmix_plots', package = 'sppmix'))
    } else
    if(whichdemo==3)# || whichdemo=="modelfits")
    {
      print(demo(topic = 'use_sppmix_DAMCMC_noedgeeffects', package = 'sppmix'))
    } else
    if(whichdemo==4)# || whichdemo=="modelcheck")
    {
      print(demo(topic = 'use_sppmix_DAMCMC_edgeeffects', package = 'sppmix'))
    } else
    if(whichdemo==5)# || whichdemo=="modelcheck")
    {
      print(demo(topic = 'use_sppmix_BDMCMC_noedgeeffects', package = 'sppmix'))
    } else
    if(whichdemo==6)# || whichdemo=="modelcheck")
    {
#      cat("Demo 6 is not currently available")
      print(demo(topic = 'use_sppmix_BDMCMC_edgeeffects', package = 'sppmix'))
    } else
    if(whichdemo==7)# || whichdemo=="modelcheck")
    {
      print(demo(topic = 'use_sppmix_modelcheck', package = 'sppmix'))
    } else
    if(whichdemo==8)
    {
      print(demo(topic = 'use_sppmix_MIPPP', package = 'sppmix'))
    } else
    {
      stop("Bad choice of demo. Enter a number from 1 to 8.")
    }
  }
}

#' General Helper Functions
#'
#' @description
#' The \code{plotstring} function plots a string in a generic plot device.
#' @param str A string to display.
#' @rdname helper_functions
#' @examples
#' \donttest{
#' plotstring()}
#'
#' @export
plotstring <- function(str="Hello World") {
  par(mai=c(1,1,2.1,1))
  plot(c(0,10),c(0,10),type="n",axes=FALSE,
       xlab="", ylab="",main=str)
  par(mai=c(1,1,1,1))
}

#' Matern covariance function
#'
#' @description
#' Computes the Matern covariance function. Used in the
#' creation of stationary and isotropic
#' Gaussian Random Fields (GRFs).
#'
#' For examples see
#'
#' \url{http://faculty.missouri.edu/~micheasa/sppmix/sppmix_all_examples.html
#' #MaternCov}
#'
#' @param grid An \code{nx2} matrix of locations
#' over which to compute the covariance matrix or
#' an \code{nxn} matrix representing the distances
#' of the \code{n} planar points.
#' @param nu,theta,sig Matern model parameters. See details.
#' @return A matrix representing the
#' covariance matrix for a random field (typically a GRF).
#' @seealso \code{\link{rGRF}}
#' @details The Matern covariance model
#' for two points with Euclidean
#' distance r, is given by
#'
#' C(r) = sig^2 2^(1-nu) gamma(nu)^(-1) (sqrt(2nu) r/theta)^nu B_nu(sqrt(2nu) r/theta)
#'
#' where sig, theta, nu>0, and B_nu is the
#' modified Bessel function of second kind. Note that
#' the Matern for nu=.5 reduces to the exponential.
#' @author Sakis Micheas
#' @examples
#' \donttest{
#' grid=cbind(seq(0,1,length=10), seq(0,1,length=10))
#' MaternCov(grid)}
#'
#' @export
MaternCov=function(grid,nu=.5,theta=1,sig=1)
{
#  mat=matern.image.cov(grid=grid,setup=TRUE,
#                   theta=theta,smoothness=v)
  if(ncol(grid)==2)
  #passed locations, need to compute distances
    r=as.matrix(dist(x=grid,diag=TRUE,upper=TRUE))
  else
    r=grid
  #Matern for nu=.5 reduces to exponential
  if(nu==.5)
    val=sig*sig*exp(-r/theta)
  else
    val=sig*sig*2^(1-nu)*(sqrt(2*nu)*r/theta)^nu*besselK(sqrt(2*nu)*r/theta,nu)/gamma(nu)
  return (val)
}

#' Generate a Gaussian Random Field
#'
#' @description
#' Generates Gaussian random fields
#' (GRFs) and related fields via transformations.
#' The spatial covariances are modeled using
#' Matern's model.
#'
#' For examples see
#'
#' \url{http://faculty.missouri.edu/~micheasa/sppmix/sppmix_all_examples.html
#' #rGRF}
#'
#' @param mu Mean of the stationary GRF.
#' @param gentype Set to 0 for Gaussian, 1 for Chi-square. Default
#' is \code{gentype=0}.
#' @param xlims,ylims Vectors defining the
#' grid limits of the x-y locations over
#' which to compute the covariance matrix.
#' @param LL Length of the side of the square grid.
#' @param df Degrees of freedom (an integer) for the
#' chi-square random field when \code{gentype=1}.
#' @param nu,theta,sig Matern model
#' parameters. See \code{\link{MaternCov}} for details.
#' @param pattern Optionally, a point pattern
#' as an object of type \code{\link[spatstat.geom]{ppp}}
#' containing locations within the window. The
#' values of the generated GRF over these
#' locations are returned as the marks of
#' the point pattern \code{pattern}.
#' @seealso \code{\link{MaternCov}},
#' \code{\link{plot_density}},
#' \code{\link[ggplot2]{ggtitle}},
#' \code{\link{add_title}}
#'
#' @return An image as an object of class \code{\link[spatstat.geom]{im.object}},
#' containing the realization of the field over the grid. If
#' argument \code{pattern} was supplied, the return
#' value is now a list contaning the realization of the field as
#' an image, augmented by the marked point pattern
#' with locations in \code{pattern} and marks the field
#' values over these locations. This capability is illustrated
#' for realizations of marked point processes
#' conditioning on continuous marks. See
#' function \code{rMIPPP_cond_loc} for more details.
#' @details The code of the \code{rGRF} function
#' uses a modification of the functions \code{\link[fields]{sim.rf}}
#' and \code{\link[fields]{matern.image.cov}} from the
#' \code{\link[fields]{fields}} package, by Douglas Nychka,
#' Reinhard Furrer, John Paige, and Stephan Sain.
#'
#' Depending on the choice of the Matern model parameters
#' we might end up having trouble with the FFT giving
#' negative values. The code accounts for this event and
#' adjusts the range of values via an increasing
#' variable \code{incr}. If it still takes a
#' long time to generate the fields try increasing the
#' domain of observation using wider \code{xlims}
#' and \code{ylims}.
#' @author Sakis Micheas
#' @examples
#' \donttest{
#' #Gaussian random field as an image
#' GRF1=rGRF()
#' p<-plot_density(as.data.frame(GRF1))
#' p_title<-expression( paste("GRF with Matern covariances, ", theta,"=1,",mu,"=0,",nu, "=.5,",
#'  sigma,"=1"))
#' p+ggplot2::ggtitle(p_title)
#' #or simply use the add_title function
#' p+add_title("GRF with Matern model covariances", mu=0,theta=1,nu=.5,sigma=1)
#' #Chi-Square random field as an image
#' ChiSqRF=rGRF(gentype=1,df=10)
#' p<-plot_density(as.data.frame(ChiSqRF))
#' p+add_title(paste(chi^{2}," random fields with Matern model covariances for the GRFs"),
#' mu=0,theta=1,nu=.5,sigma=1,df=10)
#' #Log-Gaussian random field as an image
#' GRF2=rGRF()
#' LogGRF=exp(rGRF())
#' p<-plot_density(as.data.frame(LogGRF))
#' p+add_title("Log-Gaussian random field with Matern model covariances", mu=0,theta=1,
#'  nu=.5,sigma=1)}
#'
#' @export
rGRF=function(mu=0,gentype=0,
  xlims=c(-5,5),ylims=c(-5,5),LL=128,
  df=10,nu=.5,theta=1,sig=1,pattern)
{
  #compute covariance matrix over grid
  cat("\nComputing distances and the Matern covariances...")
  grid=list(x=seq(xlims[1],xlims[2],length=LL),
            y=seq(ylims[1],ylims[2],length=LL))
  #if we passed a pattern add it to the
  #grid locations and update the LL value
  endgrid=LL*LL
  if(!missing(pattern))
  {
    grid$x=c(grid$x,pattern$x)
    grid$y=c(grid$y,pattern$y)
    LL=LL+pattern$n
  }
  dx <- grid$x[2] - grid$x[1]
  dy <- grid$y[2] - grid$y[1]
  m <- length(grid$x)
  n <- length(grid$y)
  ntot=LL*LL
  domainok=FALSE
  incr=2
  while(!domainok)
  {
    M <- ceiling2(incr * m)
    N <- ceiling2(incr * n)
    if (M%%2 != 0) {
      M <- M + 1
    }
    if (N%%2 != 0) {
      N <- N + 1
    }
    xGrid <- (1:M) * dx - (dx * M)/2
    yGrid <- (1:N) * dy - (dy * N)/2
    bigDistance <- sqrt(matrix(xGrid^2, M, N, byrow = FALSE) +
                          matrix(yGrid^2, M, N, byrow = TRUE))
    Gcov=MaternCov(grid=bigDistance,nu=nu,theta=theta,sig=sig)
    temp <- matrix(0, nrow = M, ncol = N)
    temp[M/2, N/2] <- 1
    wght <- fft(Gcov)/(fft(temp) * M * N)
    if (any(Re(wght) < 0))
    {
      incr=incr+1
      cat("\nFFT of covariance has negative values.\nWill increase the domain and try again. Expansion factor=",incr)
      M <- ceiling2(incr * m)
      N <- ceiling2(incr * n)
    }
    else
      domainok=TRUE
  }
  cat(" Done.\nGenerating from the GRF...")
  if(gentype==0)
  {
    z <- fft(matrix(mu+rnorm(N * M), ncol = N, nrow = M))
    Gfield=Re(fft(sqrt(wght) * z, inverse = TRUE))[1:m, 1:n]/sqrt(M*N)
    pxyG=matrix(Gfield,LL,LL,byrow=TRUE)
  }
  else if(gentype==1)
  {
    Gfield=rep(0,ntot)
    for(i in 1:ceiling(df))
    {
      z <- fft(matrix(mu+rnorm(N * M), ncol = N, nrow = M))
      GRF=Re(fft(sqrt(wght) * z, inverse = TRUE))[1:m, 1:n]/sqrt(M*N)
      Gfield=Gfield+GRF^2
    }
    pxyG=matrix(Gfield,LL,LL,byrow=TRUE)
  }
  cat(" Done.\n")
  if(missing(pattern))
    RETVAL<-as.im(list(x=grid$x,y=grid$y,z=pxyG))
  else
  {
#    markGRFondata=Gfield[(endgrid+1):ntot]
    npts=0
    marks=rep(0,pattern$n)
    done=FALSE
    for(pt in 1:pattern$n)
    {
      for(i in 1:length(grid$x))
      {
        if(abs(pattern$x[pt]-grid$x[i])<1e-05)
        {
          for(j in 1:length(grid$y))
          {
            if(abs(pattern$y[pt]-grid$y[j])<1e-05)
            {
              npts=npts+1
              marks[npts]=pxyG[i,j]
              if(npts==pattern$n)
              {
                done=TRUE
                break
              }
            }
          }
        }
        if(done)
          break
      }
      if(done)
        break
    }
    MPP=ppp(pattern$x,pattern$y,window=
              spatstat.geom::owin(xrange=xlims,yrange=ylims),
            check=FALSE,marks=marks)
    RETVAL<-list(dens_image=as.im(list(x=grid$x,
                        y=grid$y,z=pxyG)),
                MPP=MPP)
  }
  return (RETVAL)
}
