#' Fit a MIPPP conditionally on location
#'
#' @description
#' This function fits a Marked IPPP (MIPPP) on a marked
#' point pattern by modeling the (joint)
#' intensity surface of the locations and the marks
#' using an IPPP for the locations (independent
#' of the mark values) and for discrete marks a Gibbs model for the
#' mark distribution which is conditionally defined
#' on all the locations. NOTE: The estimation procedure for continuous
#' marks (random fields) will be implemented
#' in future versions of the \code{sppmix} package.
#'
#' For examples see
#'
#' \url{http://faculty.missouri.edu/~micheasa/sppmix/sppmix_all_examples.html
#' #est_MIPPP_cond_loc}
#'
#' @param pp Marked point pattern of class \code{ppp}.
#' @param r Radius used to define the neighborhood system. Any two locations within this distance are considered neighbors.
#' @param hyper Hyperparameter for the proposal distribution of gamma. This is currently the standard deviation for the random walk Metropolis-Hastings steps (one step for each gamma). Use a small value.
#' @param L Number of iterations for the MCMC; default is 10000.
#' @param burnin Number of initial realizations to discard. By default, it is 1/10 of the total number of iterations.
#' @param m The number of components to fit for the ground process when \code{fit_groundIPPP=TRUE}.
#' @param fit_groundIPPP Logical variable requesting to fit and return the DAMCMC results of the ground process.
#' @param truncate Logical variable indicating whether or not we
#' we only work with events within the window defined
#' in the point pattern \code{pp}.
#' @param grayscale Logical to request plots in grayscale.
#' @param discrete_mark Logical flag indicating whether the mark is a discrete (numerical value) or not.
#' For continuous marks set this to FALSE.
#' @param LL Length of the side of the square grid.
#' @param startgamma Initial value for the gamma vector. If missing the zero vector is used.
#' @param open_new_window Open a new window for a plot.
#' @param show_plots Logical variable requesting to produce the probability field plots for each mark.
#' @param ... Additional arguments for the S3 method.
#' @details We assume that the joint distribution of a
#' marked point pattern \code{N=[s,m(s)]} with \code{n}
#' events is of the form:
#'
#' \code{p(N)=lambda^n*exp(-lambda)/(n!)*f(all s|theta1)*g(all m|theta2(s),all s)}
#'
#' where \code{s} denotes a location and \code{m=m(s)}
#' a mark value at that location, lambda a parameter
#' with the interpretation as the average number of points
#' over the window of observation, and \code{f}, \code{g} are proper densities.
#'
#' The location (or ground process) \code{f(all s|theta1)}
#' can be fit using any method for unmarked
#' point patterns (for us it is modeled using
#' an IPPP with a mixture of normals
#' intensity surface). The function fits
#' the parameters of the second part of this model by default. However,
#' setting \code{fit_groundIPPP=TRUE} will fit a mixture intensity
#' surface for the ground process and return it for future processing.
#' Alternatively, simply retrieve the marked point process returned
#' and fit the ground process using \code{\link{est_mix_damcmc}}
#' or  \code{\link{est_mix_bdmcmc}} (the marks will be ignored and only
#' the locations will be used).
#'
#' Since \code{s} is observed over some
#' window and the marks are conditioned on
#' knowing the locations, then \code{g} is a
#' random field for each value of \code{m}.
#'
#' The neighborhood system is controlled by
#' \code{r} and is crucial in this case. Small values
#' tend to produce probability fields with concentrated
#' masses about observed events of the process,
#' whereas, large neighborhoods allow us to borrow
#' strength across locations and result in much smoother
#' probability fields. This parameter is currently NOT estimated
#' by the \code{\link{sppmix}} package, but will
#' be implemented in future releases.
#'
#' See Micheas (2014) for more details on special cases (2 marks) of
#' these Marked IPPP models via conditioning arguments.
#' @references Hierarchical Bayesian Modeling of Marked Non-Homogeneous Poisson Processes with finite mixtures and inclusion of covariate information. Micheas, A.C. (2014). Journal of Applied Statistics, 41, 12, 2596-2615, DOI: 10.1080/02664763.2014.922167.
#' @return An object of class \code{MIPPP_fit}, which is simply a list containing the following components:
#' \item{gen_gammas}{Posterior realizations of the gammas.}
#' \item{prob_fields}{Probability fields of marks.}
#' \item{discrete_mark}{Same logical flag as the input argument.}
#' \item{r}{Same as the input argument.}
#' \item{pp}{Same as the input argument.}
#' \item{ground_fit}{An object of type \code{damcmc_res} which contains the results of a DAMCMC fit to the ground process. If \code{fit_groundIPPP=FALSE} this is \code{NULL}.}
#' \item{condition_on_loc}{Logical variable indicating the type of conditioning used in order to produce this MIPPP fit. For this function it is set to TRUE.}
#' @author Sakis Micheas, Jiaxun Chen
#' @seealso \code{\link{GetStats}},
#' \code{\link{rMIPPP_cond_loc}}
#' @examples
#' \donttest{
#' # Create a marked point pattern
#' x <- runif(100)
#' y <- runif(100)
#' #mark distribution is discrete uniform
#' m <- sample(1:2, 100, replace=TRUE)
#' m <- factor(m, levels=1:2)
#' pp <- spatstat::ppp(x, y, c(0,1), c(0,1), marks=m)
#' # estimate the probability fields for each mark; since we have a discrete
#' # uniform for the mark distribution we should see probabilities about .5
#' # for both marks, as well as, the gamma credible sets should include 0,
#' # meaning that the marks are independent of location (probability .5 for
#' # each of the two mark values)
#' mpp_est <- est_MIPPP_cond_loc(pp, 0.1, hyper=0.2)
#' GetStats(mpp_est$gen_gammas[,1])$CredibleSet
#' GetStats(mpp_est$gen_gammas[,2])$CredibleSet
#' mpp_est <-est_MIPPP_cond_loc(pp, 0.3, hyper=0.2)
#' GetStats(mpp_est$gen_gammas[,1])$CredibleSet
#' GetStats(mpp_est$gen_gammas[,2])$CredibleSet
#' mpp_est <- est_MIPPP_cond_loc(pp, 0.5, hyper=0.2)
#' GetStats(mpp_est$gen_gammas[,1])$CredibleSet
#' GetStats(mpp_est$gen_gammas[,2])$CredibleSet
#' #Visualize the Tornado data about MO. We request to fit both the mark
#' #and ground process.
#' #plot the states, the tornado locations and the marks (strength of a tornado)
#' ret=PlotUSAStates(states=c('Iowa','Arkansas', 'Missouri','Illinois','Indiana',
#'  'Kentucky','Tennessee', 'Kansas','Nebraska','Texas','Oklahoma', 'Mississippi',
#'  'Alabama','Louisiana'), showcentroids=FALSE,shownames=TRUE, plotlevels = FALSE,
#'  main= "Tornadoes about MO, 2011")
#' #check out the mark values and their frequency
#' table(Tornadoes2011MO$marks)
#' #plot each point with a different shape according to its marks
#' ret$p+ggplot2::geom_point(data=as.data.frame( Tornadoes2011MO),ggplot2::aes(x=x,
#'  y=y,shape= as.factor(marks)))+ggplot2::guides(shape = ggplot2::guide_legend(
#'  title="Tornado power", ncol=2,byrow=TRUE))
#' #plot each point with a different color according to its marks
#' ret$p+ggplot2::geom_point(data=as.data.frame( Tornadoes2011MO),ggplot2::aes(x=x,
#'  y=y,color= as.factor(marks)))+ggplot2::guides(color = ggplot2::guide_legend(
#'  title="Tornado power", ncol=2,byrow=TRUE))
#' #plot each point with a different circle size according to its marks
#' ret$p+ggplot2::geom_point(data=as.data.frame( Tornadoes2011MO),ggplot2::aes(x=x,
#'  y=y,size=marks),shape=21)+ ggplot2::scale_size_continuous(breaks=sort(unique(
#'  Tornadoes2011MO$marks))) + ggplot2::guides(size =ggplot2::guide_legend(title=
#'  "Tornado power", ncol=2,byrow=TRUE))
#' # the marks must start from 1, recode the original
#' Tornadoes2011MO1=Tornadoes2011MO
#' Tornadoes2011MO1$marks=Tornadoes2011MO1$marks+1
#' mpp_est=est_MIPPP_cond_loc(Tornadoes2011MO1,r=1.5,hyper=0.01,
#'  startgamma = c(.1,.2,.3,.4,.5,.6),fit_groundIPPP=TRUE)
#' #Now generate an MIPPP with 3 marks
#' newMPP=rMIPPP_cond_loc(gammas=c(.1,.2,.5))
#' summary(newMPP)
#' plot(newMPP$surf,main="True IPPP intensity surface for the locations")
#' true_gammas=newMPP$gammas
#' genMPP=newMPP$genMPP
#' newMPP$r
#' mpp_est=est_MIPPP_cond_loc(genMPP,newMPP$r, hyper=0.2)
#' GetStats(mpp_est$gen_gammas[,1])$Mean
#' GetStats(mpp_est$gen_gammas[,2])$Mean
#' GetStats(mpp_est$gen_gammas[,3])$Mean
#' GetStats(mpp_est$gen_gammas[,1])$CredibleSet
#' GetStats(mpp_est$gen_gammas[,2])$CredibleSet
#' GetStats(mpp_est$gen_gammas[,3])$CredibleSet
#' summary(mpp_est)
#' plot(mpp_est)
#' plot(mpp_est,newMPP$surf)}
#'
#' @export
est_MIPPP_cond_loc <- function(pp, r, hyper=1, L = 10000,
     burnin= floor(L/10), m=3, fit_groundIPPP=FALSE,
     truncate=FALSE, grayscale=FALSE,startgamma,
     discrete_mark = TRUE, LL = 150,
     open_new_window = FALSE, show_plots = TRUE)
{
  x=y=value=NULL
  if(grayscale)
    cols <- gray.colors(100,start = 1, end = 0)
  else
    cols <- c("#00007F", "blue", "#007FFF", "cyan", "#7FFF7F", "yellow",
              "#FF7F00", "red", "#7F0000")
  if(!discrete_mark)
    stop("Option not implemented yet.")
  window = pp$window
  # test if this pattern has marks
  if(length(pp$marks) == 0) {
    stop("This point pattern doesn't have any marks.")
  }
  # get number of points and the location of points
  n <- pp$n
  pattern <- cbind(pp$x, pp$y)
  if(discrete_mark)
  {
    # get the levels of the mark
    ppmarks <- pp$marks
    if (is.factor(pp$marks)) {
      marks <- as.numeric(levels(pp$marks))
    } else {
      marks <- sort(unique(pp$marks))
    }
    nmarks <- length(marks)
    if(missing(startgamma))
      startgamma=rep(0,nmarks)
    hyperparams=rep(0,nmarks+1)
    hyperparams[1]=hyper
    hyperparams[2:(nmarks+1)]=startgamma
    MPPfit=MIPPCondLoc_sppmix(pattern,pp$marks,
      pp$window$xrange,pp$window$yrange,
      L,truncate,hyperparams,marks,discrete_mark,r)

    meangammas=colMeans(MPPfit$gen_gammas[
      -(1:burnin), ])

    probfields=GetProbFieldsCondLoc_sppmix(pattern,pp$marks,
      pp$window$xrange,pp$window$yrange,
      LL,meangammas,marks,truncate,r)

    ps_fields <- list()
    for (m in 1:nmarks)
    {
      ps_fields[[m]]<-as.data.frame(as.im(list(x=probfields$x,y=probfields$y,z=probfields$ps_fields[,, m])))
    }
  }
  else
  #continuous marks need to estimate the
  #random field parameters
  {

  }

  if(show_plots)
  {
    if(discrete_mark)
    {
      for (m in 1:nmarks)
      {
        openwin_sppmix(check2open=open_new_window)
        pp_df <- data.frame(pp$x,pp$y)
        names(pp_df) <- c("x", "y")
        p <- ggplot2::ggplot(ps_fields[[m]],aes(x, y))+
            ggplot2::labs(x = "x", y = "y")+
            ggplot2::theme_classic() +
            ggplot2::theme(panel.border =ggplot2::element_rect(fill = NA,size = 1))

        p<-p + ggplot2::geom_raster(aes(fill = value), interpolate = TRUE) +
          ggplot2::scale_fill_gradientn(
            colors = cols,limits=c(0,1)) +
          ggplot2::guides(fill = guide_colorbar(
            title = "Probability",nbin = 100, barheight = 15))+
          ggplot2::coord_cartesian(xlim =pp$window$xrange,
                                   ylim =pp$window$yrange,expand=FALSE)

        title1=paste("Probability field of observing mark",
                     marks[m],", r =",r)
        gamma1=bquote(gamma==.(meangammas[m]))
        nmarks1=bquote(paste("total number of marks is ",.(nmarks)))
        all_char=list(gamma1,nmarks1)
        cal=do.call(function(...) substitute(list(...)),all_char)

        p<-p+ggplot2::geom_point(data = pp_df,size=1)+
          #  ggplot2::ggtitle(paste("Probability of observing mark",
          #                        marks[m],", r=",r))
          ggtitle(substitute(atop(title1, cal)))
        print(p)

      }
    }
    else
    {#continuous marks

    }
  }

  if(fit_groundIPPP==TRUE)
    fitDA<-est_mix_damcmc(pp,m=m,truncate=truncate)
  else
    fitDA=NULL

  RVAL <- list(gen_gammas = MPPfit$gen_gammas,
               prob_fields = ps_fields,
               discrete_mark = discrete_mark,
               r=r,
               pp=pp,
               ground_fit=fitDA,
               condition_on_loc=TRUE)
  class(RVAL) <- c("MIPPP_fit")
  return(RVAL)
}

#' Plot the mark probabilities of a marked point pattern
#'
#' @description
#' For discrete marks only, the function
#' displays for each location of a marked
#' point pattern, the probabilities of observing
#' each of the discrete marks at that location.
#'
#' For examples see
#'
#' \url{http://faculty.missouri.edu/~micheasa/sppmix/sppmix_all_examples.html
#' #plot_MPP_probs}
#'
#' @param MPPfit Object of class \code{MIPPP_fit}.
#' @param truncate Logical variable indicating to discard points
#' if they are not within the window of observation. Default is FALSE.
#' @param open_new_window Open a new window for the plot.
#'
#' @seealso \code{\link{rMIPPP_cond_loc}},
#' \code{\link{est_MIPPP_cond_loc}}
#' @author Sakis Micheas
#' @examples
#' \donttest{
#' newMPP=rMIPPP_cond_loc(gammas=c(.1,.2,.5))
#' plot(newMPP$surf,main="True IPPP intensity surface for the locations")
#' genMPP=newMPP$genMPP
#' newMPP$r
#' mpp_est <- est_MIPPP_cond_loc(genMPP,newMPP$r, hyper=0.2)
#' plot_MPP_probs(mpp_est)}
#'
#' @export
plot_MPP_probs <- function(MPPfit,
   truncate = FALSE,open_new_window=FALSE)
{
  openwin_sppmix(check2open=open_new_window)
  point=component=probability=NULL

  marks=unique(MPPfit$pp$marks)
  nmarks=length(marks)
  pattern <- cbind(MPPfit$pp$x, MPPfit$pp$y)
  meangammas=apply(MPPfit$gen_gammas,2,mean)
  gensm=GetProbCondLoc_sppmix(
    pattern,MPPfit$pp$marks,
    MPPfit$pp$window$xrange,
    MPPfit$pp$window$yrange,
    meangammas,marks,truncate,MPPfit$r)

  labsy=1:nmarks
  labsx=as.integer(seq(1,MPPfit$pp$n,length=10))

  probs=gensm$probs
  plot_df <- data.frame(probability = as.vector(probs),
                        point = 1:MPPfit$pp$n,
                        component = rep(1:nmarks, each = MPPfit$pp$n))

  ggplot(plot_df, aes(point, component, xend = point,
                      yend = component + 1,
                      col = probability,
                      xmin=1,xmax=MPPfit$pp$n,
                      ymin=1,ymax=nmarks)) +
    geom_segment(size = I(5)) +
    ggplot2::scale_color_gradient(low = "white", high = "grey18",
                                  limits = c(0,1),breaks=c(0,.25,.5,.75,1),            guide = guide_colorbar(nbin = 100, barheight = 15)) +
    ggplot2::theme_classic() +
    ggplot2::theme(panel.border = ggplot2::element_rect(fill = NA, size = 1)) +
    ggplot2::labs(x = "Point", y = "Mark", colour = "Probability") +
    ggplot2::scale_x_discrete(limits = labsx)+
    ggplot2::scale_y_discrete(limits = labsy+.5,labels=labsy)+
    add_title("Mark probabilities (estimated)", n = MPPfit$pp$n)
}

#' Plot the mark probability fields
#'
#' @description
#' The function displays the mark probability fields
#' for each location of a marked
#' point pattern. These fields are simply the
#' probabilities of observing the corresponding
#' mark value at that location.
#'
#' For examples see
#'
#' \url{http://faculty.missouri.edu/~micheasa/sppmix/sppmix_all_examples.html
#' #plot_MPP_fields}
#'
#' @param MPP A marked point pattern as an
#' object of class \code{\link[spatstat]{ppp}}.
#' @param discrete_mark Logical flag indicating whether the mark is discrete or not.
#' Default is TRUE. For continuous marks set this to FALSE.
#' @param gammas For discrete marks (\code{discrete_mark=TRUE}), this is
#' a vector of length equal to the number of marks.
#' These parameters should typically be non-negative and
#' they represent weights affecting the probability
#' fields of each mark. For values close to 0, we
#' get higher probabilities of observing this mark.
#' Large positive values lead to small probabilities of observing
#' the corresponding mark.
#' @param r Radius used to define the
#' neighborhood system. Any two locations
#' within this distance are considered
#' neighbors.
#' @param truncate Logical variable indicating to discard points
#' if they are not within the window of observation. Default is FALSE.
#' @param grayscale Logical to request plots in grayscale.
#' @param open_new_window Logical requesting a new window
#' for the plot(s).
#' @param LL Length of the side of the square grid.
#' The larger this value is, the better the picture resolution.
#'
#' @seealso \code{\link{rMIPPP_cond_loc}}
#' @author Sakis Micheas
#' @examples
#' \donttest{
#' newMPP=rMIPPP_cond_loc(gammas=c(.1,.2,.5), r=.5)
#' plot(newMPP$surf,main="True IPPP intensity surface for the locations")
#' plot_MPP_fields(newMPP$genMPP,newMPP$gammas,newMPP$r)
#' plot_MPP_fields(newMPP$genMPP,newMPP$gammas,1)
#' plot_MPP_fields(newMPP$genMPP,newMPP$gammas,1.5)
#' plot_MPP_fields(newMPP$genMPP,newMPP$gammas,2)}
#'
#' @export
plot_MPP_fields <- function(MPP,gammas,r,
    discrete_mark=TRUE,grayscale=FALSE,
    truncate = FALSE,open_new_window=FALSE,
    LL=128)
{
  x=y=value=NULL
  xlims=MPP$window$xrange
  ylims=MPP$window$yrange
  pattern <- cbind(MPP$x, MPP$y)
  if(grayscale)
    cols <- gray.colors(100,start = 1, end = 0)
  else
    cols <- c("#00007F", "blue", "#007FFF", "cyan", "#7FFF7F", "yellow",
              "#FF7F00", "red", "#7F0000")
  if(discrete_mark)
  {
    nmarks=length(gammas)
    marks=1:nmarks

    probfields=GetProbFieldsCondLoc_sppmix(
      pattern,MPP$marks,xlims,ylims,LL,
      gammas,marks,truncate,r)
    ps_fields <- list()
    for (m in 1:nmarks)
    {
      ps_fields[[m]]<-as.data.frame(
        as.im(list(x=probfields$x,y=
                     probfields$y,z=
                     probfields$ps_fields[,, m])))
    }
    gen_marks <- factor(MPP$marks, levels=1:nmarks)

    for (m in 1:nmarks)
    {
      openwin_sppmix(check2open=open_new_window)
      pp_df <- data.frame(MPP$x,MPP$y)
      names(pp_df) <- c("x", "y")

      p <- ggplot2::ggplot(ps_fields[[m]],aes(x, y))+
        ggplot2::labs(x = "x", y = "y")+
        ggplot2::theme_classic() +
        ggplot2::theme(panel.border =ggplot2::element_rect(fill = NA,size = 1))

      p<-p + ggplot2::geom_raster(aes(fill = value), interpolate = TRUE) +
        ggplot2::scale_fill_gradientn(
          colors = cols,limits=c(0,1)) +
        ggplot2::guides(fill = guide_colorbar(
          title = "Probability",nbin = 100, barheight = 15))+
        ggplot2::coord_cartesian(xlim =xlims,
                                 ylim =ylims,expand=FALSE)

      title1=paste("Probability field of observing mark",
                   marks[m],", r =",r)
      gamma1=bquote(gamma==.(gammas[m]))
      nmarks1=bquote(paste("total number of marks is ",.(nmarks)))
      all_char=list(gamma1,nmarks1)
      cal=do.call(function(...) substitute(list(...)),all_char)

      p<-p+ggplot2::geom_point(data = pp_df,size=1)+
        #  ggplot2::ggtitle(paste("Probability of observing mark",
        #                        marks[m],", r=",r))
        ggtitle(substitute(atop(title1, cal)))

      print(p)
    }
  }
}

#' Generate a Marked Poisson point process (conditional on location)
#'
#' @description
#' This function generates realizations (point patterns) from
#' a given Marked IPPP or a generated one. See details for the choice of models
#' for the mark distribution. The location (ground) process is
#' a standard IPPP (unmarked) with mixture intensity surface, and is responsible
#' for the number of events in the point pattern.
#'
#' For examples see
#'
#' \url{http://faculty.missouri.edu/~micheasa/sppmix/sppmix_all_examples.html
#' #rMIPPP_cond_loc}
#'
#' @param surf An object of type \code{intensity_surface} representing the
#' IPPP surface for the ground process. Omit this argument to create a surface randomly.
#' @param locPP The ground IPPP (locations of the events). If missing then these
#' are generated using a call to \code{\link{rsppmix}}. Note that if
#' \code{surf} is not supplied, then it will be generated which may lead to
#' completely inappropriate locations of the events, if the supplied \code{locPP}
#' was created with a completely different surface. It is safer to supply both the surface
#' and ground locations at the same time or none of the two, so that both will be generated.
#' @param gammas For discrete marks (\code{discrete_mark=TRUE}), this is
#' a vector of length equal to the number of marks.
#' These parameters should typically be non-negative and
#' they represent weights affecting the probability
#' fields of each mark. For values close to 0, we
#' get higher probabilities of observing this mark.
#' Large positive values lead to small probabilities of observing
#' the corresponding mark. Negative values are allowed, but they can lead
#' to a mark not being present in the generated pattern.
#' If the vector \code{gammas} is not supplied,
#' then we randomly generate the number of marks from \code{1:10}
#' and the values of the vector \code{gammas}
#' from a gamma distribution.
#' @param r Radius used to define the
#' neighborhood system. Any two locations
#' within this distance are considered
#' neighbors. If missing, we randomly select
#' the radius using the generated (ground) point
#' pattern over the window parameter \code{win}.
#' @param hyper Hyperparameter for the distribution of gamma.
#' @param truncate Logical variable indicating whether or not we
#' normalize the densities of the mixture components
#' to have all their mass within the window defined
#' in the window \code{win}. This affects the mixture model for the
#' intensity surface of the ground process.
#' @param win Object of type \code{\link[spatstat]{owin}} defining the window of observation.
#' @param bigwin Object of type \code{\link[spatstat]{owin}}. If supplied, this will be the
#' window of observation, even if the pattern is generated over \code{win}. Useful if we
#' do not truncate (\code{truncate=FALSE}) and we want better presentation of the generated MIPPP.
#' @param discrete_mark Logical flag indicating whether the mark is discrete or not.
#' Default is TRUE. For continuous marks set this to FALSE.
#' @param open_new_window Open a new window for a plot.
#' @param grayscale Logical to request plots in grayscale.
#' @param show_plots Logical variable requesting to produce exploratory plots of the
#' Marked IPPP intensity surface and generated point pattern.
#' @param mark_distr_choice A number indicating which
#' mark distribution to use. Currently we have
#' only one choice in the discrete mark case, which is essentialy a Markov random field (MRF)
#' over the window. See details for more on the mark model currently used. For continuous marks,
#' we have two choices, Gaussian random field (GRF) for
#' \code{mark_distr_choice=0} or Chi-Square random field for
#' \code{mark_distr_choice=1}.
#' @param LL Length of the side of the square grid.
#' The larger this value is, the better the picture resolution.
#' @param L Number of iterations. Required when sampling from
#' the mark model conditional on locations.
#' @param GRFmu This is the mean of the
#' Gaussian random field. Only stationarity
#' is currently supported (i.e., \code{GRFmu} does
#' not depend on location). Used only if
#' \code{discrete_mark=FALSE}.
#' @param df Degrees of freedom (an integer) for the
#' chi-square random field when \code{mark_distr_choice=1}. Default is \code{df=10}. Used only if
#' \code{discrete_mark=FALSE}.
#' @param nu,theta,sig Additional arguments passed to the
#' \code{\link{MaternCov}} function in order to
#' create the spatial covariance field.
#' Default values are \code{nu=.5},
#' \code{theta=1}, and \code{sig=1}.
#' See \code{\link{MaternCov}} for details. Used only if
#' \code{discrete_mark=FALSE}.
#' @details We assume that the joint distribution of a
#' marked point pattern \code{N=[s,m(s)]} with \code{n}
#' events is of the form:
#'
#' \code{p(N)=lambda^n*exp(-lambda)/(n!)*f(all s|theta1)*g(all m|theta2(s),all s)}
#'
#' where \code{s} denotes a location and \code{m=m(s)}
#' a mark value at that location, lambda a parameter
#' with the interpretation as the average number of points
#' over the window of observation, and \code{f}, \code{g} are proper densities.
#'
#' In order to simulate from this Marked IPPP
#' we first simulate the number of events
#' and their locations from an IPPP with
#' mixture intensity surface \code{lambda*f(s|theta1)} (e.g.,
#' using \code{\link{rsppmix}}), and then generate
#' the mark at that location \code{s}.
#'
#' In the discrete mark case, the mark is modeled using
#' a mixture distribution of Dirac measures on
#' the marks with the probability \code{q(m,s)} of observing a
#' specific mark value \code{m} depending on the current location
#' \code{s} and the marks of its neighbors. Since
#' we have a window of observation, any point in there
#' can potentially be marked, which leads to \code{q(m,s)} being
#' a field. In particular, the probability \code{q(m,s)} is analogous to
#'
#' \code{exp(-gammas_(j)*(sum over all neighbors of s of their marks minus m squared))}
#'
#' and when we fit the MIPPP model, our goal
#' is to estimate the parameters \code{gammas}.
#'
#' Note that if all \code{gammas} are zero then
#' we fall back to a discrete uniform mark distribution.
#'
#' The neighborhood system is controlled by
#' \code{r} and is crucial in this case. Small values
#' tend to produce probability fields with concentrated
#' masses about observed events of the process,
#' whereas, large neighborhoods allow us to borrow
#' strength across locations and result in much smoother
#' probability fields.
#'
#' In the continuous case the mark is generated from
#' a (typically stationary) Gaussian process or chi-squared random process,
#' e.g., using function \code{\link{rGRF}}.
#'
#' See Micheas (2014) for more details on
#' Marked IPPP models via conditioning arguments.
#' @references Hierarchical Bayesian Modeling of Marked Non-Homogeneous Poisson Processes with finite mixtures and inclusion of covariate information. Micheas, A.C. (2014). Journal of Applied Statistics, 41, 12, 2596-2615, DOI: 10.1080/02664763.2014.922167.
#' @return A list containing the following components:
#' \item{surf}{The generated or supplied intensity surface object \code{surf} for the ground process.}
#' \item{gammas}{The generated or supplied parameters \code{gammas}. Returned only if \code{discrete_mark=TRUE}.}
#' \item{genMPP}{The generated point pattern as an object of class \code{\link[spatstat]{ppp}} and \code{sppmix}. The member \code{$marks} contains the marks at each of the generated locations. If the ground PP \code{locPP} was supplied, this is also the ground process for the MIPPP and only the marks are generated (at those locations).}
#' \item{r}{The generated or supplied parameter \code{r}. Returned only if \code{discrete_mark=TRUE}.}
#' \item{prob_fields}{In the continuous mark case this is the realization of the random field (as an image \code{\link[spatstat]{im}} object). For discrete marks, this is a list of size equal to the number of marks containing the probability fields for each mark value.}
#' \item{prob_field_params}{A list of the parameters used to create the continuous valued mark fields.  Returned only if \code{discrete_mark=FALSE}.}
#' @author Sakis Micheas
#' @seealso \code{\link{plot_MPP_fields}}
#' @examples
#' \donttest{
#' # Create a marked point pattern; use randomization and discrete marks (default values)
#' newMPP=rMIPPP_cond_loc()
#' plot(newMPP$surf,main="True IPPP intensity surface for the locations")
#' newMPP$gammas
#' newMPP$genMPP
#' newMPP$r
#' print(table(newMPP$genMPP$marks))
#' #we can reproduce the random field plots anytime using the following call
#' plot_MPP_fields(newMPP$genMPP,newMPP$gammas,newMPP$r)
#' #Now generate continuous marks according to a Gaussian process
#' newMPP=rMIPPP_cond_loc(discrete_mark = FALSE)
#' plot(newMPP$surf,main="True IPPP intensity surface for the locations")
#' #now the marks are taken from a chi-square field
#' newMPP=rMIPPP_cond_loc(mark_distr_choice=1, discrete_mark = FALSE)
#' plot(newMPP$surf,main="True IPPP intensity surface for the locations")}
#'
#' @export
rMIPPP_cond_loc <- function(surf,locPP, gammas, r,
            hyper=0.01, truncate=FALSE,
            win=owin(c(-3,3),c(-3,3)),
            bigwin,
            discrete_mark = TRUE,
            open_new_window = FALSE,
            grayscale=FALSE,
            show_plots = TRUE,LL=128,L=50000,
            mark_distr_choice=0,GRFmu=0,
            df=10,nu=.5,theta=1,sig=1)
{
  x=y=value=chi=NULL
  if(grayscale)
    cols <- gray.colors(100,start = 1, end = 0)
  else
    cols <- c("#00007F", "blue", "#007FFF", "cyan", "#7FFF7F", "yellow",
              "#FF7F00", "red", "#7F0000")
  xlims=win$xrange
  ylims=win$yrange
  if(missing(surf))
    surf=rmixsurf(5,200,rand_m = TRUE,
                  xlim=xlims, ylim=ylims)
  #generate the locations if not supplied
  if(missing(locPP))
    genPP=suppressWarnings(
      rsppmix(surf,truncate = truncate))
  else
    genPP=locPP
  if(!missing(bigwin))
  {
    surf=to_int_surf(surf, win = bigwin)
    win=bigwin
    genPP$window$xrange=win$xrange
    genPP$window$yrange=win$yrange
    xlims=win$xrange
    ylims=win$yrange
  }
  if(missing(r))
    r=runif(1,minnndist(genPP),maxnndist(genPP))
  n <- genPP$n
  pattern <- cbind(genPP$x, genPP$y)
  if(discrete_mark)
  {
    if(missing(gammas))
      gammas=rgamma(sample(1:10,1,TRUE),10,.1)
    nmarks=length(gammas)
    marks=1:nmarks
    #generate from the MRF
    #start with random marks at each location
    gensm=GenMarksProbCondLoc_sppmix(
      pattern,L,xlims,ylims,gammas,
      marks,truncate,r)
    gen_marks=gensm$marks
    print(table(gen_marks))
    probfields=GetProbFieldsCondLoc_sppmix(
      pattern,gen_marks,xlims,ylims,LL,
      gammas,marks,truncate,r)
    ps_fields <- list()
    for (m in 1:nmarks)
    {
#      cat("\nm=",m,"\n")
#      print(table((probfields$ps_fields[,,m]>=0
#             & probfields$ps_fields[,,m]<=1)
#            ,useNA= "ifany"))
      ps_fields[[m]]<-as.data.frame(
        as.im(list(x=probfields$x,y=
        probfields$y,z=
        probfields$ps_fields[,, m])))
        #,na.replace=0))
    }
    gen_marks <- factor(gen_marks, levels=1:nmarks)
  }
  else
  {
    gammas=NULL
    genGRF<-rGRF(mu=GRFmu,
      gentype=mark_distr_choice,xlims=xlims,
      ylims=ylims,LL=LL,df=df,nu=nu,
      theta=theta,sig=sig,pattern=genPP)
    prob_fields<-genGRF$dens_image
    gen_marks <- genGRF$MPP$marks
    if(mark_distr_choice==0)
      prob_field_params=list(GRFmu=GRFmu,
          nu=nu,theta=theta,sigma=sig)
    else
      prob_field_params=list(GRFmu=GRFmu,
          nu=nu,theta=theta,sigma=sig,df=df)
  }

  genMPP <- ppp(genPP$x, genPP$y,window=win,
                marks=gen_marks,check=FALSE)

  if(show_plots)
  {
    p<-plotmix_2d(intsurf=surf,open_new_window = open_new_window)
    if(discrete_mark)
    {
      p <-p+ ggplot2::geom_point(
        data=as.data.frame(genMPP),
        aes(color =marks))+guides(
          color=guide_legend(
            title="Mark",ncol=2,
            byrow=TRUE))+
        add_title("Locations and ground intensity surface of a marked IPPP",
                       lambda = surf$lambda,
                       m = surf$m,nmarks=nmarks)
      print(p)

      for (m in 1:nmarks)
      {
        openwin_sppmix(check2open=open_new_window)
        pp_df <- data.frame(genPP$x,genPP$y)
        names(pp_df) <- c("x", "y")

        p <- ggplot2::ggplot(ps_fields[[m]],aes(x, y))+
          ggplot2::labs(x = "x", y = "y")+
          ggplot2::theme_classic() +
          ggplot2::theme(panel.border =ggplot2::element_rect(fill = NA,size = 1))

        p<-p + ggplot2::geom_raster(aes(fill = value), interpolate = TRUE) +
          ggplot2::scale_fill_gradientn(
            colors = cols,limits=c(0,1)) +
          ggplot2::guides(fill = guide_colorbar(
            title = "Probability",nbin = 100, barheight = 15))+
          ggplot2::coord_cartesian(xlim =xlims,
                                   ylim =ylims,expand=FALSE)

        title1=paste("Probability field of observing mark",
              marks[m],", r =",r)
        gamma1=bquote(gamma==.(gammas[m]))
        nmarks1=bquote(paste("total number of marks is ",.(nmarks)))
        all_char=list(gamma1,nmarks1)
        cal=do.call(function(...) substitute(list(...)),all_char)

        p<-p+ggplot2::geom_point(data = pp_df,size=1)+
        #  ggplot2::ggtitle(paste("Probability of observing mark",
         #                        marks[m],", r=",r))
          ggtitle(substitute(atop(title1, cal)))
        print(p)
      }
    }
    else
    {
      p <-p+ ggplot2::geom_point(
        data=as.data.frame(genMPP))+
        add_title("Locations and ground intensity surface of a marked IPPP",
                       lambda = surf$lambda,
                       m = surf$m)
      print(p)

      openwin_sppmix(check2open=open_new_window)
      p<-plot_density(as.data.frame(prob_fields))+
        guides(fill = guide_colorbar(
        title = "Mark value",nbin = 100, barheight = 15))

      if(mark_distr_choice==0)
        p<-p+add_title("Marks from a GRF with Matern model covariances",
                    mu=GRFmu,theta=theta,
                    nu=nu,sigma=sig)
      else if(mark_distr_choice==1)
        p<-p+add_title(paste("Marks from a ",chi^{2}," random field (Matern model for the GRF covariances)"),
                  mu=GRFmu,theta=theta,
                  nu=nu,sigma=sig,df=df)
      print(p)
    }
  }
  if(discrete_mark)
    RVAL <- list(surf = surf,
               gammas = gammas,
               genMPP= genMPP,
               r= r,
               prob_fields=ps_fields)
  else
    RVAL <- list(surf = surf,
               genMPP= genMPP,
               prob_fields=prob_fields,
               prob_field_params=prob_field_params)
  return(RVAL)
}

#' @rdname est_MIPPP_cond_loc
#' @param x An object of class \code{MIPPP_fit} (the result of a request from \code{\link{est_MIPPP_cond_loc}}).
#' @param surf An object of type \code{intensity_surface} representing the
#' IPPP surface for the ground process (conditioning on location only). This can be the surface of posterior means
#' or the MAP from a \code{damcmc_res} object or the MAP number of components
#' surface from a \code{bdmcmc_res} object.
#' @param main Main title for the plot.
#'
#' @export
#' @method plot MIPPP_fit
plot.MIPPP_fit<- function(x,surf,
                    open_new_window = FALSE,
                    grayscale=FALSE,
                    main="Locations and ground intensity surface of a marked IPPP",...)
{
  y=value=NULL
  MPPfit<-x
  if(MPPfit$condition_on_loc)
  {
    ps_fields=MPPfit$prob_fields
    pp=MPPfit$pp
    r=MPPfit$r
    win = pp$window
    # test if this pattern has marks
    if(length(pp$marks) == 0) {
      stop("This point pattern doesn't have any marks.")
    }
    if(grayscale)
      cols <- gray.colors(100,start = 1, end = 0)
    else
      cols <- c("#00007F", "blue", "#007FFF", "cyan", "#7FFF7F", "yellow",
                "#FF7F00", "red", "#7F0000")
    # get the levels of the mark
    ppmarks <- pp$marks
    if (is.factor(pp$marks)) {
      marks <- as.numeric(levels(pp$marks))
    } else {
      marks <- sort(unique(pp$marks))
    }
    nmarks <- length(marks)
    xlims=win$xrange
    ylims=win$yrange
    if(!missing(surf))
    {
      p<-plotmix_2d(intsurf=surf,open_new_window = open_new_window)
      p <-p+ ggplot2::geom_point(
        data=as.data.frame(pp),
        aes(color =marks))+guides(color =
                                    guide_legend(title="Mark",ncol=2,
                                                 byrow=TRUE))
      all_char=list(bquote(paste(lambda == .(surf$lambda))),
      bquote(paste(m == .(surf$m), " components")),
      bquote(paste(.(nmarks)," marks")))
      cal <- do.call(function(...) substitute(list(...)), all_char)
      p<-p + ggtitle(substitute(atop(main, cal)))
      print(p)
    }

    for (m in 1:nmarks)
      {
        openwin_sppmix(check2open=open_new_window)
        pp_df <- data.frame(pp$x,pp$y)
        names(pp_df) <- c("x", "y")
        p <- ggplot2::ggplot(ps_fields[[m]],aes(x, y))+
          ggplot2::labs(x = "x", y = "y")+
          ggplot2::theme_classic() +
          ggplot2::theme(panel.border =ggplot2::element_rect(fill = NA,size = 1))

        p<-p + ggplot2::geom_raster(aes(fill = value), interpolate = TRUE) +
          ggplot2::scale_fill_gradientn(
            colors = cols,limits=c(0,1)) +
          ggplot2::guides(fill = guide_colorbar(
            title = "Probability",nbin = 100, barheight = 15))+
          ggplot2::coord_cartesian(xlim =xlims,
                                   ylim =ylims,expand=FALSE)

        gammas=apply(MPPfit$gen_gammas,2,mean)
        title1=paste("Probability field of observing mark",
                     marks[m],", r =",r)
        gamma1=bquote(gamma==.(gammas[m]))
        nmarks1=bquote(paste("total number of marks is ",.(nmarks)))
        all_char=list(gamma1,nmarks1)
        cal=do.call(function(...) substitute(list(...)),all_char)

        p<-p+ggplot2::geom_point(data = pp_df,size=1)+
          #  ggplot2::ggtitle(paste("Probability of observing mark",
          #                        marks[m],", r=",r))
          ggtitle(substitute(atop(title1, cal)))

        print(p)

    }
  }
  else
  {
    MPP=MPPfit$pp
    win = MPP$window
    # test if this pattern has marks
    if(length(MPP$marks) == 0) {
      stop("This point pattern doesn't have any marks.")
    }
    if(grayscale)
      cols <- gray.colors(100,start = 1, end = 0)
    else
      cols <- c("#00007F", "blue", "#007FFF", "cyan", "#7FFF7F", "yellow",
                "#FF7F00", "red", "#7F0000")
    # get the levels of the mark
    ppmarks <- MPP$marks
    if (is.factor(ppmarks)) {
      marks <- as.numeric(levels(ppmarks))
    } else {
      marks <- sort(unique(ppmarks))
    }
    nmarks <- length(marks)
    xlims=win$xrange
    ylims=win$yrange
    for (mark in 1:nmarks)
    {
      if(MPPfit$fit_DAMCMC)
      {
        openwin_sppmix(check2open=open_new_window)
        print(plot_density(MPPfit$ground_fitsAoS[[mark]],main=paste("AoS surface for mark",mark,", x denotes the estimated posterior means"),
                   pp=MPPfit$pp[MPPfit$pp$marks==mark]
                   ,surf=MPPfit$post_surf[[mark]],grayscale = grayscale))
        plotmix_3d(MPPfit$ground_fitsAoS[[mark]],
                   title1=paste("AoS surface for mark",mark),grayscale = grayscale)
      }
      else
      {
        openwin_sppmix(check2open=open_new_window)
        print(plot_density(MPPfit$ground_fitsAoS[[mark]],main=paste("BMA surface for mark",mark,", x denotes the estimated posterior means"),
                     pp=MPPfit$pp[MPPfit$pp$marks==mark]
                     ,surf=MPPfit$post_surf[[mark]],grayscale = grayscale))
        plotmix_3d(MPPfit$ground_fitsAoS[[mark]],
                   title1=paste("BMA surface for mark",mark),grayscale = grayscale)
      }
    }
  }

}

#' @rdname est_MIPPP_cond_loc
#' @param object An object of class \code{MIPPP_fit} (the result of a request from \code{\link{est_MIPPP_cond_loc}}).
#'
#' @export
#' @method summary MIPPP_fit
summary.MIPPP_fit <- function(object,...)
{
  MPPfit=object
  if(MPPfit$condition_on_loc)
  {
    cat('\nNeighborhood r=',MPPfit$r)
    pp=MPPfit$pp
    # test if this pattern has marks
    if(length(pp$marks) == 0) {
      stop("This point pattern doesn't have any marks.")
    }
    # get the levels of the mark
    ppmarks <- pp$marks
    if (is.factor(pp$marks)) {
      marks <- as.numeric(levels(pp$marks))
    } else {
      marks <- sort(unique(pp$marks))
    }
    nmarks <- length(marks)
    # get number of points and the location of points
    summary(pp)
    cat('\nTotal number of marks:',nmarks)
    if(MPPfit$discrete_mark)
      cat('\nMarks are discrete.')
    else
      cat('\nMarks are continuous.')
    cat('\nPossible mark values:',1:nmarks)
    cat('\nFrequency table of the point pattern marks:\n')
    print(table(ppmarks))
  }
  else
  {
    pp=MPPfit$pp
    # test if this pattern has marks
    if(length(pp$marks) == 0) {
      stop("This point pattern doesn't have any marks.")
    }
    # get the levels of the mark
    ppmarks <- pp$marks
    if (is.factor(pp$marks)) {
      marks <- as.numeric(levels(pp$marks))
    } else {
      marks <- sort(unique(pp$marks))
    }
    nmarks <- length(marks)
    cat('\nNumber of marks is ',nmarks)
    cat('\nMark distribution: ',MPPfit$mark_dist,'\n')
    # get number of points and the location of points
    summary(pp)
    cat('\nEvent allocation to mark values:')
    print(table(pp$marks))
  }
}

#' Generate a Marked Poisson point process (conditional on mark)
#'
#' @description
#' This function generates realizations (point patterns) from
#' a given Marked IPPP via conditioning of the joint intensity surface
#' on its marked component. See details for the choice of models
#' for the mark distribution. For each mark value we obtain a ground
#' process. There processes are
#' standard IPPP (unmarked) with mixture intensity surfaces. The mark
#' distribution is responsible for the number of events in the point pattern.
#'
#' For examples see
#'
#' \url{http://faculty.missouri.edu/~micheasa/sppmix/sppmix_all_examples.html
#' #rMIPPP_cond_mark}
#'
#' @param lambda Average number of mark values
#' observed over the window. This is the
#' total number of points observed (on the average).
#' @param params Parameters for the mark
#' distribution. The value depends on the \code{mark_distr_choice}
#' parameter, e.g., \code{params} is a vector
#' of probabilities if the mark
#' distribution is discrete
#' (\code{mark_distr_choice=0}).
#' @param mark_distr_choice A number indicating which
#' mark distribution to use. In the
#' discrete mark case, the mark distribution
#' is discrete over the marks \code{1:length(params)} with
#' corresponding probabilities in \code{params}.
#' The continuous mark case has not been implemented yet.
#' @param truncate Logical variable indicating whether or not we
#' normalize the densities of the mixture components
#' to have all their mass within the window defined
#' in the window \code{win}. This affects the mixture model for the
#' intensity surface of the ground process.
#' @param discrete_mark Logical flag indicating whether the mark is discrete or not.
#' Default is TRUE. For continuous marks set this to FALSE.
#' @param win Object of type \code{\link[spatstat]{owin}} defining the window of observation.
#' @param bigwin Object of type \code{\link[spatstat]{owin}}. If supplied, this will be the
#' window of observation, even if the pattern is generated over \code{win}. Useful if we
#' do not truncate (\code{truncate=FALSE}) and we want better presentation of the generated MIPPP.
#' @param open_new_window Open a new window for a plot.
#' @param grayscale Logical to request plots in grayscale.
#' @param show_plots Logical variable requesting to produce exploratory plots of the
#' Marked IPPP intensity surface and generated point pattern for each mark.
#' @details For discrete marks, we assume that the joint intensity function of a
#' marked point pattern \code{N=[s,m]} with \code{n}
#' events is of the form:
#'
#' \code{intensity(s,m)=lambda*M(m|theta1)*g(s(m)|theta2(m))}
#'
#' where \code{m} denotes a mark and \code{s=s(m)}
#' a location with mark \code{m}, lambda a parameter
#' with the interpretation as the average number of events
#' over the window of observation, and \code{M} the mark distribution and
#' \code{g} the ground intensity are proper densities.
#'
#' In order to simulate from this Marked IPPP
#' we first simulate the number of events
#' and their marks from an IPPP with
#' intensity \code{lambda*M(m|theta1)}, and then generate
#' the ground intensities for each mark. Marks are assumed to be
#' independnet of each other and the mixture parameters
#' describing each ground process are also assumed to be independent
#' of each other.
#'
#' The continuous mark case will be implemented in future releases.
#'
#' See Micheas (2014) for more details on
#' Marked IPPP models via conditioning arguments.
#' @references Hierarchical Bayesian Modeling of Marked Non-Homogeneous Poisson Processes with finite mixtures and inclusion of covariate information. Micheas, A.C. (2014). Journal of Applied Statistics, 41, 12, 2596-2615, DOI: 10.1080/02664763.2014.922167.
#' @return A list containing the following components:
#' \item{groundsurfs}{A list of \code{intensity_surface} objects containing the surfaces of the ground processes (one for each discrete mark value).}
#' \item{groundPPs}{A list of \code{\link[spatstat]{ppp}} objects containing the locations of the ground processes (one for each discrete mark value).}
#' \item{genMPP}{The generated point pattern as an object of class \code{\link[spatstat]{ppp}} and \code{sppmix}. The member \code{$marks} contains the marks at each of the generated locations.}
#' \item{mark_distr_choice}{The choice of mark distribution. Same as the supplied parameter.}
#' \item{params}{The default or supplied parameter \code{params}.}
#' @author Sakis Micheas
#' @seealso \code{\link{plotmix_2d}}
#' @examples
#' \donttest{
#' # Create a marked point pattern; use randomization and 2 discrete uniform
#' # marks (default values)
#' newMPP=rMIPPP_cond_mark(bigwin = spatstat::owin(c(-10,10),c(-10,10)))
#' newMPP$params
#' plot(newMPP$genMPP, showmarks=TRUE)+add_title("Marked Poisson point pattern",
#'  n=newMPP$genMPP$n, nmarks=2)
#' plotmix_2d(newMPP$groundsurfs[[1]], newMPP$groundPPs[[1]])+ add_title(
#'  "Poisson point pattern for mark 1", n=newMPP$genMPP$n, m=newMPP$groundsurfs[[1]]$m)
#' plotmix_2d(newMPP$groundsurfs[[2]], newMPP$groundPPs[[2]])+ add_title(
#'  "Poisson point pattern for mark 2", n=newMPP$genMPP$n, m=newMPP$groundsurfs[[2]]$m)}
#'
#' @export
rMIPPP_cond_mark<- function(
      lambda=500,params=c(.5,.5),
      mark_distr_choice=0,
      truncate=FALSE,discrete_mark = TRUE,
      win=owin(c(-3,3),c(-3,3)),
      bigwin,open_new_window = FALSE,
      grayscale=FALSE,show_plots = TRUE)
{
  xlims=win$xrange
  ylims=win$yrange
  #generate the marks first
  #depending on the mark distribution
  #params has different meaning
  if(discrete_mark)
  #discrete uniform 1:length(params)
  {
    if(sum(params)!=1)
      stop("Probabilities do not sum up to 1")
    nmarks=length(params)
    cat("\nMark distribution is discrete over",nmarks,"marks.\n")
    cat("Probabilities:",params,"\n")
    marks=1:nmarks
    n=rpois(1,lambda=lambda)
    gen_marks<-sample(x=1:nmarks,
            size=n, replace=TRUE,
            prob=params)
    countmarks=table(gen_marks)
    cat("Mark allocation for",n,"generated points:\n")
    print(table(gen_marks))
    #generate the ground process for this
    #mark value, its just a normal mixture
    groundsurfs=vector("list",nmarks)
    groundPPs=vector("list",nmarks)
    pat=NULL
    patmarks=NULL
    for(i in 1:nmarks)
    {
      truemix<-rnormmix(m=5,df=10,
        rand_m = TRUE,xlim=xlims,
        ylim=ylims)
      groundsurfs[[i]]=to_int_surf(truemix,
          lambda = countmarks[i],
          win=win)
      genpts=gen_n_from_mix(
        countmarks[i],truemix)
      pat=rbind(pat,genpts)
      patmarks=c(patmarks,
                 rep(i,countmarks[i]))
      groundPPs[[i]] <-ppp(x=genpts[,1],
          y=genpts[,2],window=win,check= truncate)
      groundPPs[[i]]$comp <- genpts[, 3]
      class(groundPPs[[i]]) <- c("sppmix", "ppp")
      if(!missing(bigwin))
      {
        groundsurfs[[i]]$window=bigwin
        groundPPs[[i]]$window=bigwin
      }
    }
    if(!missing(bigwin))
    {
      win=bigwin
    }
    genMPP <- ppp(x=pat[,1],y=pat[,2],
                  window=win,marks=patmarks,
                  check= truncate)
    genMPP$comp <- pat[, 3]
    class(genMPP) <- c("sppmix", "ppp")

    RVAL <- list(groundsurfs = groundsurfs,
                 groundPPs=groundPPs,
                 genMPP=genMPP,
                 mark_distr_choice=mark_distr_choice,
                 params=params)
    if(show_plots)
    {
      print(plot(genMPP,open_new_window=open_new_window,showmarks=TRUE)+
        add_title("Marked Poisson point pattern",n=genMPP$n,nmarks=nmarks))
      for(i in 1:nmarks)
      {
        print(plotmix_2d(groundsurfs[[i]],groundPPs[[i]],
                   open_new_window=open_new_window,grayscale = grayscale)+
          ggtitle(paste("Poisson point pattern for mark",i,
           "\n",countmarks[i],"points,",groundsurfs[[i]]$m,"components")))
      }
    }
  }
  else
    stop("Option not implemented.")

  return(RVAL)
}

#' Fit a MIPPP conditionally on mark
#'
#' @description
#' This function fits a Marked IPPP (MIPPP) on a marked
#' point pattern by modeling the (joint)
#' intensity surface of the locations and the marks
#' using an IPPP for the marks (independent
#' of the locations) and an IPPP with mixture intensity
#' for the corresponding ground process, where the
#' mixture parameters depend on the mark value.
#' NOTE: The estimation procedure for continuous
#' marks will be implemented
#' in future versions of the \code{sppmix} package.
#'
#' For examples see
#'
#' \url{http://faculty.missouri.edu/~micheasa/sppmix/sppmix_all_examples.html
#' #est_MIPPP_cond_mark}
#'
#' @param pp Marked point pattern of class \code{\link[spatstat]{ppp}}.
#' @param m A vector representing the number
#' of components to fit for the ground
#' process corresponding to each mark. Since in
#' real applications we don't know these numbers
#' we can specify an integer so that
#' the routine will fit a BDMCMC
#' with this \code{m} as the maximum
#' number of components. Then we use
#' the MAP number of components for
#' each ground process with a mixture intensity function
#' of this many components. If not supplied
#' the default is \code{m=10}.
#' @param L Number of iterations for the MCMC; default is 50000.
#' @param burnin Number of initial realizations to discard. By default, it is 1/10 of the total number of iterations.
#' @param hyper_da A list of hyperparameters for
#' \code{\link{est_mix_damcmc}}. Each element of
#' this list should contain 3 values (hyperparameters) and
#' the number of elements should be the same as the
#' number of marks. If this parameter is omitted
#' the default hyperparameters of \code{\link{est_mix_damcmc}} will be used.
#' @param hyper Hyperparameter for the mark distribution. Must be a vector of positive real numbers. If omitted the vector of one's is used.
#' @param fit_markdist Logical variable requesting to fit and return the parameter estimates of the mark distribution.
#' @param truncate Logical variable indicating whether or not we
#' we only work with events within the window defined
#' in the point pattern \code{pp}.
#' @param grayscale Logical to request plots in grayscale.
#' @param discrete_mark Logical flag indicating whether the mark is discrete or not.
#' For continuous marks set this to FALSE.
#' @param LL Length of the side of the square grid.
#' @param open_new_window Open a new window for a plot.
#' @param show_plots Logical variable requesting to produce the ground fits and probability field plots for each mark. If label switching is present, the MAPE surface is computed and returned, otherwise the PME.
#' @param compute_surfaces Logical to request computation of the Average of Surfaces (if \code{m} is a vector) or the Bayesian Model Average (if \code{m} is an integer or missing). Default is TRUE. This is a SLOW operation.
#' @return An object of class \code{MIPPP_fit}, which is simply a list containing the following components:
#' \item{gen_mark_ps}{The posterior realizations of the discrete mark distribution probabilities.}
#' \item{mark_dist}{The posterior means of the discrete mark distribution probabilities.}
#' \item{discrete_mark}{Same logical flag as the input argument.}
#' \item{pp}{Same as the input argument.}
#' \item{ground_fits}{A List of objects of type \code{damcmc_res} which contain the results of the DAMCMC (or the BDMCMC for MAP number of components) fits to the ground process for each discrete mark value.}
#' \item{ground_fitsAoS}{A List of objects of type \code{\link[spatstat]{im}} which contain the AoS (average of surfaces) surface based on the DAMCMC (or the BMA from BDMCMC) fits to the ground process for each discrete mark value.}
#' \item{post_surf}{A List of \code{intensity_surface} objects, one for each mark, representing the surface of posterior means, after fixing label switching using SEL permutation.}
#' \item{condition_on_loc}{Logical variable indicating the type of conditioning used in order to produce this MIPPP fit. For this function it is set to FALSE.}
#' \item{fit_DAMCMC}{Logical variable indicating whether or not a DAMCMC or BDMCMC fit was requested.}
#' \item{m}{Same as input.}
#' @author Sakis Micheas, Jiaxun Chen
#' @seealso \code{\link{rMIPPP_cond_mark}},
#' \code{\link{GetStats}}
#' @examples
#' \donttest{
#' #Create a marked point pattern; use randomization and 3 discrete marks
#' newMPP=rMIPPP_cond_mark( params=c(.2,.5,.3),bigwin = spatstat::owin(c(-10,10),c(-10,10)))
#' newMPP$params
#' #supply the true number of components for each ground process
#' m=c(newMPP$groundsurfs[[1]]$m, newMPP$groundsurfs[[2]]$m, newMPP$groundsurfs[[3]]$m)
#' MIPPPfit=est_MIPPP_cond_mark(newMPP$genMPP,m=m,compute_surfaces=FALSE)
#' #check out the mark distribution parameters
#' #posterior means
#' MIPPPfit$mark_dist
#' #credible sets
#' GetStats(MIPPPfit$gen_mark_ps[,1])$CredibleSet#should contain .2
#' GetStats(MIPPPfit$gen_mark_ps[,2])$CredibleSet#should contain .5
#' GetStats(MIPPPfit$gen_mark_ps[,3])$CredibleSet#should contain .3
#' #now pretend we do not know the truth as is usually the case. Supply an integer
#' #for m so that the routine will fit a BDMCMC with this as the max number of
#' #components and use the MAP number of components
#' MIPPPfit=est_MIPPP_cond_mark(newMPP$genMPP,m=7,compute_surfaces=FALSE)
#' #check out the mark distribution parameters
#' MIPPPfit$mark_dist
#' GetStats(MIPPPfit$gen_mark_ps[,1])$CredibleSet#should contain .2
#' GetStats(MIPPPfit$gen_mark_ps[,2])$CredibleSet#should contain .5
#' GetStats(MIPPPfit$gen_mark_ps[,3])$CredibleSet#should contain .3}
#'
#' @export
est_MIPPP_cond_mark <- function(pp,m=10, L = 50000,
                                burnin = floor(L/10), hyper_da, hyper,
                                fit_markdist=TRUE,
                                truncate = FALSE,grayscale = FALSE,
                                discrete_mark = TRUE, LL = 256,
                                open_new_window = FALSE, show_plots = TRUE,
                                compute_surfaces = TRUE)
{
  if(!discrete_mark) {
    stop("Option not implemented yet.")
  }
  # test if this pattern has marks
  if(length(pp$marks) == 0) {
    stop("This point pattern doesn't have any marks.")
  }
  # get the levels of the mark
  if(is.factor(pp$marks)) {
    marks <- as.numeric(levels(pp$marks))
  } else {
    marks <- sort(unique(pp$marks))
  }

  if(is.vector(m,mode="integer"))
  {
    if(any(m<=0))
      stop("Bad m parameters.")
    else
      fit_DAMCMC=TRUE
  }
  else if (m>0)
  {
    m=floor(m)
    fit_DAMCMC=FALSE
  }
  else
    stop("Bad m parameter.")
  nmarks <- length(marks)
  fitsDA <- list()
  fitsAoS <- list()
  post_surf <- list()
  if(missing(hyper_da))
  {
    hyper_da <- list()
    for (k in 1:nmarks)
    {
      if(fit_DAMCMC)
        hyper_da[[k]] <- c(3, 1, 1)
      else
        hyper_da[[k]] <- c(m/2, 1/36, 3, 2, 1, 1)
    }
  }

  if(missing(hyper))
  {
    hyper <- rep(1, nmarks)
  }
  else if (any(hyper<=0))
  {
    hyper <- rep(1, nmarks)
  }

  for(k in 1:nmarks)
  {
    pattern <- pp[pp$marks == marks[k]]
    if(fit_DAMCMC)
    {
      fitsDA[[k]]<-est_mix_damcmc(pattern, m[k], hyper_da = hyper_da[[k]], L =L,
                                  truncate = truncate)
      fitsDA[[k]]=drop_realization(fitsDA[[k]],drop=burnin)
      if(compute_surfaces)
        fitsAoS[[k]]=plot_avgsurf(
          fitsDA[[k]],win = fitsDA[[k]]$data$window,
          LL = 100,burnin =0,zlims = c(0, 0),
          grayscale = grayscale,showplot = FALSE)
      else
        fitsAoS[[k]]=NULL
    }
    else
    {
      BDMCMCfit=est_mix_bdmcmc(pp=pattern,m=m,L=L,
                                 hyper = hyper_da[[k]],truncate = truncate)
      if(compute_surfaces)
        fitsAoS[[k]]=GetBMA(BDMCMCfit,LL=100)
      else
        fitsAoS[[k]]=NULL
      BDMCMCfit=drop_realization(BDMCMCfit,drop=burnin)
      #now we drop the bad realizations
      BDMCMCfit=drop_realization(BDMCMCfit,(BDMCMCfit$Badgen==1))
      cat("\nDropped bad BDMCMC realizations. Number left:",BDMCMCfit$L,"realizations\n")
      BDtab=GetBDTable(BDMCMCfit,TRUE)
      MAPm=BDtab$MAPcomp
      cat("\nThe ground process for mark value",marks[k],"has",MAPm," (MAP) number of components\n")
      BDMCMCfitMAPcomp=GetBDCompfit(BDMCMCfit,MAPm,burnin=0)
      fitsDA[[k]]=BDMCMCfitMAPcomp$BDgens
    }

    if(show_plots)
    {
      test_label <- check_labels(fitsDA[[k]], showmessages = FALSE,lagnum = 10)
      if(test_label)
      {
        cat("\nLabel switching present for mark",k,
            "(computing the MAPE surface, please wait...)\n")#(permuting labels, please wait...)\n")
        #fixed_real <- FixLS_da(fitsDA[[k]], burnin = 0,approx = FALSE, run_silent = TRUE)
        #post_surf[[k]] <- GetPMEst(fixed_real, burnin = 0)
        post_surf[[k]] <- GetMAPEst(fitsDA[[k]], burnin = 0)
      } else
      {
        cat("\nNo label switching detected for mark",k,"\n")
        post_surf[[k]] <- GetPMEst(fitsDA[[k]], burnin = 0)
      }
      print(plotmix_2d(post_surf[[k]], pattern, open_new_window = open_new_window,
                       grayscale = grayscale, L = LL)+
              ggtitle(paste("PME of the IPPP intensity surface for mark",
                            marks[k],"(SEL)\nm =",post_surf[[k]]$m,", n =", pattern$n)))
#      MAP_surf <- GetMAPEst(fitsDA[[k]], burnin = 0)
#      print(plotmix_2d(MAP_surf, pattern, open_new_window = open_new_window,
#                       grayscale = grayscale, L = LL)+
#              add_title("IPPP intensity surface of MAP estimates",
#                        lambda = MAP_surf$lambda, m = MAP_surf$m, n = pattern$n))
    }
  }
  if(fit_markdist)
  {
    markcounts <- as.numeric(table(pp$marks))
    post_param <- markcounts + hyper
    gen_markps <- matrix(0, nmarks, L)
    for (i in 1:L)
    {
      gen_markps[, i] <- rDirichlet_sppmix(post_param)
    }
    gen_mark_ps <- t(gen_markps)
    post_temp <- colMeans(gen_mark_ps)
    post_mark_ps <- post_temp/sum(post_temp)
  }
  else
  {
    gen_mark_ps <- NULL
    post_mark_ps <- NULL
  }
  RVAL <- list(gen_mark_ps = gen_mark_ps,
               mark_dist = post_mark_ps,
               discrete_mark = discrete_mark,
               pp = pp,
               ground_fits = fitsDA,
               ground_fitsAoS = fitsAoS,
               post_surf=post_surf,
               condition_on_loc=FALSE,
               fit_DAMCMC=fit_DAMCMC,
               m=m)
  class(RVAL) <- c("MIPPP_fit")
  return(RVAL)
}
