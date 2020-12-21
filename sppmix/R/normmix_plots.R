#' @title Plots a normal mixture intensity in 3d
#' @description
#' Plot the 3d intensity surface of a Poisson point process
#' with mixture intensity of normal components.
#'
#' For examples see
#'
#' \url{http://faculty.missouri.edu/~micheasa/sppmix/sppmix_all_examples.html
#' #plot.intensity_surface}
#'
#' @param x Object of class \code{intensity_surface} or \code{normmix}.
#' @param truncate Requests to truncate the components
#' of the mixture intensity to have all their mass
#' within the window of the intensity object intsurf. Default is TRUE.
#' @param L Length of the side of the square grid.
#' The intensity is calculated on an L * L grid.
#' The larger this value is, the better the picture resolution.
#' @param zlims The limits of the z axis. Defaults to [0,1.1*max(intensity)].
#' @param main Title for the plot.
#' @param grayscale Logical flag to request a gray scale plot.
#' @param ... Additional parameters passed to \code{to_int_surf()}.
#'
#' @author Jiaxun Chen, Sakis Micheas, Yuchen Wang
#' @seealso \code{\link{normmix}},
#' \code{\link{to_int_surf}}
#' @examples
#' \donttest{
#' truemix <- rnormmix(m = 5, sig0 = .1, df = 5, xlim= c(-1, 5), ylim =c(2, 5))
#' intsurf=to_int_surf(truemix, lambda = 200, win =spatstat::owin( c(-1, 5),c(2, 5)))
#' plot(intsurf,main = "True Poisson intensity surface (mixture of normal components)")
#' #use the demo intensity surface
#' demo_intsurf
#' summary(demo_intsurf)
#' #3d plot of the intensity surface
#' plot(demo_intsurf,main = "True Poisson intensity surface (mixture of normal components)")}
#'
#' @export
#' @method plot intensity_surface
plot.intensity_surface <- function(x, truncate = TRUE, L = 256,
                                   zlims = c(0, 0),
                                   main = "Poisson intensity surface (mixture of normal components)",
                                   grayscale = FALSE, ...) {

  intsurf <- to_int_surf(x, ...)
#  cat("\nTo save all open rgl graphs use Save_AllOpenRglGraphs.\n")

  win <- intsurf$window
  xlims <- win$xrange
  ylims <- win$yrange

  est_intensity <- dnormmix(intsurf, xlim = xlims, ylim = ylims,
                            L = L, truncate = truncate)

  x <- est_intensity$xcol
  y <- est_intensity$yrow
  z <- t(est_intensity$v)

  jet.colors <- colorRampPalette(c("#00007F", "blue", "#007FFF", "cyan",
                                   "#7FFF7F", "yellow", "#FF7F00", "red",
                                   "#7F0000"))
  if (grayscale == TRUE) {
    col <- gray.colors(100,start = 1, end = 0)[findInterval(z, seq(min(z),
                                                   max(z), length = 100))]
  } else {
    col <- jet.colors(100)[findInterval(z, seq(min(z),
                                                  max(z), length = 100))]
  }

  if (zlims[1] == 0 && zlims[2] == 0) {
    zlims=c(0,1.1*max(z))
  }

  rgl::open3d(windowRect=c(50,50,1000,800),
              zoom=1.2)
  U=rgl::par3d("userMatrix")
  rgl::par3d(userMatrix=
               rgl::rotate3d(U,pi/4,0,0,1))
  rgl::persp3d(x = x, y = y, z = z,
               color = col, xlab="x",ylab="y",zlab="",
               zlim=c(zlims[1]-0.01,zlims[2]),
               box = FALSE, axes = FALSE)
  rgl::axis3d('x')
  rgl::axis3d('y')
  rgl::axis3d('z-+', pos = c(xlims[1], ylims[2], 0))
  rgl::rgl.lines(c(xlims[1], xlims[1]),
                 c(ylims[2], ylims[2]),
                 c(0,max(z)),
                 color = 'black')
  rgl::title3d(main = NULL)
  rgl::text3d(xlims[2], ylims[2], zlims[2],
              texts = main)

  if (grayscale == TRUE) {
    rgl::bgplot3d(suppressWarnings(
      fields::image.plot(legend.only = TRUE,
                         zlim = zlims,
                         smallplot= c(.8,.82,0.05,.7),
                         col = gray.colors(100,start = 1, end = 0))))
  } else {
    rgl::bgplot3d(suppressWarnings(
      fields::image.plot(legend.only = TRUE,
                         zlim = zlims,
                         smallplot= c(.8,.82,0.05,.7),
                         col = jet.colors(100))))

  }
}


#' Plot a spatial point pattern
#'
#' @description
#' Plot a spatial point pattern generated from a Poisson with
#' mixture intensity surface. Alternatively, the function can
#' plot a spatstat \code{\link[spatstat]{ppp}} object.
#'
#' For examples see
#'
#' \url{http://faculty.missouri.edu/~micheasa/sppmix/sppmix_all_examples.html
#' #plot.sppmix}
#'
#' @param x A point pattern of class \code{\link{sppmix}} or
#' \code{\link[spatstat]{ppp}}.
#' @param mus An optional list of the theoretical means of the mixture components.
#' @param estcomp The estimated component label should be a vector whose length
#' should be the same as number of points. If \code{estcomp} is not missing, the function will plot the
#' points using different colors according to \code{estcomp}. See the example on how
#' to calculate \code{estcomp} from a DAMCMC fit. If this variable is missing and we
#' pass a point pattern generated using \code{rsppmix}, then the true component labels will be used, otherwise,
#' the function will not plot the points with different colors to indicate the different components.
#' @param open_new_window Open a new window for the plot.
#' @param colors Logical flag requesting to use different colors for the points based
#' on which component they belong to.
#' @param showmarks Logical flag requesting
#' to plot each point with a different circle size
#' according to its mark value. If the mark is a data.frame object (multivariate marks), the first column (mark) is used as the the marks and displayed.
#' @param whichmark If multivariate marks, choose to display this one.
#' @param ... Additional parameters to the \code{add_title} function. Valid choices
#' are m, n and L. To add a different title than the default, use \code{add_title} after the plot call (see examples below).
#' @author Jiaxun Chen, Sakis Micheas, Yuchen Wang
#' @seealso \code{\link{normmix}},
#' \code{\link{to_int_surf}},
#' \code{\link[spatstat]{owin}},
#' \code{\link{rsppmix}},
#' \code{\link{est_mix_damcmc}},
#' \code{\link{GetMAPLabels}},
#' \code{\link{rMIPPP_cond_mark}}
#' @examples
#' \donttest{
#' mix1 <- rnormmix(5, sig0 = .01, df = 5, xlim=c(0, 5), ylim=c(0, 5))
#' intsurf1=to_int_surf(mix1, lambda = 40, win =spatstat::owin( c(0, 5),c(0, 5)))
#' pp1 <- rsppmix(intsurf1)
#' plot(pp1)
#' plot(pp1, mus=intsurf1$mus)
#' plot(pp1,mus=intsurf1$mus)+add_title(
#'  "Poisson point pattern along with the true component means", m=intsurf1$m,n=pp1$n)
#' plot(pp1, mus = intsurf1$mus, lambda = intsurf1$lambda)
#' plot(pp1, mus = intsurf1$mus)+ add_title(
#'  "Poisson point pattern along with the true component means", lambda = intsurf1$lambda,
#'  m=intsurf1$m,n=pp1$n)
#' #use the demo intensity surface
#' demo_intsurf
#' pp2 <- rsppmix(demo_intsurf,marks = 1:3)
#' plot(pp2)
#' plot(pp2, mus = demo_intsurf$mus)#plot the mixture means as well
#' #plot the points with different colors depending on the true component label
#' plot(pp2, colors = TRUE)
#' #plot the points with different colors depending on the estimated component label
#' fit <- est_mix_damcmc(pp2, 2)
#' est_comp <- GetMAPLabels(fit)
#' plot(pp2, estcomp = est_comp, colors = TRUE)
#' #generate and plot a marked point pattern
#' newMPP=rMIPPP_cond_mark()
#' plot(newMPP$genMPP, showmarks=TRUE)}
#'
#' @export
#' @method plot sppmix
plot.sppmix <- function(x, mus, estcomp,
      open_new_window=FALSE, colors = FALSE,showmarks=TRUE,whichmark=1, ...)
{
  pattern<-x
  if(pattern$markformat=="none")
    showmarks=FALSE
  else
  {
    if(class(pattern$marks)=="data.frame")
      pattern$marks=pattern$marks[,whichmark]
  }

  y=NULL
  n <- pattern$n
  openwin_sppmix(check2open=open_new_window)
  if (colors == TRUE) {
    if (!missing(estcomp)) {
      comp <- factor(estcomp)
    } else {
      if (length(pattern$comp)>0) {
        comp <- factor(pattern$comp)
      } else {
        colors =FALSE#stop("Missing true component labels or estimated component labels.")
      }
    }
    if(colors == TRUE)
      p <- ggplot2::ggplot(as.data.frame(pattern),
                           aes(x, y, color = comp))
    else
      p <- ggplot2::ggplot(as.data.frame(pattern),
                           aes(x, y))

    if(!showmarks)
      p<-p+ geom_point()
    p<-p+ ggplot2::labs(x = "x", y = "y") + scale_color_hue(
        name="Component") +
      ggplot2::coord_cartesian(xlim = pattern$window$xrange,
                               ylim =pattern$window$yrange,expand=FALSE)+
      ggplot2::theme_classic() +
      ggplot2::theme(panel.border = ggplot2::element_rect(fill = NA, size = 1)) +
      #    ggplot2::limits(pattern$window$xrange,"x")+
      #    ggplot2::limits(pattern$window$yrange,"y")+
      add_title("Point pattern from a Poisson with mixture intensity", n = pattern$n, ...)

  } else {
  p <- ggplot2::ggplot(as.data.frame(pattern),
        aes(x, y))
  if(!showmarks)
    p<-p+ geom_point()
  p<-p+
    ggplot2::labs(x = "x", y = "y") +
    ggplot2::coord_cartesian(xlim = pattern$window$xrange,
                         ylim =pattern$window$yrange,expand=FALSE)+
   ggplot2::theme_classic() +
    ggplot2::theme(panel.border = ggplot2::element_rect(fill = NA, size = 1)) +
#    ggplot2::limits(pattern$window$xrange,"x")+
#    ggplot2::limits(pattern$window$yrange,"y")+
    add_title("Point pattern from a Poisson with mixture intensity", n = pattern$n, ...)
  }

  if(showmarks)
      p<-p+geom_point(data=as.data.frame(pattern),aes(x=x,y=y,size=marks),shape=21)+
    scale_size_continuous(breaks=sort(unique(pattern$marks)))+
    guides(size =guide_legend(title="Mark",ncol=2,byrow=TRUE))

  if (!missing(mus)) {
    mean_df <- data.frame(do.call(rbind, mus))
    names(mean_df) <- c("x", "y")
    p <- p + ggplot2::geom_point(data = mean_df, color = "red",shape = "x", size = 8) +
      add_title("Point pattern from a Poisson with mixture intensity", n = pattern$n,
                m = nrow(mean_df), ...)
  }
  p
}

#' Plot a spatial point pattern
#'
#' @description
#' Standard 2d plot for a spatial point pattern.
#'
#' For examples see
#'
#' \url{http://faculty.missouri.edu/~micheasa/sppmix/sppmix_all_examples.html
#' #plot2dPP}
#'
#' @param pp A point pattern of class sppmix or
#' \code{\link[spatstat]{ppp}}.
#' @param mus An optional list of the theoretical means of the mixture components.
#' @param add2plot Logical variable to indicate if the function should add the points to an existing plot.
#' @param title1 Title for the plot.
#' @param open_new_window Open a new window for the plot.
#' @author Sakis Micheas
#' @seealso \code{\link{normmix}},
#' \code{\link{to_int_surf}},
#' \code{\link[spatstat]{owin}},
#' \code{\link{rsppmix}}
#' @examples
#' \donttest{
#' mix1 <- rnormmix(5, sig0 = .01, df = 5, xlim=c(0, 5), ylim=c(0, 5))
#' intsurf1=to_int_surf(mix1, lambda = 40, win =spatstat::owin( c(0, 5),c(0, 5)))
#' pp1 <- rsppmix(intsurf1)
#' plot2dPP(pp1)
#' plot2dPP(pp1, mus = intsurf1$mus)}
#'
#' @export
plot2dPP <- function(pp, mus,add2plot=FALSE,
    title1="Spatial point pattern",
    open_new_window=FALSE ) {
  n =pp$n
  if(add2plot)
    points(pp$x,pp$y,cex=0.8,pch=20)
  else
  {
    openwin_sppmix(check2open=open_new_window)

    plot(pp$x,pp$y,pch=20,cex=0.8,main="",xlab="x",
       ylab="y")
  }
  if(missing(mus))
  {
    titleLines <- list(
    bquote(paste("n=",.(n)," points")),
      title1)
    mtext(do.call(expression, titleLines),side=3,line=0:1)
  }
  else
  {
    m=length(mus)
    for(j in 1:m)
    {
      points(mus[[j]][1],mus[[j]][2],cex=1.5,pch="x",col="red")
    }
    titleLines <- list(
      bquote(paste("n=",.(n)," points, m=",.(m)," components")),
      title1)
    mtext(do.call(expression, titleLines),side=3,line=0:1)
  }
}

#' 2d exploratory plots for mixture intensity surfaces
#'
#' @description
#' Create a 2d image or contour plot of the
#' intensity surface, with the option to display
#' a point pattern.
#'
#' For examples see
#'
#' \url{http://faculty.missouri.edu/~micheasa/sppmix/sppmix_all_examples.html
#' #plotmix_2d}
#'
#' @inheritParams rsppmix
#' @param pattern Optional spatial point pattern
#' to add to the plot. This is an object of
#' class \code{\link[spatstat]{ppp}}.
#' @param estcomp The estimated component label should be a vector whose length
#' should be the same as number of points. If \code{estcomp} is not missing, the function will plot the
#' points using different colors according to \code{estcomp}. See the example on how
#' to calculate \code{estcomp} from a DAMCMC fit. If this variable is missing and we
#' pass a point pattern generated using \code{rsppmix}, then the true component labels will be used, otherwise,
#' the function will not plot the points with different colors to indicate the different components.
#' @param contour Logical flag requesting the countour plot only.
#' @param L Length of the side of the square grid.
#' The intensity is calculated on an L * L grid.
#' The larger this value is, the better the picture resolution.
#' @param open_new_window Open a new window for the plot.
#' @param grayscale Plot in gray scale. Default is FALSE (use colors).
#' @param colors Logical flag requesting to use different colors for the points based
#' on which component they belong to.
#' @param ... Additional parameters passed to \code{to_int_surf()}.
#' @import ggplot2
#' @seealso \code{\link{normmix}},
#' \code{\link{to_int_surf}},
#' \code{\link[spatstat]{owin}},
#' \code{\link{rsppmix}},
#' \code{\link{GetMAPLabels}},
#' \code{\link{est_mix_damcmc}},
#' \code{\link{PlotUSAStates}}
#' @author Jiaxun Chen, Sakis Micheas, Yuchen Wang
#' @examples
#' \donttest{
#' # plot normmix density
#' truemix<- rnormmix(m = 3, sig0 = .1, df = 5, xlim= c(0, 5), ylim = c(0, 5))
#' summary(truemix)
#' intsurf=to_int_surf(truemix, lambda = 100, win =spatstat::owin( c(0, 5),c(0, 5)))
#' #plot the intensity surface
#' plotmix_2d(intsurf)
#' plotmix_2d(intsurf,contour = TRUE)
#' pp1 <- rsppmix(intsurf = intsurf)# draw points
#' plotmix_2d(intsurf, pp1)
#' plotmix_2d(intsurf, pp1,contour = TRUE)
#' #fit a Poisson with mixture intensity surface
#' CAgens=est_mix_damcmc(pp = CAQuakes2014.RichterOver3.0, m = 5)
#' #retrieve the surface of posterior means
#' CAfit=GetPMEst(CAgens)
#' #plot the surface and the point pattern
#' plotmix_2d(CAfit,CAQuakes2014.RichterOver3.0)
#' #to include the state boundaries use function PlotUSAStates
#' ret=PlotUSAStates(states=c('California','Nevada','Arizona'), showcentroids=FALSE,
#'  shownames=TRUE, main="Earthquakes in CA, 2014", pp=CAQuakes2014.RichterOver3.0,
#'  surf=CAfit, boundarycolor="white", namescolor="white")
#' #plotting the points with different colors depending on the component they belong to
#' truemix <- rnormmix(m = 5, sig0 = .1, df = 5, xlim=c(-2,2), ylim=c(-2,2))
#' intsurf=to_int_surf(truemix, lambda = 100, win = spatstat::owin(c(-2,2),c(-2,2)))
#' pp1 <- rsppmix(intsurf)
#' #plot the points with different colors depending on the true component label
#' plotmix_2d(intsurf,pp1, colors = TRUE)
#' #plot the points with different colors depending on the estimated component label
#' fit <- est_mix_damcmc(pp1, 5)
#' est_comp <- GetMAPLabels(fit)
#' plotmix_2d(intsurf,pp1, estcomp = est_comp, colors = TRUE)
#' plotmix_2d(intsurf,pp1, estcomp = est_comp, contour = TRUE,colors = TRUE)}
#'
#' @export
plotmix_2d <- function(intsurf, pattern,estcomp, contour = FALSE, truncate = TRUE,
      L = 256,open_new_window=FALSE,grayscale = FALSE,
      colors = FALSE, ...) {
  x=y=value=..level..=NULL
  openwin_sppmix(check2open=open_new_window)

  intsurf <- to_int_surf(intsurf, ...)
  win <- intsurf$window

  est_intensity <- dnormmix(intsurf, xlim = win$xrange, ylim = win$yrange,
                            L = L, truncate = truncate)

#  p <- plot_density(as.data.frame(est_intensity), contour = contour,grayscale = grayscale) +
#    labs(fill = "Intensity")
  density_df=as.data.frame(est_intensity)
  xrange=range(density_df$x)
  yrange=range(density_df$y)
  p <- ggplot2::ggplot(density_df,
                       ggplot2::aes(x, y))+
    ggplot2::labs(x = "x", y = "y")+
    ggplot2::theme_classic() +
    ggplot2::theme(panel.border = ggplot2::element_rect(fill = NA, size = 1))

  if(grayscale)
    cols <- gray.colors(100,start = 1, end = 0)
  else
    cols <- c("#00007F", "blue", "#007FFF", "cyan", "#7FFF7F", "yellow",
               "#FF7F00", "red", "#7F0000")

  if (!contour) {
    p<-p + ggplot2::geom_raster(aes(fill = value), interpolate = TRUE) +
      ggplot2::scale_fill_gradientn(colors = cols) +
      ggplot2::guides(
        fill = guide_colorbar(title = "Intensity",nbin = 100,
          barheight = 15))+
      ggplot2::coord_cartesian(xlim =xrange,
                               ylim =yrange,expand=FALSE)
  } else {
    p<-p + ggplot2::stat_contour(aes(z = value, color = ..level..),bins=10) +
      ggplot2::scale_color_gradientn(colors = cols) +
      ggplot2::guides(color = guide_colorbar(
        title = "Elevation",nbin = 100, barheight = 15))+
      ggplot2::coord_cartesian(xlim =xrange,
                               ylim =yrange,expand=FALSE)
  }

  #show means
  mean_df <- data.frame(do.call(rbind, intsurf$mus))
  names(mean_df) <- c("x", "y")
  if(!contour)
   p<-p + ggplot2::geom_point(data = mean_df, color = "grey",shape = "x", size = 8)
  else
    p<-p + ggplot2::geom_point(data = mean_df, color = "red",shape = "x", size = 8)

  if (!missing(pattern))
  {
    if (colors == TRUE)
    {
      if (!missing(estcomp))
      {
        comp <- factor(estcomp)
      }
      else
      {
        if (length(pattern$comp)>0)
        {
          comp <- factor(pattern$comp)
        }
        else
        {
          colors =FALSE#stop("Missing true component labels or estimated component labels.")
         }
      }
      if(colors == TRUE)
        if(!contour)
        {
          p <-p+ ggplot2::geom_point(data=as.data.frame(pattern)
                ,aes(color =comp)
                )+guides(
                  color = guide_legend(
                    title="Component",nrow=2,
                    byrow=TRUE))
        }else{
          p <-p+ ggplot2::geom_point(data=as.data.frame(pattern)
                                     ,aes(shape =comp)
          )+guides(
            shape = guide_legend(
              title="Component",nrow=2,
              byrow=TRUE))

        }

        else
        p <-p+ ggplot2::geom_point(data=as.data.frame(pattern))

      p<-p+
        add_title(
          "Poisson Intensity Surface",
          lambda = intsurf$lambda,
          m = intsurf$m, n = pattern$n)
    }
    else
    {
      p<-p + geom_point(data=as.data.frame(pattern)) +
       add_title("Poisson Intensity Surface",
                lambda = intsurf$lambda, m = intsurf$m, n = pattern$n)
    }
  }
  else
  {
    p<-p + add_title("Poisson Intensity Surface",
                  lambda = intsurf$lambda, m = intsurf$m)
  }

  p
}

#' Plot a mixture of normal components
#'
#' @description
#' Create a 3d plot and 2d image or contour plots of the
#' density of a mixture of normal components.
#'
#' For examples see
#'
#' \url{http://faculty.missouri.edu/~micheasa/sppmix/sppmix_all_examples.html
#' #plot.normmix}
#'
#' @param x Object of class \code{normmix}.
#' @param xlim,ylim The observation window.
#' @param contour Logical flag requesting the countour plot only.
#' @param truncate Logical flag requesting that the components are truncated within the window.
#' @param open_new_window Open a new window for the plot.
#' @param grayscale Plot in gray scale. Default is FALSE (use colors).
#' @param L Length of the side of the square grid.
#' The intensity is calculated on an L * L grid.
#' The larger this value is, the better the picture resolution.
#' @param title1 Optional title for the 3d plot.
#' @param whichplots Requests plots of the normal mixture density (surface). To get only the 2d plot
#' set \code{whichplots}=0, only the 3d plot set \code{whichplots}=1,
#' or for both the 2d and 3d plots set \code{whichplots}=2. Default action is to produce both plots.
#' @param ... Additional arguments for the S3 method.
#' @author Jiaxun Chen, Sakis Micheas, Yuchen Wang
#' @seealso \code{\link{normmix}},
#' \code{\link{to_int_surf}},
#' \code{\link[spatstat]{owin}},
#' \code{\link{rsppmix}}
#' @examples
#' \donttest{
#' # plot normmix density
#' truemix<- rnormmix(m = 3, sig0 = .1, df = 5, xlim= c(-1, 2), ylim = c(-1, 2))
#' summary(truemix)
#' #plot the normal mixture
#' plot(truemix, xlim= c(-1, 2), ylim = c(-1,2),
#'  title1="True mixture density in 3d")+add_title(
#'  "True mixture of normals density")
#' plot(truemix,xlim= c(-1, 2), ylim = c(-1, 2),contour = TRUE)+add_title(
#'  "Contour plot of the true mixture of normals density")
#' #build a mixture intensity surface for the Poisson point process
#' trueintsurf=to_int_surf(truemix, lambda = 100, win=
#'  spatstat::owin( c(-1, 2),c(-1, 2)))
#' plot(trueintsurf)#plot the surface, it is lambda*normmix}
#'
#' @export
#' @rdname density_plots
#' @method plot normmix
plot.normmix <- function(x, xlim, ylim, contour = FALSE,
    truncate = FALSE,open_new_window=FALSE,
    grayscale=FALSE,L = 256,
    title1="Mixture with normal components",
    whichplots=2,...) {
  mix <- x
  stopifnot(is.normmix(mix))

  if (missing(xlim) || missing(ylim)) {

    retlims=GetMixtureLimitsList(mix=mix)
    if(missing(xlim))xlim=retlims[[1]]
    if(missing(ylim))ylim=retlims[[2]]
  }

  est_density <- dnormmix(mix, xlim = xlim, ylim = ylim,
                          L = L, truncate = truncate)

  if(whichplots>=1)
    plotmix_3d(est_density,title1=title1)

  if(whichplots!=1)
  {
    openwin_sppmix(check2open=open_new_window)
    p<-plot_density(as.data.frame(est_density),
     contour = contour, grayscale=grayscale) +
    labs(fill = "Density") +
    add_title("Density of a mixture with normal components", m = mix$m)
    #show means
    mean_df <- data.frame(do.call(rbind, mix$mus))
    names(mean_df) <- c("x", "y")
    if(!contour)
      p<-p + ggplot2::geom_point(data = mean_df, color = "grey",shape = "x", size = 8)
    else
      p<-p + ggplot2::geom_point(data = mean_df, color = "red",shape = "x", size = 8)
    p
  }
}

#' Plots a density or image
#'
#' @description
#' Create a 2d image or contour plot of the density, intensity or any
#' image object.
#'
#' For examples see
#'
#' \url{http://faculty.missouri.edu/~micheasa/sppmix/sppmix_all_examples.html
#' #plot_density}
#'
#' @param density_df A data frame. Typically density_df=as.data.frame(imdens), where imdens an \code{\link[spatstat]{im}} object.
#' @param contour Logical flag requesting the countour plot only.
#' @param grayscale Plot in gray scale. Default is FALSE (use colors).
#' @param main A title for the 2d plot.
#' @param pp Optional point pattern to display (a \code{\link[spatstat]{ppp}} or \code{sppmix} object).
#' @param surf Optional \code{intensity_surface} object containing means to be displayed in the plot.
#' @param ppsize Size of the points in the plot.
#' @details This function does not open a new window for the plot.
#' @seealso \code{\link{dnormmix}}
#' @author Sakis Micheas
#' @examples
#' \donttest{
#' # plot a mixture of normals density
#' truemix <- rnormmix(m = 3, sig0 = .1, df = 5, xlim= c(0, 5), ylim = c(0, 5))
#' summary(truemix)
#' normdens=dnormmix(truemix, xlim = c(0, 5), ylim = c(0, 5))
#' #2d plots
#' plot_density(normdens, main="2d mixture density plot\nWindow=[0,5]x[0,5]")
#' #Contour plot
#' plot_density(normdens, contour=TRUE, main="2d mixture contour plot\nWindow=[0,5]x[0,5]")}
#'
#' @export
plot_density <- function(density_df,
  contour = FALSE, grayscale = FALSE,
  pp=NULL,surf=NULL,ppsize=1.0,
  main="2d surface (density or intensity)")
{
  x=y=value=..level..=NULL
  density_df=as.data.frame(density_df)
  xrange=range(density_df$x)
  yrange=range(density_df$y)
  p <- ggplot2::ggplot(density_df,
   ggplot2::aes(x, y))+
    ggplot2::labs(x = "x", y = "y")+
    ggplot2::theme_classic() +
    ggplot2::theme(panel.border = ggplot2::element_rect(fill = NA, size = 1))

  if(grayscale)
    color <- gray.colors(100,start = 1, end = 0)
  else
  color <- c("#00007F", "blue", "#007FFF", "cyan", "#7FFF7F", "yellow",
             "#FF7F00", "red", "#7F0000")

  if (!contour) {
    p<-p + ggplot2::geom_raster(aes(fill = value), interpolate = TRUE) +
      ggplot2::scale_fill_gradientn(colors = color) +
      ggplot2::guides(fill = guide_colorbar(
        title = "Intensity",nbin = 100, barheight = 15))+
      ggplot2::coord_cartesian(xlim =xrange,
                               ylim =yrange,expand=FALSE)
  } else {
    p<-p + ggplot2::stat_contour(aes(z = value, color = ..level..)) +
      ggplot2::scale_color_gradientn(colors = color) +
      ggplot2::guides(color = guide_colorbar(
        title = "Elevation",nbin = 100, barheight = 15))+
      ggplot2::coord_cartesian(xlim =xrange,
                           ylim =yrange,expand=FALSE)
  }
  if(!is.null(pp))
  {
    #show the point pattern points
    pp_df <- data.frame(pp$x,pp$y)
    names(pp_df) <- c("x", "y")
    p<-p + ggplot2::geom_point(data = pp_df,size=ppsize)
  }
  if(!is.null(surf))
  {
    #show the true means
    mean_df <- data.frame(do.call(rbind, surf$mus))
    names(mean_df) <- c("x", "y")
    p<-p + ggplot2::geom_point(data = mean_df, color = "gray",shape = "x", size = 8)
  }
  p<-p+ggplot2::ggtitle(main)
  p
}

#' @title Add a title to an existing ggplot2 plot
#'
#' @description The function
#' adds extra titles to a ggplot2 plot.
#'
#' For examples see
#'
#' \url{http://faculty.missouri.edu/~micheasa/sppmix/sppmix_all_examples.html#add_title}
#'
#' @param title The title to use for the plot.
#' @param lambda,m,n,t,L,nmarks,mu,theta,nu,gamma,sigma,df Optional info
#' to display on the second row of
#' the title: average number of points,
#' number of components in the mixture,
#' number of points in the point pattern,
#' time frame, number of iterations, total number
#' of marks, a mean value, parameters
#' for a Matern covariance model, a parameter gamma, and
#' degrees of freedom, respectively.
#'
#' @author Sakis Micheas, Yuchen Wang
#' @seealso \code{\link{rnormmix}},
#' \code{\link{dnormmix}},
#' \code{\link{MaternCov}}
#' @examples
#' \donttest{
#' truemix = rnormmix(m = 5, sig0 = .1, df = 5,xlim= c(0, 3), ylim = c(0, 3))
#' intsurf=to_int_surf(truemix,lambda = 100,win=spatstat::owin( c(0, 5),c(0, 5)))
#' #plot the intensity surface
#' plotmix_2d(intsurf)+add_title("A pretty projection of the 3d surface to 2 dimensions")}
#'
#' @export
add_title <- function(title, lambda = "",
  m = "", n = "", t="", L = "",nmarks="",
  mu="",theta="",nu="",gamma="",sigma="",df="")
{
  if (!missing(lambda)) lambda <- bquote(paste(lambda == .(lambda)))
  if (!missing(m)) m <- bquote(paste(m == .(m), " components"))
  if (!missing(n)) n <- bquote(paste(n == .(n), " points"))
  if (!missing(t)) t <- bquote(paste("time ",.(t)))
  if (!missing(L)) L <- bquote(paste(L == .(L), " iterations"))
  if (!missing(nmarks)) nmarks <- bquote(paste(.(nmarks)," marks"))
  if (!missing(mu)) mu <- bquote(mu == .(mu))
  if (!missing(theta)) theta <- bquote(theta == .(theta))
  if (!missing(nu)) nu <- bquote(nu == .(nu))
  if (!missing(sigma)) sigma <- bquote(sigma == .(sigma))
  if (!missing(gamma)) gamma <- bquote(gamma == .(gamma))
  if (!missing(df)) df <- bquote(paste(.(df)," degrees of freedom"))

  all_char <- list(lambda = lambda, m = m, n = n, t=t, L = L,nmarks=nmarks,
                   mu =mu,theta=theta,nu=nu,gamma=gamma,sigma=sigma,df=df)
  non_empty_char <- all_char[nchar(all_char) > 0]

  cal <- do.call(function(...) substitute(list(...)), non_empty_char)

  ggplot2::ggtitle(substitute(atop(title, cal)))
}

#' Plot the density or intensity of a
#' normal mixture in 3d over a fine grid
#'
#' @description
#' When a \code{normmix} object is given, this
#' function calculates the mixture density over
#' a fine grid for the given window. When an \code{intensity_surface} object
#' is given, the function multiplies the density with
#'  the surface lambda, and returns the Poisson
#'  mixture intensity function over the grid. Used for plotting.
#'
#' For examples see
#'
#' \url{http://faculty.missouri.edu/~micheasa/sppmix/sppmix_all_examples.html
#' #plotmix_3d}
#'
#' @param dens_image An image as an object of class \code{\link[spatstat]{im}}.
#' @param title1 A title for the 3d plot.
#' @param zlims The limits of the z axis. Defaults to [0,1.1*max(dens_image)].
#' @param grayscale Plot in gray scale. Default is FALSE (use colors).
#'
#' @author Jiaxun Chen, Sakis Micheas, Yuchen Wang
#' @seealso \code{\link{rnormmix}},
#' \code{\link{dnormmix}}
#' @examples
#' \donttest{
#' truemix <- rnormmix(m = 5, sig0 = .1, df = 5, xlim= c(0, 3), ylim = c(0, 3))
#' normdens=dnormmix(truemix, xlim= c(0, 3), ylim = c(0, 3))
#' plotmix_3d(normdens)
#' plotmix_3d(normdens, title1="Density of a normal mixture")
#' #use the demo_mix and demo_truemix3comp objects; the windows are found in the
#' #corresponding demo demo_intsurf and demo_intsurf3comp
#' demo_intsurf$window
#' normdens1=dnormmix(demo_mix, xlim= c(0, 1), ylim = c(0, 1))
#' plotmix_3d(normdens1, title1="Density of a normal mixture, 2 components")
#' #change the window
#' normdens1=dnormmix(demo_mix, xlim= c(-1, 1.5), ylim = c(-1, 1.5))
#' plotmix_3d(normdens1, title1="Density of a normal mixture, 2 components")
#' demo_intsurf3comp$window
#' normdens2=dnormmix(demo_truemix3comp, xlim= c(-1, 1), ylim = c(-2, 3))
#' plotmix_3d(normdens2, title1="Density of a normal mixture, 3 components")}
#'
#' @export
plotmix_3d<-function(dens_image,
          title1="3d Surface (Density or Intensity)",zlims=NULL,
          grayscale=FALSE)
{
  if (!is.im(dens_image)) {
    stop("dens_image must be of class im (a pixel image)")
  }
#  cat("\nTo save all open rgl graphs use Save_AllOpenRglGraphs.\n")
  xlims=dens_image$xrange
  ylims=dens_image$yrange
  LL=dens_image$dim[1]
  xcoord=seq(xlims[1],xlims[2],length=dens_image$dim[1])
  ycoord=seq(ylims[1],ylims[2],length=dens_image$dim[2])
  zcoord=t(dens_image$v)
  if(grayscale)
    cols <- gray.colors(100,start = 1, end = 0)[findInterval(zcoord, seq(min(zcoord), max(zcoord), length = 100))]
  else
  {
    jet.colors <- colorRampPalette(c("#00007F", "blue", "#007FFF", "cyan",
                                     "#7FFF7F", "yellow", "#FF7F00", "red",
                                     "#7F0000"))
    cols <- jet.colors(100)[findInterval(zcoord, seq(min(zcoord), max(zcoord), length = 100))]
  }
  if(is.null(zlims))
    zlims=c(0,1.1*max(zcoord))
  rgl::open3d(windowRect=c(50,50,1000,800),
              zoom=1.2)
  #rotation about the x-axis, 45 degrees
  U=rgl::par3d("userMatrix")
  rgl::par3d(userMatrix=
               rgl::rotate3d(U,pi/4,0,0,1))
  zmax=max(zcoord)
#  Rangez=zmax-min(zcoord);
  rgl::persp3d(x = xcoord, y = ycoord, z = zcoord,
               color = cols, xlab="x",ylab="y",zlab="",
               zlim=c(zlims[1]-0.01,zlims[2]),
               box = FALSE, axes = FALSE)
  rgl::axis3d('x')
  rgl::axis3d('y')
  rgl::axis3d('z-+',pos = c(xlims[1], ylims[2], 0))
  rgl::rgl.lines(c(xlims[1], xlims[1]),
                 c(ylims[2], ylims[2]),
                 c(0,max(zcoord)),
                 color = 'black')
  rgl::title3d(main=NULL)
  rgl::text3d(xlims[2],ylims[2],
              zlims[2]#+0.2*Rangez
              ,texts=title1)
  if(grayscale)
    rgl::bgplot3d(suppressWarnings(
    fields::image.plot(legend.only = TRUE,
                       smallplot= c(.8,.82,0.05,.7),
                       zlim = zlims,
                       col = gray.colors(100,start = 1, end = 0))))
  else
    rgl::bgplot3d(suppressWarnings(
      fields::image.plot(legend.only = TRUE,
                         smallplot= c(.8,.82,0.05,.7),
                         zlim = zlims,
                         col = jet.colors(100))))
}


#' Plot the true membership indicators
#'
#' @description
#' The function plots the true membership
#' indicators (or allocation variables)
#' of each point to one of the mixture components,
#' based on a generated \code{sppmix} object.
#' These are the true probabilities
#' of a point belonging to a component.
#'
#' For examples see
#'
#' \url{http://faculty.missouri.edu/~micheasa/sppmix/sppmix_all_examples.html
#' #plot_true_labels}
#'
#' @param pattern Object of class \code{sppmix}.
#' @param open_new_window Open a new window for the plot.
#'
#' @seealso \code{\link{rnormmix}},
#' \code{\link{to_int_surf}},
#' \code{\link[spatstat]{owin}},
#' \code{\link{rsppmix}}
#' @author Sakis Micheas
#' @examples
#' \donttest{
#' truemix <- rnormmix(m = 5, sig0 = .1, df = 5, xlim= c(-3, 3), ylim = c(-3, 3))
#' intsurf=to_int_surf(truemix, lambda = 100, win =spatstat::owin( c(-3, 3),c(-3, 3)))
#' pp1 <- rsppmix(intsurf,FALSE)
#' plot_true_labels(pp1)}
#'
#' @export
plot_true_labels <- function(pattern,open_new_window=FALSE)
{
  point=component=probability=NULL
  openwin_sppmix(check2open=open_new_window)

  if (length(pattern$comp)>0)
  {
    comp <- factor(pattern$comp)
  }
  else
  {
    stop("Missing true component labels.")
  }
  m=max(as.numeric(comp))
  n=pattern$n
  probs=matrix(0,n,m)
  for(i in 1:n)
    probs[i,comp[i]]=1
  labsy=1:m
  labsx=as.integer(seq(1,n,length=10))

  plot_df <- data.frame(probability = as.vector(probs),
                        point = 1:n,
                        component = rep(1:m, each =n))

  p<-ggplot(plot_df, aes(point, component, xend = point,
                      yend = component + 1,
                      col = probability,
                      xmin=1,xmax=n,
                      ymin=1,ymax=m)) +
    geom_segment(size = I(5)) +
    ggplot2::scale_color_gradient(low = "white", high = "grey18",
                                  guide = guide_colorbar(nbin = 100, barheight = 15)) +
    ggplot2::theme_classic() +
    ggplot2::theme(panel.border = ggplot2::element_rect(fill = NA, size = 1)) +
    ggplot2::labs(x = "Point", y = "Component", colour = "Probability") +
    ggplot2::scale_x_discrete(limits = labsx)+
    ggplot2::scale_y_discrete(limits = labsy+.5,labels=labsy)+
    add_title("True membership indicators", m = m, n = n)
  p
}



