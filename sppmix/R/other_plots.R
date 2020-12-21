get_states_bbox <- function(x)
{
  n=length(x$Names)
  xlims=rep(0,2)
  ylims=rep(0,2)
  xlims[1]=ylims[1]=100000000000
  xlims[2]=ylims[2]=-100000000000
  for(i in 1:n)
  {
    polies=x$Polies[[i]]
    for(j in 1:length(polies))
    {
      #get x-y range
      xlim=range(x$Polies[[i]][[j]][,1])
      ylim=range(x$Polies[[i]][[j]][,2])
      if(xlim[1]<xlims[1])xlims[1]=xlim[1]
      if(ylim[1]<ylims[1])ylims[1]=ylim[1]
      if(xlim[2]>xlims[2])xlims[2]=xlim[2]
      if(ylim[2]>ylims[2])ylims[2]=ylim[2]
    }
  }
  return(list(xlims=xlims,ylims=ylims))
}

get_states_info <- function(states,stateonly=TRUE)
{
  shp=sppmix::USAStatesCounties2016
  #retrieves the states requested
  #and all county info on those states
  if (is.null(states))
    stateids=rep(TRUE,length(shp$StateNames) )
  else
    #find the state indices corresponding to these states
    stateids=tolower(shp$StateNames) %in% tolower(states)
  info=NULL
  if(stateonly)
  {
    #build a list with all state selected
    info$Names=shp$StateNames[stateids]
    info$Polies=shp$StatePolygons[stateids]
  }
  else
    #build a list with all counties for
    #the selected states
  {
    dat=shp$CountiesbyState[stateids]
    info$Names=NULL
    info$Polies=NULL
    for(i in 1:length(dat))
    {
      info$Names=c(info$Names,dat[[i]]$CountyName)
      info$Polies=c(info$Polies,dat[[i]]$CountyPolies)
    }
  }
  return(info)
}

#' Visualization of USA states and their counties
#'
#' @description
#' The function plots the requested USA state
#' or county boundaries and additional information if requested or
#' if certain parameters are supplied. We use this function
#' for visualization of geostatistical data, in particular,
#' (Marked) IPPPs.
#'
#' For examples see
#'
#' \url{http://faculty.missouri.edu/~micheasa/sppmix/sppmix_all_examples.html
#' #PlotUSAStates}
#'
#' @param showcounties Logical to denote that we want a plot of counties.
#' Default is FALSE. Setting this to TRUE will show all the counties
#' for the states passed in the \code{states} parameter.
#' @param states A vector of state names. Set to \code{NULL} to
#' request all states or \code{ContinentalUSA_state_names} to
#' show only the continental USA states.
#' @param showcentroids Logical requesting to
#' show centroids for each state or county. These centroids are returned in
#' a \code{\link[spatstat.geom]{ppp}} object. The centroid is chosen
#' so that it is always within the state or county boundaries.
#' @param typecentroid If \code{showcentroids=TRUE}
#' we can display either the average of the boundary (\code{typecentroid=0})
#' or the "marker point" of the state or county (\code{typecentroid=1}). For convex
#' states or counties, the latter point
#' is the most south-western point of the state or county.
#' @param shownames Logical to display the
#' names of the states for \code{showcounties=FALSE}
#' or counties for \code{showcounties=TRUE}.
#' @param showmarks Logical to display the
#' mark values given to each state for \code{showcounties=FALSE}
#' or county for \code{showcounties=TRUE}.
#' @param grayscale Logical to request plots in grayscale.
#' @param open_new_window Logical to request plotting in a new graphics window.
#' @param main A character string to serve as the main title for the plot.
#' @param guidemain A character string to be used as the title for
#' the guide used (legend or colorbar).
#' @param discretelevels Logical indicating that the marks are discrete valued.
#' @param levels When \code{discretelevels=TRUE}, the parameter \code{levels}
#' contains all the possible discrete levels (marks). This is a
#' vector of integers or strings. Default is \code{1:3}.
#' @param showplot Logical requesting to show the plot. Set to FALSE if
#' you want to simply retrieve the centroids of the states or counties,
#' in which case the plot will not be created.
#' @param plotlevels Logical requesting that
#' the levels (marks) of each state or county
#' are displayed. If \code{marks} is not
#' supplied, then for \code{discretelevels=TRUE}
#' the mark of each state or county is
#' uniformly generated over the values of
#' \code{levels}, otherwise the marks are
#' uniform in (0,1) (probabilities). If
#' \code{marks} are given, then they are
#' used to appropriately paint a state or county.
#' @param marks A vector of length equal to
#' the number of states or counties requested,
#' containing the mark values for each state or county.
#' A mark is an integer pointing to an element from the vector
#' \code{levels} for \code{discretelevels=TRUE},
#' otherwise a real number.
#' @param pp Optionally, a point pattern as an object of type \code{\link[spatstat.geom]{ppp}}
#' to be displayed over the created plot. The window of this point pattern will be used as
#' the window of observation (overrides the window in the
#' \code{surf} parameter).
#' @param surf Optionally, an intensity surface
#' as an object of type \code{intensity_surface}
#' or an image object of class \code{\link[spatstat.geom]{im}}
#' to be plotted first and then the map will be
#' displayed over this field. Supplying this
#' parameter sets the flag \code{plotlevels=FALSE}
#' automatically. The window of this intensity surface will be used as
#' the window of observation.
#' @param boundarycolor A specific color to use for drawing boundaries.
#' Default is "black". Set to \code{NULL} if you do not want boundaries drawn.
#' @param namescolor A specific color to use
#' for drawing the state or county names when
#' \code{plotnames=TRUE}. Default is "black".
#' @param ppsize Size used in plotting the points. Default is 1.
#'
#' @details Note that we use the state and county longitude and latitude boundaries in
#' the \code{\link{USAStatesCounties2016}} object.
#'
#' @return A list containing the following components:
#' \item{PPPcent}{The centroids of the states or counties requested, returned as a marked point pattern.}
#' \item{PPPMarker}{The marker points of the states or counties requested, returned as a marked point pattern.}
#' \item{itemnames}{Vector of strings containing all items processed (i.e., either all state names or all county names).}
#' \item{p}{The created plot, otherwise NULL.}
#'
#' @author Sakis Micheas and Jiaxun Chen
#'
#' @seealso \code{\link{est_MIPPP_cond_loc}},
#' \code{\link{est_mix_damcmc}},
#' \code{\link{est_mix_bdmcmc}},
#' \code{\link{plot_CompDist}},
#' \code{\link{drop_realization}},
#' \code{\link{GetBDTable}},
#' \code{\link{GetBDCompfit}},
#' \code{\link{plotmix_2d}},
#' \code{\link{GetBMA}},
#' \code{\link{plot_MPP_probs}},
#' \code{\link{GetMAPEst}}
#' @examples
#' \donttest{
#' #plot the continental USA with uniformly sampled discrete marks from 10 different levels
#' ret=PlotUSAStates(states=ContinentalUSA_state_names, levels=1:10, grayscale = FALSE,
#'  shownames=TRUE, plotlevels =TRUE, discretelevels=TRUE, main="Continental USA (generated levels)")
#' #now use continuous marks
#' ret=PlotUSAStates(states=ContinentalUSA_state_names, shownames=FALSE, discretelevels=FALSE,
#'  main="Continental USA (generated probabilities)", guidemain="Probability", showcentroids = FALSE)
#' #Fit an IPPP to the California Earthquake data
#' fitDA=est_mix_damcmc(CAQuakes2014.RichterOver3.0, 8, L = 20000)
#' #get the surface of Maximum a Posteriori estimates
#' MAPsurf=GetMAPEst(fitDA)
#' #plot the states and the earthquake points along with the fitted MAP IPPP intensity surface
#' ret=PlotUSAStates(states=c('California','Nevada','Arizona'), showcentroids=FALSE,
#'  shownames=TRUE, main="Earthquakes in CA, 2014", pp=CAQuakes2014.RichterOver3.0, surf=MAPsurf,
#'  boundarycolor="white", namescolor="white")
#' #Visualize the Tornado data about MO
#' #plot the states and the tornado points
#' ret=PlotUSAStates(states=c('Iowa','Arkansas','Missouri','Illinois','Indiana','Kentucky',
#'  'Tennessee','Kansas','Nebraska','Texas','Oklahoma','Mississippi','Alabama','Louisiana'),
#'  showcentroids=FALSE, shownames=TRUE, plotlevels = FALSE, main="Tornadoes about MO, 2011",
#'  pp=Tornadoes2011MO)
#' #Visualize aggregate income levels in MO by county using data from the American Community
#' #Survey (ACS)
#' #plot in the original scale first; here we pass the marks vector which contains the aggregate
#' #income values of Missourian counties
#' ret=PlotUSAStates(showcounties=TRUE, states=c('Missouri'), showcentroids=TRUE, typecentroid=1,
#'  discretelevels=FALSE, shownames=TRUE, plotlevels=TRUE, marks=MOAggIncomeLevelsPerCounty,
#'  main="Aggregate Income in MO, 2014", guidemain = "Income level", namescolor="gray",
#'  boundarycolor="gray")
#' #plot in the log scale
#' ret=PlotUSAStates(showcounties=TRUE, states=c('Missouri'), showcentroids=TRUE, typecentroid=1,
#'  discretelevels=FALSE, shownames=TRUE, plotlevels=TRUE, marks=log(MOAggIncomeLevelsPerCounty),
#'  main="Aggregate Income in MO, 2014", guidemain = "Income level\n(log scale)", namescolor="gray",
#'  boundarycolor="gray")
#' #plot the marker points, county boundaries and names
#' ret=PlotUSAStates(showcounties=TRUE, states=c('Missouri'), showcentroids=TRUE, typecentroid = 1,
#'  discretelevels=FALSE, shownames=TRUE, plotlevels=FALSE, marks=log(MOAggIncomeLevelsPerCounty),
#'  main="Marker points for Missouri counties")
#' #now plot only the marker points, we treat this as a marked IPPP
#' ret=PlotUSAStates(showcounties=TRUE, states=c('Missouri'), showcentroids=TRUE, typecentroid = 1,
#'  discretelevels=FALSE, shownames=FALSE, plotlevels=FALSE, marks=log(MOAggIncomeLevelsPerCounty),
#'  main="Marker points for Missouri counties", boundarycolor = NULL)
#' #let us discretize log(income) to 3 levels; low if <=20, average if >20 and <=23, and high if >23
#' newmarks=rep(0,length(MOAggIncomeLevelsPerCounty))
#' newmarks[log(MOAggIncomeLevelsPerCounty)<=20]=1
#' newmarks[log(MOAggIncomeLevelsPerCounty)>20 & log(MOAggIncomeLevelsPerCounty)<=23]=2
#' newmarks[log(MOAggIncomeLevelsPerCounty)>23]=3
#' table(newmarks)
#' levels=c("low","average","high")
#' ret=PlotUSAStates(showcounties=TRUE, states=c('Missouri'), showcentroids=TRUE, typecentroid=1,
#'  discretelevels=TRUE, shownames=TRUE, plotlevels=TRUE, main="Aggregate Income in MO, 2014",
#'  marks=newmarks, levels=levels, guidemain = "Income level", namescolor="gray",
#'  boundarycolor="gray")
#' #now fit a marked IPPP model, use the PP of marker points
#' MPP=ret$PPPMarker
#' mpp_est <- est_MIPPP_cond_loc(MPP,r=1, hyper=0.2)
#' plot_MPP_probs(mpp_est)
#' #now obtain a BDMCMC fit for the ground process this way we can cluster the data
#' BDMCMCfit <- est_mix_bdmcmc(MPP,m=10,L = 50000)
#' plot_CompDist(BDMCMCfit)
#' #use the original output of BDMCMC and apply 10% burnin (default)
#' BDMCMCfit=drop_realization(BDMCMCfit)
#' #get the realizations corresponding to the MAP number of components
#' BDtab=GetBDTable(BDMCMCfit,FALSE)#retrieve frequency table and MAP estimate for
#' #the number of components
#' MAPm=BDtab$MAPcomp
#' BDMCMCfitMAPcomp=GetBDCompfit(BDMCMCfit,MAPm)
#' BDMCMCfitMAPcompgens=BDMCMCfitMAPcomp$BDgens
#' MAPsurf=GetMAPEst(BDMCMCfitMAPcompgens)
#' plotmix_2d(MAPsurf,MPP)+add_title(
#'  "IPPP intensity surface of MAP estimates (MAP number of components)",
#'  lambda =MAPsurf$lambda, m=MAPsurf$m, n=MPP$n, L=MAPsurf$L)
#' plot_ind(BDMCMCfitMAPcompgens)
#' ret=PlotUSAStates(showcounties=TRUE, states=c('Missouri'),
#'  showcentroids=TRUE, typecentroid=1, discretelevels=TRUE, shownames=TRUE,
#'  main="Ground surface of MAP estimates", marks=newmarks, levels=levels,
#'  guidemain = "Income level", namescolor="gray", boundarycolor="gray",
#'  pp=MPP, surf=MAPsurf)
#' #obtain and plot the Bayesian model average; first drop the bad realizations
#' BDMCMCfit=drop_realization(BDMCMCfit,(BDMCMCfit$Badgen==1))
#' BMAest=GetBMA(BDMCMCfit)
#' ret=PlotUSAStates(showcounties=TRUE, states=c('Missouri'),
#'  showcentroids=TRUE, typecentroid=1, discretelevels=TRUE, shownames=TRUE,
#'  main="Bayesian model average ground intensity surface", marks=newmarks,
#'  levels=levels, guidemain = "Income level", namescolor="gray",
#'  boundarycolor="gray", pp=MPP, surf=BMAest)}
#'
#' @export
PlotUSAStates=function(showcounties=FALSE,
                       states="Missouri",showcentroids=TRUE,
                       typecentroid=0,shownames=FALSE,
                       showmarks=FALSE,grayscale=FALSE,
                       open_new_window = FALSE,
                       main="States (true levels)",
                       guidemain="Level",discretelevels=TRUE,
                       levels=1:3,showplot=TRUE,plotlevels=TRUE,
                       marks,pp,surf,boundarycolor="black",
                       namescolor="black",ppsize=1.0)
{
  x=y=id=value=statenames=colval=NULL
  all_info=get_states_info(
    states=states,stateonly=!showcounties)
  all_names <- factor(as.character(all_info$Names))
  n_info=length(all_names)
  if(showcounties)
    cat("\nNumber of counties requested is",n_info,"\n")
  else
    cat("\nNumber of states requested is",n_info,"\n")
  if(missing(surf) & missing(pp))
  {
    limits=get_states_bbox(all_info)
    xlims=limits$xlims
    ylims=limits$ylims
  }
  else
  {
    if(!missing(surf))
    {
      xlims=surf$window$xrange
      ylims=surf$window$yrange
    }
    if(!missing(pp))
    {
      xlims=pp$window$xrange
      ylims=pp$window$yrange
    }
  }
  win=spatstat.geom::owin(xlims,ylims)
  centers=matrix(0,n_info,2)
  #lower-left points
  ll_pt=matrix(0,n_info,2)
  nlevels=length(levels)
  levelnumbers=1:nlevels
  allstatecoords=NULL
  #length of ids is the same as the number
  #of polygons to draw
  #  ids=factor(1:nstate_info)
  #value is the mark to use for each polygon
  genMARKS=FALSE
  if(missing(marks))
  {
    marks=rep(0,n_info)
    marks_names=marks
    genMARKS=TRUE
  }
  else
  {
    if(length(marks)!=n_info)
    {
      warning("The number of marks passed is not the same as the number of states passed. The mark vector should have ",n_info," elements. Will generate the marks.")
      marks=rep(0,n_info)
      marks_names=marks
      genMARKS=TRUE
    }
  }
  rep_ids=NULL
  ids=NULL
  plotmarkvals=NULL
  newid = 0
  for (i in 1:n_info)
  {
    if(genMARKS)
    {
      if(discretelevels)
      {
        marks[i]=sample(x=levelnumbers,1,replace=TRUE)
        marks_names[i]=levels[marks[i]]
      }
      else
        marks[i]=runif(1)
    }
    statepolies=all_info$Polies[[i]]
    numpolies=length(statepolies)
    #    if(numpolies>1){print(state_names[i]);cat("\n",numpolies,"\n")}
    statecoords=NULL
    curid=NULL
    for(j in 1:numpolies)
    {
      statepoly=statepolies[[j]]
      statecoords<-rbind(statecoords,statepoly)

      nscoord=nrow(statepoly)
      newid=newid + 1
      ids=c(ids,newid)
      curpolyid=rep(newid,nscoord)
      curid=c(curid,curpolyid)
      plotmarkvals=c(plotmarkvals,marks[i])
    }
    #    nscoord=nrow(statecoords)
    #    curid=rep(ids[i],nscoord)
    rep_ids=c(rep_ids,curid)
    allstatecoords<-rbind(allstatecoords,statecoords)
    centers[i,]<-apply(statecoords,2,mean)

    gotit=FALSE
    for(j in 1:numpolies)
    {
      statecoords1=statepolies[[j]]
      if(CheckInPoly(statecoords1,
                     centers[i,]))
      {
        gotit=TRUE
        break
      }
    }
    #if the center is not in the polygon
    #sample it randomly from the poly
    while(!gotit)
    {
      centers[i,1]=runif(1,min(statecoords[,1])
                         ,max(statecoords[,1]))
      centers[i,2]=runif(1,min(statecoords[,2])
                         ,max(statecoords[,2]))
      for(j in 1:numpolies)
      {
        statecoords1=statepolies[[j]]
        if(CheckInPoly(statecoords1,
                       centers[i,]))
        {
          gotit=TRUE
          break
        }
      }
    }
    minlatind=which.min(statecoords[,2])
    allpts=statecoords[(statecoords[,2]==statecoords[minlatind,2]),]
    if(is.matrix(allpts))
    {
      minlonind=which.min(allpts[,1])
      ll_pt[i,]=allpts[minlonind,]
    }
    else
      ll_pt[i,]=allpts
  }
  values<-data.frame(id=ids,value=plotmarkvals)
  positions <- data.frame(id = rep_ids,
                          x = allstatecoords[,1],
                          y = allstatecoords[,2])
  if(showplot)
  {
    openwin_sppmix(check2open=open_new_window)
    if(missing(surf))
      p <- ggplot()
    else
    {
      plotlevels=FALSE
      #      p <- plotmix_2d(surf,open_new_window = FALSE,truncate=FALSE,grayscale = grayscale)
      if(grayscale)
        cols <- gray.colors(100,start = 1, end = 0)
      else
        cols <- c("#00007F", "blue", "#007FFF", "cyan", "#7FFF7F", "yellow",
                  "#FF7F00", "red", "#7F0000")

      if(is.im(surf))
      {
        density_df=as.data.frame(surf)
        p <- ggplot(data=density_df,aes(x=x, y=y))+
          geom_raster(aes(fill = value), interpolate = TRUE) +
          scale_fill_gradientn(colors = cols) +
          guides(fill = guide_colorbar(
            title = "Elevation",nbin = 100, barheight = 15))
      }
      else
      {
        est_intensity <- dnormmix(surf,
                                  xlim =xlims,ylim = ylims,
                                  L = 128, truncate = FALSE)

        density_df=as.data.frame(est_intensity)
        p <- ggplot2::ggplot(data=density_df,
                             aes(x=x, y=y))
        p<-p + ggplot2::geom_raster(aes(fill = value), interpolate = TRUE) +
          ggplot2::scale_fill_gradientn(colors = cols) +
          ggplot2::guides(
            fill = guide_colorbar(title = "Elevation",nbin = 100,
                                  barheight = 15))
      }
    }
    datapoly <- merge(values, positions, by=c("id"), sort=FALSE)
    if(plotlevels)
    {
      if(is.null(boundarycolor))
        p <-p+ geom_polygon(data=datapoly,
                            aes(x=x, y=y,fill=value, group=id))
      else
        p <-p+ geom_polygon(data=datapoly,
                            aes(x=x, y=y,fill=value, group=id),
                            color=boundarycolor)
    }
    else
    {
      if(!is.null(boundarycolor))
        p<-p+geom_path(data=datapoly,
                       aes(x=x,y=y, group=id),
                       color=boundarycolor)
    }
    #    cat("pass1");# return()
    p<-p+ggtitle(main)+
      labs(x = "Longitude",y = "Latitude")+
      theme_classic() +
      theme(panel.border = element_rect(fill = NA, size = 1))+
      coord_cartesian(xlim =xlims,
                      ylim =ylims)#,expand=FALSE)
    #    cat("pass2");# return()
    if(plotlevels)
    {
      if(discretelevels)
      {
        if(grayscale)
        {
          lowcol="grey"
          highcol="black"
        }
        else
        {
          #          lowcol="#00007F"
          #          highcol="#7F0000"
          lowcol= "#FFFF00"
          highcol="#FF7F00"
        }
        p<-p+
          scale_fill_continuous(low=lowcol,
                                high=highcol,
                                breaks=sort(unique(marks)),
                                labels=sort(unique(marks)),
                                guide = guide_legend(
                                  reverse=TRUE,title = guidemain)
          )
        #       cat("pass3"); return()
      }
      else
      {
        if(grayscale)
          cols <- gray.colors(100,start = 1, end = 0)
        else
          cols <- c("#00007F", "blue", "#007FFF", "cyan", "#7FFF7F", "yellow",
                    "#FF7F00", "red", "#7F0000")
        p<-p+scale_fill_gradientn(colors = cols
                                  #,limits = c(0,1),breaks=c(0,.25,.5,.75,1)
        ) +
          guides(fill = guide_colorbar(
            title =guidemain,nbin = 100,
            barheight = 15))
      }
    }
    #print(p); cat("pass4");# return()
    cents=data.frame(x=centers[,1],
                     y=centers[,2],
                     statenames=as.character(all_names))
    llpts<-data.frame(x=ll_pt[,1],
                      y=ll_pt[,2],
                      colval=marks)
    if(shownames)
      p<-p+geom_text(data=cents,
                     color=namescolor,#check_overlap=TRUE,
                     size=2.5,
                     aes(x=x,y=y,label=statenames,
                         show.legend=FALSE))
    #print(p);cat("pass5");# return()
    if(showmarks)
      p<-p+geom_text(data=llpts,
                     color=namescolor,#check_overlap=TRUE,
                     aes(x=x,y=y,label=colval,# sort(unique(marks)),
                         show.legend=FALSE))
    #   print(p);cat("pass6");# return()
    if(showcentroids)
    {
      if(typecentroid==0)
        dfpoint=cents
      else
        dfpoint=llpts
      if(grayscale)
        p<-p+geom_point(#inherit.aes=FALSE,
          data=dfpoint,aes(x=x,y=y),
          #data=llpts,aes(x=llpts[,1],y=llpts[,2]),
          color="red",show.legend=FALSE)
      else
        p<-p+geom_point(#inherit.aes=FALSE,
          data=dfpoint,aes(x=x,y=y),
          #data=llpts,aes(x=llpts[,1],y=llpts[,2]),
          color="grey",show.legend=FALSE)
    }
    if(!missing(pp))
    {
      p<-p+geom_point(#inherit.aes=FALSE,
        data=as.data.frame(pp),
        aes(x=pp$x,y=pp$y),size=ppsize,
        color="black",show.legend=FALSE)
    }
    if(!missing(pp) | !missing(surf))
      p<-p+coord_cartesian(xlim =xlims,
                           ylim =ylims,expand=FALSE)
    p<-p+ggplot2::theme(plot.margin=ggplot2::unit(c(1,1,1,0), "lines"))
    print(p)
  }
  else
    p=NULL
  PPPcent<-ppp(centers[,1],centers[,2],
               window=win,marks=marks, check=FALSE)
  PPPMarker<-ppp(ll_pt[,1],ll_pt[,2],
                 window=win,marks=marks, check=FALSE)
  return(list(PPPcent=PPPcent,
              PPPMarker=PPPMarker,
              itemnames=all_names,
              p=p))
}
