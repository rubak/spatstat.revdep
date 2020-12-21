#' @title Function: Map a bipolar theme broken around a neutral value
#'
#' @description \code{mapBiPolar} generates a map of a bipolar theme variable
#'
#' @details The function \code{mapBiPolar} maps a \emph{bipolar variable} with a
#'   \emph{divergent color ramp} around a specific break values. A legend is generated.
#'   Below values are coded blue and above values red. Each branch is broken
#'   into 'quantiles'. Therefore, the number of class in each branch should be
#'   proportional to the number of observations in each class. NA's are permitted.
#'
#' @usage mapBiPolar(var.name,shape, break.value=0, neg.breaks=4, pos.breaks=neg.breaks,
#'                   map.title="", legend.title=deparse(substitute(var.name)),
#'                   legend.pos="bottomleft", legend.cex=1, add.to.map=FALSE)
#'
#' @param var.name A variable to be mapped in a bipolar theme. If it is in a data-frame,
#'   then the data-frame must be refered to, e.g., \code{df$var}
#' @param shape An existing spatial polygon or spatial polygon data-frame
#' @param break.value Neutral value separating the negative branch from the
#'   positive branch of the variable
#' @param neg.breaks Number of classes in the negative branch (default=4)
#' @param pos.breaks Number of classes in the positive branch (default=\code{neg.breaks})
#' @param map.title Character string with map title
#' @param legend.title Character string with legend title (default=\code{var.name})
#' @param legend.pos Location of legend in the map frame (default=\code{"bottomleft"})
#' @param legend.cex Relative font size of the legend
#' @param add.to.map Logical to start a new map frame if \code{FALSE} or overlay onto an
#'   existing map frame if \code{TRUE}
#' @export
#' @return \code{NULL}
#' @author Michael Tiefelsdorf <tiefelsdorf@@utdallas.edu>
#' @examples
#' library(maptools)
#' validTractShp <- tractShp[!is.na(tractShp$BUYPOW), ]         # Remove 2 tracts with NA's
#' mapColorQual(validTractShp$CITYPERI, validTractShp,
#'              map.title="Cities and Peripherie in Dallas County",
#'              legend.title="Regions")
#'
#' mapColorRamp(validTractShp$bad1500D, validTractShp, breaks=9,
#'              map.title="Density of Convenience Stores in Dallas County\nbw=1500 meters",
#'              legend.title="Junk Food")
#'
#' hist(tractShp$LRRmedD)
#' mapBiPolar(validTractShp$LRRmedD, validTractShp, break.value=0,
#'            neg.breaks=5, pos.breaks=5,
#'            map.title="LRR: log(f(junk food),f(healthy food))\nbw=medium",
#'            legend.title="log relative risk")
#'
mapBiPolar <- function(var.name,shape,
                       break.value=0,neg.breaks=4,pos.breaks=neg.breaks,
                       map.title="",legend.title=deparse(substitute(var.name)),
                       legend.pos="bottomleft",legend.cex=1,
                       add.to.map=FALSE) {
   ##
   ## Plot bipolar map theme for variable "var.name"
   ##
   #require(RColorBrewer); require(classInt); require(maptools)

   ## define quantile breaks and color assignment
   q.neg.breaks <- classInt::classIntervals((var.name[var.name < break.value]),
                                            n=neg.breaks, style="quantile")
   q.pos.breaks <- classInt::classIntervals((var.name[var.name > break.value]),
                                            n=pos.breaks, style="quantile")
   q.breaks <- c(q.neg.breaks$brks[-(neg.breaks+1)],break.value,q.pos.breaks$brks[-1])     # combine neg and pos over zero

   pal.neg <- RColorBrewer::brewer.pal(neg.breaks, "Blues")
   pal.pos <- RColorBrewer::brewer.pal(pos.breaks, "Reds")
   pal <- c(rev(pal.neg),pal.pos)                                                # combine palettes
   map.col <- pal[findInterval(var.name,q.breaks,rightmost.closed=T)]

   legendLabs <- maptools::leglabs(round(q.breaks,digits=3))

   if (anyNA(var.name)){
      map.col[is.na(pal)] <- grDevices::grey(0.96)                      # Set NA's to light grey
      pal <- c(pal[1:(pos.breaks+neg.breaks)],grDevices::grey(0.96))    # Augment legend color
      legendLabs <- c(legendLabs, "NA")        # Augment legend name
   }

   ## generate choropleth map
   sp::plot(shape, col=map.col, border=grDevices::grey(0.9), axes=TRUE, add=add.to.map)
   graphics::legend(legend.pos, title=legend.title, legend=legendLabs,
          cex=legend.cex, fill=pal, bty="n", ncol=1)
   graphics::title(map.title)
   graphics::box()
} # end:mapBiPolar
