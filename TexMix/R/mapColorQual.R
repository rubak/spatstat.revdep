#' @title Function: Maps a qualitative theme for a maximum of 12 categories
#'
#' @description \code{mapColorQual} generates a map of a qualitative variable
#'
#' @details The function \code{mapColorQual} maps a \emph{categorical variable} with a
#'   \emph{set of distinct colors}. A legend is generated. The maximum number of
#'   valid factor levels should not exceed 12. NA's are permitted.
#'
#' @usage mapColorQual(var.name, shape, map.title="", legend.title=deparse(substitute(var.name)),
#'                     legend.pos="bottomleft", legend.cex=1, add.to.map=FALSE)
#'
#' @param var.name A factor, perhaps with NA's, to be mapped with a maximum of 12 categories.
#'   If the factor is in a data-frame, then the data-frame must be explicitly
#'   referred to, e.g., \code{df$var}
#' @param shape An existing spatial polygon or spatial polygon data-frame
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

mapColorQual <- function(var.name, shape,
                         map.title="", legend.title=deparse(substitute(var.name)),
                         legend.pos="bottomleft", legend.cex=1, add.to.map=FALSE) {
  ##
  ## Plot a qualitative colors for factor "var.name"
  ##
  #require(maptools); require(RColorBrewer); require(classInt)
  if (!is.factor(var.name)) stop("plotColorQual: Not a factor.")
  if (length(levels(var.name)) > 12) stop("The maximum number of factor levels exceeds 12")

  qualVal <- as.numeric(unclass(var.name))
  qualName <- levels(var.name)
  pal.Qual <- RColorBrewer::brewer.pal(12,"Set3")
  map.col <- pal.Qual[qualVal]
  if (anyNA(qualVal)){
    map.col[is.na(map.col)] <- grDevices::grey(0.96)                   # Set NA's to light grey
    pal.Qual <- c(pal.Qual[1:length(qualName)],grDevices::grey(0.96))  # Augment legend color
    qualName <- c(qualName,"NA")                                      # Augment legend name
  }

  ## generate choropleth map
  sp::plot(shape, col=map.col, border=grDevices::grey(0.9), axes=TRUE, add=add.to.map)
  graphics::legend(legend.pos, title=legend.title, legend=qualName,
                   cex=legend.cex, fill=pal.Qual[1:length(qualName)],bty="n",
                   ncol=1)
  graphics::title(map.title)
  graphics::box()
} # end::mapColorQual
