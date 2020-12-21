#' @title Layout of box-plots: variables by factor
#'
#' @description \code{plotBoxesByFactor} generates box-plots of several variable by a factor.
#'
#' @details
#' This function organizes several box-plots of metric variables broken
#' by a factor variable in to a layout. By default 3 box-plots are shown
#' per row and the variables are standardized by the z-transformation. As the number
#' of factor levels increases the number of plots per row should be decreased to maintain
#' a reasonable visual resolution.
#'
#' @usage plotBoxesByFactor(xVars, groups, ncol=3, zTrans=TRUE, varwidth=FALSE)
#'
#' @param xVars data-frame of metric variables.
#' @param groups factor variable.
#' @param ncol number of layout columnes.
#' @param zTrans logical. \code{xTrans=TRUE} performs a z-transformation of all
#' variables.
#' @param varwidth logical. \code{varwidth=TRUE} makes the box width proportional
#' to the number of observation used to generate the box for a specific factor
#' level.
#' @export
#' @return NULL
#' @author Michael Tiefelsdorf <tiefelsdorf@@utdallas.edu>
#' @examples
#' varsKeep <- c("PCTWHITE","PCTBLACK","PCTASIAN","PCTHISPAN","MEDAGE","MEDVALHOME")
#' myData <- tractShp@data
#' plotBoxesByFactor(myData[,varsKeep], tractShp$CITYPERI, ncol=2, zTrans=TRUE, varwidth=FALSE)

plotBoxesByFactor <- function(xVars, groups, ncol=3, zTrans=TRUE, varwidth=FALSE){
  xVars <- xVars[ , sapply(xVars, is.numeric)] # only metric variables
  namesVars <- names(xVars)
  if (zTrans) xVars <- as.data.frame(scale(xVars))
  nVars <- length(namesVars)
  nrow <- nVars%/%ncol+1                      # calculate number of rows
  graphics::layout(matrix(c(1:(nrow*ncol)), nrow=nrow, ncol=ncol, byrow=T))
  for (i in 1:nVars){                         # cycle over all variables
    graphics::boxplot(xVars[,i]~groups, varwidth=varwidth, xlab="Group", ylab=namesVars[i],
                      main=paste("Feature:", namesVars[i]))
    graphics::abline(h=mean(xVars[,i], na.rm=T), lty=5, col="red")
  }
  graphics::layout(c(1))                                # reset layout
} ## end::plotBoxesByFactor

