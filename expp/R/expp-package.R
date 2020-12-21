

#' Study area boundary.
#' 
#' \code{SpatialPolygonsDataFrame} Study area boundary of two Blue Tit
#' populations: Kolbeterberg, Vienna, Austria (1998 through 2004) and
#' Westerholz, Bavaria, Germany (2007 through 2011) .
#' 
#' 
#' @name bluetit_boundary
#' @docType data
#' @format A SpatialPolygonsDataFrame with 12 SpatialPolygons.  
#' \describe{
#'  \item{list("year_")}{numeric. The year of the observation.} 
#'  }
#' @keywords datasets
#' @examples
#' 
#' data(bluetit_boundary)
#' summary(bluetit_boundary)
#' 
NULL





#' Blue Tit breeding data.
#' 
#' Breeding data recorded for two Blue Tit populations in Kolbeterberg, Vienna,
#' Austria (1998 through 2004) and Westerholz, Bavaria, Germany (2007 through
#' 2011) . The data set contains breeding attempts locations, the respective
#' social pair, and several individual and nest parameters.
#' 
#' 
#' @name bluetit_breeding
#' @docType data
#' @format A data frame with 1025 observations on the following 10 variables.
#' \itemize{ 
#'      \item year_       numeric. The year of the observation.
#'      \item id          numeric. The identity of the nest box in which the
#'                        breeding attempt took place. 
#'      \item x           numeric. The east-west location of the nest box. 
#'      \item y           numeric. The north-south location of the nest box. 
#'      \item female      character. The identity of the female. 
#'      \item male        character. The identity of the male. 
#'      \item layingDate  numeric. The day of the year when the first egg was produced. 
#'      \item male_age    character. The age class of the male ('juv' = 1st year breeder; 'adult' = older)
#'      \item male_tarsus numeric. tarsus length (mm))
#'      \item study_area  character. The study area name.  
#'      }
#' @keywords datasets
#' @examples
#' 
#' data(bluetit_breeding)
#' head(bluetit_breeding)
#' 
NULL





#' Blue tit extra-pair paternity data.
#' 
#' \code{data.frame} Extra-pair paternity data recorded for two Blue Tit
#' populations in Kolbeterberg, Vienna, Austria (1998 through 2004) and
#' Westerholz, Bavaria, Germany (2007 through 2011) .
#' 
#' 
#' @name bluetit_epp
#' @docType data
#' @format A data frame with 425 observations on the following 3 variables.
#' \itemize{ 
#'      \item year_   numeric. The year of the observation.
#'      \item female   character. The female involved in the respective EPP event.
#'      \item male     character. The male involved in the respective EPP event.
#'      }
#' @keywords datasets
#' @examples
#' 
#' data(bluetit_epp)
#' head(bluetit_epp)
#' 
NULL


#' Tools and data to accompany Schlicht, Valcu and Kempenaers "Spatial patterns
#' of extra-pair paternity: beyond paternity gains and losses"
#' 
#' The expp package provides classes and functions for the investigation of the
#' probability of having extra-pair young within local networks of breeding
#' pairs including both realized and potential extra-pairings.
#' 
#' \tabular{ll}{ Package: \tab expp\cr Type: \tab Package\cr Version: \tab
#' 1.0\cr Date: \tab 2014-Aug-04\cr License: \tab GPL-3\cr } See
#' \code{help(epp)} and \code{vignette('expp') }
#' 
#' @name expp-package
#' @aliases expp-package expp
#' @docType package
#' @section Functions: \tabular{ll}{ \code{\link{epp}} \tab Final
#' data-transformation to male-female combinations and their extra-pair levels
#' \cr \code{\link{eppSimDat}} \tab "Toy"-dataset creation to investigate
#' potential Type I error rate inflation for models where the datapoints are
#' male-female combinations \cr
#' 
#' \cr ---------------------------------------- \cr
#' \code{\link{DirichletPolygons}} \tab Territory calculation via Dirichlet
#' tesselation\cr \code{\link{eppMatrix}} \tab \code{data.frame} to
#' \code{eppMatrix} object\cr \code{\link{neighborsDataFrame}} \tab \code{nb}
#' object to \code{data.frame}\cr \code{\link{SpatialPointsBreeding}} \tab
#' \code{data.frame} to \code{SpatialPointsBreeding object} \cr
#' 
#' }
#' @author Mihai Valcu and Lotte Schlicht \cr Maintainer: Mihai Valcu
#' <valcu@@orn.mpg.de>
#' @references 
#' 
#' Schlicht, Lotte, Mihai Valcu, and Bart Kempenaers. 
#' "Spatial patterns of extra-pair paternity: beyond paternity gains and losses." 
#'  Journal of Animal Ecology 84.2 (2015): 518-531.
#' @keywords package
NULL



