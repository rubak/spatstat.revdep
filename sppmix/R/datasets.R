#' @title Datasets
#' @aliases Datasets
#' @description
#' Several datasets have been included in the
#' \code{sppmix} package. They are all open source
#' datasets that have been processed into
#' \code{\link[spatstat.geom]{ppp}} objects. Many examples and tutorials
#' use these datasets.
#'
#' For basic examples see
#'
#' \url{http://faculty.missouri.edu/~micheasa/sppmix/sppmix_all_examples.html
#' #Datasets}
#'
#' @format All datasets are objects of type \code{\link[spatstat.geom]{ppp}} and \code{sppmix}, except for object \code{USAStatesCounties2016}.
#' @author Sakis Micheas
#' @examples
#' \donttest{
#' ChicagoCrime2015
#' summary(ChicagoCrime2015)
#' plot(ChicagoCrime2015)+add_title("Chicago Crime, 2015")
#' CAQuakes2014.RichterOver3.0
#' summary(CAQuakes2014.RichterOver3.0)
#' plot(CAQuakes2014.RichterOver3.0)+add_title("Earthquakes in California, 2014")
#' Tornadoes2011MO
#' summary(Tornadoes2011MO)
#' plot(Tornadoes2011MO)+add_title("Tornado events about Missouri, 2011")
#' }
#'
#' @export
Datasets <- function()
{
  cat("\nDatasets included in the sppmix package:")
  cat("\nChicagoCrime2015, CAQuakes2014 and Tornadoes2011MO.")
}

#' @rdname Datasets
#' @description
#'
#' \code{ChicagoCrime2015}
#'
#' The 2015 Chicago crime dataset (\url{http://data.cityofchicago.org/Public-Safety/Crimes-2001-to-present/ijzp-q8t2})
#' contains the reported locations of homicide cases in Chicago during the year 2015. The \code{sppmix}
#' package object \code{ChicagoCrime2015} is a marked point pattern containing the crime locations
#' for two types of crimes; kidnappings coded as 0, and homicides coded as 1.
#' In particular, there are 654 reported incident locations, 188 kidnappings, and 466 homicides.
"ChicagoCrime2015"

#' @rdname Datasets
#' @description
#'
#' \code{ChicagoArea}
#'
#' The \code{ChicagoArea} object contains the coordinates for the different neighborhoods of the city of Chicago.
"ChicagoArea"

#' @rdname Datasets
#' @description
#'
#' \code{CAQuakes2014}
#'
#' The Southern California Earthquake Data Center
#' (SCEDC, Southern California Earthquake Center, Caltech,
#' Dataset: doi:10.7909/C3WD3xH1), operates at the Seismological
#' Laboratory at Caltech and is the primary archive of
#' seismological data for southern California, recording the
#' seismological activity from 1932 to present. Here we concentrate
#' on the year 2014, and all this data is
#' contained in a marked point process object
#' named \code{CAQuakes2014}. There are two
#' continuous variables containing marks for the events,
#' including variable \code{mag} which represents
#' the magnitude of the earthquake and variable \code{depth}
#' representing the depth.
"CAQuakes2014"

#' @rdname Datasets
#' @description
#'
#' \code{CAQuakes2014.RichterOver3.0}
#'
#' In many examples we further restrict to events with magnitudes
#' over 3.0 in the Richter scale. The latter marked point pattern
#' is aplty named \code{CAQuakes2014.RichterOver3.0}.
"CAQuakes2014.RichterOver3.0"

#' @rdname Datasets
#' @description
#'
#' \code{TornadoesAll}
#'
#' The National Oceanic and Atmospheric Administration (NOAA, \url{http://www.noaa.gov/}) is a
#' U.S.A. agency tasked with the dissemination of daily weather forecasts
#' and severe storm warnings, as well as, conducting climate monitoring
#' among many other important tasks. The Storm Prediction Center of NOAA contains
#' important information on tornado occurrences throughout the U.S., starting from 1950
#' all the way to the present. All this information is contained in
#' the data.frame object \code{TornadoesAll}. The variables (columns) included are as follows:
#' 1="RecordNumber", 2="Year", 3="Month",
#' 4="Day", 5="Date:yyyy-mm-dd", 6="Time:HH:MM:SS",
#' 7="State", 8="Fscale", 9="Injuries"
#' 10="Fatalities", 11="Estimated property loss",
#' 12="Estimated crop loss", 13="Starting latitude",
#' 14="Starting longitude", 15="Length in miles",
#' and 16="Width in yards".
#'
#' Note that each event (row) is marked using
#' one of 6 levels (variable Fscale), each denoting the strength of a tornado
#' in the Fujita scale (Enhanced Fujita scale
#' after January, 2007), with 0 denoting minimal damage and
#' 5 indicating complete destruction.
"TornadoesAll"

#' @rdname Datasets
#' @description
#'
#' \code{Tornadoes2011MO}
#'
#' In many examples we use the marked point pattern for year
#' 2011 (a single time snap-shot), contained in
#' the object \code{Tornadoes2011MO}, in order to study the behavior
#' of the events of that year; this year is of particular interest
#' since on Sunday, May 22, 2011, Joplin
#' Missouri, USA, was struck by a destructive tornado resulting in over 150
#' deaths and $2.8 billion in damages.
"Tornadoes2011MO"

#' @rdname Datasets
#' @description
#'
#' \code{ContinentalUSA_state_names}
#'
#' In several examples we use the names
#' of the USA states. This data is contained
#' in the object \code{ContinentalUSA_state_names}.
"ContinentalUSA_state_names"

#' @rdname Datasets
#' @description
#'
#' \code{MOAggIncomeLevelsPerCounty}
#'
#' In some examples we use the aggregate income levels
#' per county for the state of Missouri, USA.
#' This data is contained in the
#' object \code{MOAggIncomeLevelsPerCounty}.
"MOAggIncomeLevelsPerCounty"

#' @rdname Datasets
#' @description
#'
#' \code{USAStatesCounties2016}
#'
#' In many examples we use the boundaries of the states and counties of the USA.
#' The Cartographic Boundary Shapefiles (boundary data) is provided by the
#' USA Census Bureau at \url{https://www.census.gov/geo/maps-data/data/tiger-cart-boundary.html}.
#' This is a list containing \code{StateNames} (the state and territory names, 56 total),
#' \code{StatePolygons} (a list of size 56, containing a list of matrices describing the boundaries of the states/territories),
#' and \code{CountiesbyState} (a list of size 56, with each list element a list containing element
#' \code{CountyName} and \code{CountyPolies}, describing the county name of the specific state and polygons used).
#' For example, USAStatesCounties2016$StateNames[1] is Alabama, and
#' USAStatesCounties2016$CountiesbyState[[1]]$CountyName[1] corresponds to Escambia county
#' which can be plotted using the boundary coordinates in USAStatesCounties2016$CountiesbyState[[1]]$CountyPolies[[1]].
"USAStatesCounties2016"
