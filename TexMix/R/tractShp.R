#' @title Polygon layer of Dallas County, TX, census tracts
#' @description Census tracts in Dallas County, Texas, in the longitude and
#'   latitude format (see \code{proj4string=CRS("+proj=longlat +ellps=WGS84")}).
#'
#'   \strong{Attention}: Two census tracts do not have a night time population and therefore
#'   NA's in most of their variables. These tracts are the airports of Love Field  and DFW
#'   (see \code{tractShp$TRACT=="48113980000" | tractShp$TRACT=="48113980100"}).
#' @docType data
#' @name tractShp
#' @source Based on 2018 ACS data, which were retrieved from Maptitude
#' (\url{https://www.caliper.com/}). The store density statistics by census
#' tract were calculated by the authors.
#' @format Spatial polygon data-frame with 529 census tracts. The variables are
#' as follows:
#' \describe{
#'   \item{ID}{Internal ID.}
#'   \item{SeqId}{Sequence ID of tracts from 1 to 529.}
#'   \item{TRACT}{Factor with the Census Bureau's tract numbers.}
#'   \item{AREA}{Calculated area of the census tract in square miles. \strong{
#'               Proportional to the denominator of the average store weighted
#'               kernel density.}}
#'   \item{LANDAREA}{Land area of the census tract in square miles.}
#'   \item{WATERAREA}{Water area of the census tract in square miles.}
#'   \item{nCells}{Number of 100x100 meters raster cells within each census
#'                 tract. \strong{Denominator of the average store weighted kernel
#'                 density in each census tract.}}
#'   \item{all1500D}{Average weighted kernel density of all stores with a
#'                   bandwith of 1500 meters. Weights are based on store's food sales volume.}
#'   \item{all2250D}{Average weighted kernel density of all stores with a
#'                   bandwith of 2250 meters. Weights are based on store's food sales volume.}
#'   \item{all3000D}{Average weighted kernel density of all stores with a
#'                   bandwith of 3000 meters. Weights are based on store's food sales volume.}
#'   \item{good2000D}{Average weighted kernel density of grocery stores with a
#'                    bandwith of 2000 meters. Weights are based on store's food sales volume.}
#'   \item{good3000D}{Average weighted kernel density of grocery stores with a
#'                    bandwith of 3000 meters. Weights are based on store's food sales volume.}
#'   \item{good4000D}{Average weighted kernel density of grocery stores with a
#'                    bandwith of 4000 meters. Weights are based on store's food sales volume.}
#'   \item{bad1000D}{Average weighted kernel density of convenience stores with a
#'                    bandwith of 1000 meters. Weights are based on store's food sales volume.}
#'   \item{bad1500D}{Average weighted kernel density of convenience stores with a
#'                   bandwith of 1500 meters. Weights are based on store's food sales volume.}
#'   \item{bad2000D}{Average weighted kernel density of convenience stores with a
#'                   bandwith of 2000 meters. Weights are based on store's food sales volume.}
#'   \item{LRRlowD}{Average low Log-relative risk \eqn{log(bad1000D/good2000D)}.}
#'   \item{LRRmedD}{Average medium Log-relative risk \eqn{log(bad1500D/good3000D)}.}
#'   \item{LRRhiD}{Average high Log-relative risk \eqn{log(bad2000D/good4000D)}.}
#'   \item{DES3NEIG}{\strong{Factor} distinguishing the three putative food desert
#'                   neighborhood against the remaining census tracts.}
#'   \item{FOODDES}{\strong{Factor} distinguishing the 14 putative food desert
#'                  census tracts against the remaining 515 census tracts.}
#'   \item{CITYPERI}{\strong{Factor} distinguishing the census tracts in the city
#'                   of Dallas, the Park cities, North, East, South and West
#'                   census tracts.}
#'   \item{DAYPOP}{Caliper's estimate of the absolute day time population.}
#'   \item{NIGHTPOP}{Census's night time population counts in a census tract.}
#'   \item{PCTDAYPOP}{\% day time population: \eqn{DAYPOP/(DAYPOP+NIGHTPOP)}.}
#'   \item{POPDEN}{Population density: \eqn{(0.4*DAYPOP+0.6*NIGHTPOP)/LANDAREA}.
#'                 \strong{Relative measure for potential food demand.}}
#'   \item{BUYPOW}{Absolute buying power in a census tract in $. Source: IRS 2016
#'                 records. \strong{Absolute measure for potential food demand.}}
#'   \item{MALE}{Census's male population count in a census tract.}
#'   \item{FEMALE}{Census's female population count in a census tract.}
#'   \item{MEDAGE}{Median population age in a census tract.}
#'   \item{PCTWHITE}{\% white population in a census tract.}
#'   \item{PCTBLACK}{\% black population in a census tract.}
#'   \item{PCTASIAN}{\% asian population in a census tract.}
#'   \item{PCTHISPAN}{\% hispanic population in a census tract.}
#'   \item{PCTMINOR}{\% of population belonging to a minority.}
#'   \item{PCTBADENG}{\% of population, which does not speak English well.}
#'   \item{PCTNOHIGH}{\% of population 25+ without a high school degree.}
#'   \item{PCTUNIVDEG}{\% of population 25+ with a univeristy degrees.}
#'   \item{PCTNOVEH}{\% of households not owning a car.}
#'   \item{PCTPUB2WRK}{\% of employed population taking public transportation
#'                     to work.}
#'   \item{TIME2WORK}{Average travel time to work. \strong{Attention}: 110 NA's.}
#'   \item{PCTNOHINS}{\% of civilian population without health insurance.}
#'   \item{PCTUNEMP}{\% of population in the labor force, which is unemployed.}
#'   \item{PCTFAMPOV}{\% of family below the poverty threshold.}
#'   \item{PCTPOPPOV}{\% of population below the poverty threshold.}
#'   \item{HHMEDINC}{Median household income in $.}
#'   \item{MEDFAMINC}{Median family income in $.}
#'   \item{PERCAPINC}{Per capita income in $.}
#'   \item{PCTHUVAC}{\% of vacant housing units.}
#'   \item{MEDVALHOME}{Median home value. \strong{Attention}: 27 NA's.}
#'   \item{PCTB2010}{\% of homes built from 2010 to 2018.}
#'   \item{PCTB2000}{\% of homes built from 2000 to 2009.}
#'   \item{PCTB1990}{\% of homes built from 1990 to 1999.}
#'   \item{PCTB1980}{\% of homes built from 1980 to 1989.}
#'   \item{PCTB1970}{\% of homes built from 1970 to 1979.}
#'   \item{PCTB1960}{\% of homes built from 1960 to 1969.}
#'   \item{PCTB1950}{\% of homes built from 1950 to 1959.}
#'   \item{PCTB1940}{\% of homes built from 1940 to 1949.}
#'   \item{PCTBPRE}{\% of homes built before 1940.}
#' }
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
NULL
