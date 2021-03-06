#' Frogs sighting locations.
#'
#' Dataset containing sighting locations of L. raniformis in New Zealand.
#'
#' \itemize{
#'   \item year: year of the sighting
#'   \item y: latitude in the New Zealand Transverse Mercatore coordinate system (length unit of measure: meter).
#'   \item x: longitude in the New Zealand Transverse Mercatore coordinate system (length unit of measure: meter).
#'   \item species: the name of the species.
#' }
#'
#' @docType data
#' @keywords datasets
#' @name frogs
#' @usage data(frogs)
#' @format A dataframe with 194 rows and 4 columns
#' @source Selection of sightings locations form the Herpetofauna Database of New Zealand. Kindly provided by Benno Kappers, Department of Conservation.
NULL


#' Frogs sighting locations and the output of the EM() function.
#'
#' Version of 'frogs' Dataset containing sighting locations of L. raniformis in New Zealand as returned by function EM().
#'
#' \itemize{
#'   \item year: year of the sighting
#'   \item y: latitude in the New Zealand Transverse Mercatore coordinate system (length unit of measure: meter).
#'   \item x: longitude in the New Zealand Transverse Mercatore coordinate system (length unit of measure: meter).
#'   \item species: the name of the species.
#'   \item Pnat: probability of being of natural origin [0;1].
#'   \item Dist: distance from nearest point of natural origin or form anchor point.
#' }
#'
#' @docType data
#' @keywords datasets
#' @name frogsEM
#' @usage data(frogsEM)
#' @format A dataframe with 194 rows and 6 variables
NULL


#' Four Jackknife re-samplings done on dataset 'frogs' generated by function jackKnife().
#'
#' List of four Jackknife re-samplings of dataframe 'frogs' as returned by function jackKnife()(see examples in ?jackKnife).
#'
#' \itemize{
#'   \item year: year of the sighting
#'   \item y: latitude in the New Zealand Transverse Mercatore coordinate system (length unit of measure: meter).
#'   \item x: longitude in the New Zealand Transverse Mercatore coordinate system (length unit of measure: meter).
#'   \item species: the name of the species.
#' }
#'
#' @docType data
#' @keywords datasets
#' @name frogsJK
#' @usage data(frogsJK)
#' @format A list of four dataframes with 164 rows and 6 columns each.
NULL


#' Simulated datasets for comparison with 'frogs' dataset generated with function simulacro() (see examples in ?simulacro).
#'
#' List of eight lists (one per different Alpha value used) each containing ten dataframes of simulated sigthing locations built to simulate 'frogs' dataset.
#'
#' \itemize{
#'   \item year: year of the sighting
#'   \item y: latitude in the New Zealand Transverse Mercatore coordinate system (length unit of measure: meter).
#'   \item x: longitude in the New Zealand Transverse Mercatore coordinate system (length unit of measure: meter).
#'   \item species: the name of the species.
#'   \item Pnat: probability of being of natural origin [0;1].
#'   \item Dist: distance from nearest point of natural origin or form anchor point.
#' }
#'
#' @docType data
#' @keywords datasets
#' @name frogsLacro
#' @usage data(frogsLacro)
#' @format list of 8 lists of 10 dataframes each with 164 rows and 6 columns each.
NULL


#' Summary of the dissimilarity values of each of the dataframes in 'frogsLacro' with dataframe 'frogs'.
#'
#' Dataframe with dissimilarity values and respective Alpha value as obtained by function modSel() (see examples in ?modSel).
#'
#' \itemize{
#'   \item dissimilarity: dissimilarity value
#'   \item compAlpha: alpha value of the comparison dataset.
#' }
#'
#' @docType data
#' @keywords datasets
#' @name frogsSum
#' @usage data(frogsSum)
#' @format list of 8 lists of 10 dataframes each with 164 rows and 6 columns each.
NULL


#' Landmasses of the two main islands of New Zealand as SpatialPolygons.
#'
#' Object of class sp::SpatialPolygons. North and South Islands of New Zealand in the New Zealand Transverse Mercatore projection.
#'
#' @docType data
#' @keywords datasets
#' @name nzp
#' @usage data(nzp)
#' @format object of class sp::SpatialPolygons.
#' @source mapdata::worldHires()
NULL


#' Landmasses of the two main islands of New Zealand as Window object.
#'
#' Object of class spatstat.geom::owin of the North and South Islands of New Zealand in the New Zealand Transverse Mercatore projection obtained from object 'nzp' (see ?nzp).
#'
#' @docType data
#' @keywords datasets
#' @name nzw
#' @usage data(nzw)
#' @format object of class spatstat::owin.
NULL
