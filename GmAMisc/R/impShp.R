#' R function to easily import a vectorial dataset (shapefile) into R
#'
#' The function is a wrapper for the 'shapefile()' function out of the 'raster' package. It provides
#' the facility to import
#' a vectorial dataset (of shapefile type) by means of a window that allows the user to navigate
#' through the computer's folders and to select the appropriate file.
#'
#' @keywords impShp
#'
#' @export
#'
#' @importFrom raster shapefile
#'
#' @examples
#' \dontrun{
#' #a window will pop up allowing the user to select the shapefile
#' my.shapefile <- impShp()
#'}
#'
impShp <- function (){
  my.shapef <- raster::shapefile(file.choose())
  return(my.shapef)
}
