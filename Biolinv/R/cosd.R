#' Computes the cosine of an angle expressed in decimal degrees.
#'
#' @param degrees either a number or a numeric vector of angles in decimal degrees.
#'
#' @return number or numeric vector.
#'
#' @author Luca Butikofer
#'
#' @export
#'
#' @examples
#' cosd(90)
#' plot(seq(0,360,1), cosd(seq(0,360,1)), type='l')

cosd<-function(degrees) {
  sine<-cos(degrees*pi/180)
  return(sine)
}
