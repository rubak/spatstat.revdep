#' Computes the sine of an angle expressed in decimal degrees.
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
#' sind(90)
#' plot(seq(0,360,1), sind(seq(0,360,1)),type='l')

sind<-function(degrees) {
  sine<-sin(degrees*pi/180)
  return(sine)
}
