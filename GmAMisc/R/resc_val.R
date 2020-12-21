#' R function to rescale the values of a dataset between a minimum and a maximum set by the user
#'
#' The function allows to rescale the values of a dataset between a minimum and a maximum that are
#' specified by the user. In doing that, it allows to preserve the shape of the distribution of the
#' original data.
#'
#' The function produces two density charts representing the distribution of the original and of the
#' rescaled dataset. It also returns a dataframe storing the original and rescaled values.
#'
#' @param x Vector containing the values to be rescaled.
#' @param min.v Minimum value of the rescaled dataset (0 by default).
#' @param max.v Maximum value of the rescaled dataset (100 by default).
#'
#' @keywords resc.val
#'
#' @export
#'
#' @examples
#' #generate a random dataset of size 30, normally distributed with mean 1000 and
#' #standard deviation 10
#' dataset <- rnorm(30, 1000,10)
#'
#' #rescale the dataset to be constrained between 10 and 100
#' resc.val(dataset, min.v=10, max.v=100)
#'
resc.val <- function (x, min.v=0, max.v=100) {
  resc.v <- ((x-min(x))*(max.v-min.v)/(max(x)-min(x)))+min.v

  par(mfrow=c(2,1))

  d.orig <- stats::density(x)
  graphics::plot(d.orig, xlab = "",
       main = "Density distribution of the original data",
       cex.main = 0.9)
  graphics::polygon(d.orig, col = "#BCD2EE88", border = "blue")
  rug(x, col="#0000FF") #hex code for 'blue'; last two digits set the transparency

  d.resc <- stats::density(resc.v)
  graphics::plot(d.resc, xlab = "",
       main = "Density distribution of the rescaled data",
       cex.main = 0.9)
  graphics::polygon(d.resc, col = "#BCD2EE88", border = "blue")
  rug(resc.v, col="#0000FF") #hex code for 'blue'; last two digits set the transparency

  df <- data.frame(original.data=x, rescaled.data=resc.v)

  #restore the default graphic settings
  par(mfrow = c(1,1))

  return(df)
}
