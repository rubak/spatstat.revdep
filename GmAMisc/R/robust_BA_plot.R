#' R function to plot a robust version of the Bland-Altman plot
#'
#' The function allows to plot a robust version of the Bland-Altman plot.\cr
#'
#' The function returns a chart based on robust (i.e. resistant to outlying values) measures of
#' central tendency and variability: median and Median Absolute Deviation (MAD) (Wilcox R R. 2001.
#' "Fundamentals of modern statistical methods: Substantially improving power and accuracy". New
#' York: Springer) instead of mean and standard deviation.\cr
#'
#' The x-axis displays the median of the
#' two variables being compared, while the y-axis represents their difference. A solid horizontal
#' line represents the bias, i.e. the median of the differences reported on the y-axis. Two dashed
#' horizontal lines represent the region in which 95percent of the observations are expected to lie;
#' they are set at the median plus or minus z*(MAD/0.6745).
#'
#' @param a Vector storing the first set of measurements to be compared.
#' @param b Vector storing the second set of measurements to be compared.
#' @param z Value for the confidence interval for the median difference; set by default to 1.96
#'   (corresponding to 95 percent CI).
#' @param methAlab Label to be used in the returned chart to refer to the method yielding the
#'   first set of measurements.
#' @param methBlab Label to be used in the returned chart to refer to the method yielding the
#'   second set of measurements.
#' @param cex Size of the data points.
#'
#' @keywords robustBAplot
#'
#' @export
#'
#' @examples
#' #create a first toy vector
#' a <- rnorm(30,10,1)
#'
#' #create a second toy vector
#' b <- a*runif(30,1,1.5)
#'
#' robustBAplot(a,b)
#'
robustBAplot <- function(a,b,z=1.96,methAlab="Method A", methBlab="Method B",cex=0.60){
  x <- apply(rbind(a,b),2,median)
  y <- a-b
  median_diff <- median(y)
  median_diff_abs_dev <- median(abs(median_diff-y))
  lower_l <- median_diff-z*(median_diff_abs_dev/0.6745)
  upper_l <- median_diff+z*(median_diff_abs_dev/0.6745)
  max_y_lim <- max(y,upper_l)+5
  min_y_lim <- min(y,lower_l)-5
  graphics::plot(x,y, pch=19, cex=cex, ylim=c(min_y_lim, max_y_lim), xlab=paste("Median of", methAlab, "and", methBlab), ylab=paste(methAlab, "-", methBlab))
  abline(h=0, lty=3)
  abline(h=median_diff, lty=1)
  abline(h=lower_l, lty=2)
  abline(h=upper_l, lty=2)
  title(main=paste("Robust Bland-Altman plot","\nMedian (bias)=", round(median_diff,3), ", Median Absolute Deviation=", round(median_diff_abs_dev,3),"\nLower limit=", round(lower_l,3) ,", Upper limit=", round(upper_l, 3), "\nlimits set at: median +/-", z, "* (MAD/0.6745)"), cex.main=0.8)
}
