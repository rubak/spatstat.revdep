#' Plot Calibration Simplex
#'
#' @param x Object of class \code{calibration_simplex}
#' @param true_error Logical, specifying whether to use true miscalibration errors or approximate miscalibration errors.
#' @param error_scale A number specifying the magnitude of the miscalibration errors (greater 0, usually should be less than 1,
#' cf. note below).
#' @param min_bin_freq A number. Lower bound for (absolute) frequencies, i.e. how many observations have to lie in a bin
#' for it to be plotted.
#' @param plot_error_scale Logical, specifying whether to plot a scale showing the magnitude of miscalibration errors.
#' @param scale_area Optional. A number by which the areas of the points are scaled. Use if points are to small or to big.
#' @param indicate_bins Logical, specifying whether to connect points to their respective bin (center of hexagon).
#' @param category_labels A vector of length 3 containing the category names, e.g. \code{c("1","2","3")} (default)
#' @param use_pvals Logical, determines whether multinomial p-values are used for uncertainty quantification, see details.
#' @param alphas Vector of length 2 with values 1 > \code{alphas[1]} > \code{alphas[2]} >= 0.0001. Only relevant if \code{use_pvals = TRUE}. 
#' @param ... Arguments concerning the title (e.g. \code{main}, \code{cex.main}, \code{col.main} and \code{font.main})
#'            and subtitle (e.g. \code{sub}, \code{cex.sub}, \code{col.sub} and \code{font.sub}) may be passed here.
#' @details If multinomial p-values are used (\code{use_pvals = TRUE}), the dots are colored in the following way:
#' \itemize{
#'   \item Blue: p-value greater \code{alphas[1]} (0.1 by default).
#'   \item Orange: p-value between \code{alphas[1]} and \code{alphas[2]} (0.1 and 0.01 by default)
#'   \item Red: p-value less than \code{alphas[2]} (0.01 by default)
#'   \item Black: p-value is exactly 0. This only happens if a category which is assigned 0 probability realizes.
#' }
#' Many small p-values (orange and red dots) indicate miscalibrated predictions, whereas many blue dots indicate that the predictions 
#' may in fact be calibrated. WARNING: The use of the multinomial p-values is more of an experimental feature and may not yield reliable 
#' p-values, especially if \code{n} is small.
#' For details regarding the calculation of the p-values see also \code{\link{calibration_simplex}}.
#' 
#' @note For details on the meaning of the error scale, cf. Wilks, 2013, especially Fig. 2. Note that the miscalibration error in
#' each category is in "probability units" (as it is the average difference in relative frequency and forecast probability
#' in each bin).
#'
#'
#' @rdname plot.calibration_simplex
#' @export
#'
#' @importFrom graphics arrows axis par plot segments symbols text title
#' @importFrom spatstat.geom coords hexgrid hextess owin square

plot.calibration_simplex = function(x,
                                    #alpha = 0.05,
                                    
                                    true_error = TRUE,
                                    error_scale = 0.3,
                                    min_bin_freq = 10,
                                    plot_error_scale = TRUE,
                                    scale_area = NULL,
                                    indicate_bins = TRUE,
                                    category_labels = c("1","2","3"),
                                    use_pvals = FALSE,
                                    alphas = c(0.1,0.01),
                                    ...) {
  par_old <- par(no.readonly = TRUE)
  on.exit(par(par_old))

  n = x$n
  n_bins = x$n_bins
  error = error(x,true_error)
  rel_freq = rel_freq(x)

  if(max(x$freq)<min_bin_freq) stop("Nothing to plot here. Try reducing min_bin_freq (default = 10).")

  if(is.null(scale_area)) scale_area = (n/10)^2

  # Testing

  # tests = x$pvals > alpha
  # test_col = ifelse(tests,"blue","red")
  # test_col[x$pvals == -1] = "violet"

  #tests = x$pvals > alpha
  test_col = ifelse(x$pvals >= alphas[1],"blue",ifelse(x$pvals >= alphas[2],"orange","red"))
  test_col[x$pvals == -1] = "black"
  test_col[is.na(test_col)] = "black"

  # Plotting

  triangle= function(centers = F) {
    a = n - 1
    if(a == 0) a= 0.01 # allows to plot one bin (n = 1)
    s = 1/sqrt(3)
    sin60 = sqrt(3)/2
    x = c(0,sin60*a,0)
    y = c(0,a/2,a)
    nodes = data.frame(x,y)
    polygon = owin(poly=nodes)
    if(centers) return(hexgrid(polygon,s,trim=FALSE,origin=c(0,0)))
    else return(hextess(polygon,s,trim=FALSE,origin=c(0,0)))
  }

  H = triangle()
  H_centers = triangle(T)
  par(pty='s',mar=c(2,0,2,0),fig=c(0,1,0,1),cex.main = 1.5)
  plot(c(-2,n),c(-1,n),col="white",asp=1,bty="n",axes = F,xlab='',ylab='')
  if(is.element("main",names(list(...)))) title(...,line =-1)
  else title(main = "Calibration Simplex", line = -1,...)
  plot(H,border="darkgrey",add = T)

  centers = coords(H_centers)
  ordered_centers = centers[order(centers$x,centers$y),]

  displacement = cbind(-error[,3]-error[,1],0.577*error[,3]-0.577*error[,1])
  # displacement = cbind(-error$c3-error$c1,
  #                      0.577*error$c3-0.577*error$c1)

  shifted_centers = ordered_centers + 1/sqrt(3)/error_scale*displacement

  if(use_pvals){
    symbols(shifted_centers[x$freq>=min_bin_freq,],
            circles=sqrt(scale_area*rel_freq[x$freq>=min_bin_freq]/pi),
            inches = F,ann=F,fg = test_col[x$freq>=min_bin_freq],bg = test_col[x$freq>=min_bin_freq],add=T)
    if(indicate_bins) segments(shifted_centers$x[x$freq>=min_bin_freq],
                               shifted_centers$y[x$freq>=min_bin_freq],
                               ordered_centers$x[x$freq>=min_bin_freq],
                               ordered_centers$y[x$freq>=min_bin_freq],
                               col=test_col[x$freq>=min_bin_freq])
  }
  else{
    symbols(shifted_centers[x$freq>=min_bin_freq,],
            circles=sqrt(scale_area*rel_freq[x$freq>=min_bin_freq]/pi),
            inches = F,ann=F,bg = "black",add=T)
    if(indicate_bins) segments(shifted_centers$x[x$freq>=min_bin_freq],
                               shifted_centers$y[x$freq>=min_bin_freq],
                               ordered_centers$x[x$freq>=min_bin_freq],
                               ordered_centers$y[x$freq>=min_bin_freq],
                               col='red')
  }

  centers_extremes = ordered_centers[c(1,n,n_bins),]
  c1 = 3/4
  c2 = c1 *sqrt(3)
  shift_start = cbind(c(-c2,c2,0),c(c1,c1,-2*c1))
  shift_end = -shift_start[c(2,3,1),]
  start = centers_extremes+shift_start
  end = centers_extremes[c(2,3,1),]+shift_end
  arrows(start$x,start$y,end$x,end$y,length = 0.1)
  text((start+end+0.5*(shift_start+shift_end))/2,
       labels=as.expression(c(bquote(p[.(category_labels[3])]),
                                                                        bquote(p[.(category_labels[2])]),
                                                                        bquote(p[.(category_labels[1])]))),
       adj = 0.2)
  label_coords = rbind(centers_extremes+0.67*shift_end,end-0.33*shift_end)
  rot = c(30,0,-30)
  for(i in 1:3) {
    text(label_coords$x[c(i,i+3)],label_coords$y[c(i,i+3)],
         labels = paste(c(0,n-1),'/',n-1,sep=""),
         cex = 0.8,srt=rot[i])
  }

  if(plot_error_scale) {
    par(fig = c(0.65,0.95,0.05,0.35),new=T,cex.main=0.8,cex.axis=0.8,mgp=c(3,0.4,0))
    L = square(0.1)
    hexagon = hextess(L,1,trim=F,origin = c(0.0,0.0))
    plot(c(-1,1),c(-2,1),col="white",asp=1,main = "Error Scale",
         bty="o",xaxt="n",
         yaxt="n", xlab='',ylab='')
    plot(hexagon,add=T)
    axis(side = 1,at=c(-1,0,1),pos=-1.1,labels = paste(c(-error_scale,0,error_scale)))
  }
}
