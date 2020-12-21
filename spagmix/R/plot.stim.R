plot.stim <- function(x, fix.range = FALSE, sleep = 0.2, override.par = TRUE, ...){
  ellip <- list(...)
  if(is.null(ellip)) ellip <- list()
  if(is.null(ellip$box)) ellip$box <- FALSE
  if(is.null(ellip$ribargs)) ellip$ribargs <- list(box=TRUE)
  if(is.null(ellip$log)) ellip$log <- FALSE
  mn <- is.null(ellip$main)

  if(override.par) par(mfrow=c(1,1),mar=rep(2,4))

  lst <- x$v
  zlimeq <- c(0,min(sapply(lst,max)[sapply(lst,max)>0]))
  zlimconstant <- range(sapply(lst,range))
  if(ellip$log&&fix.range) zlimconstant <- log(zlimconstant)

  grt <- x$tlay

  for(i in 1:length(lst)){
    dev.hold()
    ellip$x <- lst[[i]]
    if(mn) ellip$main <- paste("t =",round(grt[i],5))
    if(diff(range(lst[[i]]))==0&&is.null(ellip$zlim)&&!fix.range) ellip$zlim <- zlimeq
    if(fix.range) ellip$zlim <- zlimconstant
    do.call("plot.im",ellip)
    if(!is.null(x$W)) plot(x$W,add=TRUE)
    axis(1)
    axis(2)
    box(bty="l")
    dev.flush()
    Sys.sleep(sleep)
  }
  invisible(NULL)

}
