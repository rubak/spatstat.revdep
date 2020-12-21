#
#  locmincon.R
#
#  Local minimum-contrast estimation
#
#   $Revision: 1.10 $  $Date: 2014/06/28 09:44:38 $
#

locmincon <- function(...,
                      sigma=NULL, f = 1/4, verbose=TRUE,
                      localstatargs=list(), LocalStats=NULL,
                      tau=NULL) {
  starttime <- proc.time()
  if(verbose) cat("Fitting template model...")
  kfit <- kppm(...)
  if(verbose) cat("Done.\n")

  argh <- list(...)
  
  X <- kfit$X
  nX <- npoints(X)
  pfit <- as.ppm(kfit)

  if(is.null(getglmfit(pfit)))
    pfit <- update(pfit, forcefit=TRUE)
  
  # determine bandwidth
  if(is.null(sigma)) {
    sigma <- bw.frac(X, f=f)
    if(verbose)
      cat(paste("sigma = ", sigma, "\n"))
  }

  stationary <- kfit$stationary   # logical
  isPCP      <- kfit$isPCP
  homlambda  <- kfit$lambda
  hompar     <- kfit$par 

  kFIT <- kfit$Fit
  HomStat <- kFIT$Stat # e.g. K-function object
  FitFun  <- kFIT$FitFun         # e.g. "thomas.estK"
  StatFun <- kFIT$StatFun         # e.g. "Kest"

  localmap <- c(Kest="localK",
                Kinhom="localKinhom",
                pcf="localpcf",
                pcfinhom="localpcfinhom")
  if(!(StatFun %in% names(localmap)))
    stop(paste("A local version of the function", sQuote(StatFun),
               "is not available"))
  LocalStatFun <- localmap[[StatFun]]

  # Compute predicted intensity of local model at each data point 
  if(is.stationary(pfit)) {
    # shortcut: equivalent to kernel smoothing
    if(verbose) cat("Estimating intensity by kernel smoothing.\n")
    lambda <- density(X, sigma=sigma, at="points")
  } else {
    # fit local Poisson models
    if(verbose) cat("Fitting local Poisson models.\n")
    coefs <- locppmEngine(pfit, sigma, X,
                          verbose=verbose, Vname="data points")
    # Compute predicted intensity of local model at each data point 
    gfit <- getglmfit(pfit)
    gdat <- getglmdata(pfit)
    lambda <- numeric(nX)
    for(i in seq_len(nX)) 
      lambda[i] <- GLMpredict(gfit, gdat[i,, drop=FALSE], coefs[i,],
                              changecoef=TRUE)
  }

  # compute local statistics
  if(is.null(LocalStats)) {
    if(verbose)
      cat(paste("Computing local statistics for", nX, "points..."))
    argh <- append(if(stationary) list(X) else list(X, lambda=lambda),
                   localstatargs)
    LocalStats <- do.call(LocalStatFun, argh)
    if(verbose) cat("Done.\n")
  } else {
    # validate
    stopifnot(is.fv(LocalStats))
    if(ncol(LocalStats) != nX + 2)
      stop(paste("Argument LocalStats has",
                 ncol(LocalStats), "columns",
                 "but", nX + 2, "columns were expected"))
    if(length(localstatargs) > 0)
      warning("Argument localstatargs ignored, because LocalStats was given")
  }
  # Allocate space to store all the weighted statistics
  WeightedStatistics <- LocalStats
  # It is assumed that the first nX columns contain the local functions.
  
  # Derive a 'template' fv object for the local weighted average statistic
  # at the current point, for minimum contrast fitting
  whichr <- which(fvnames(LocalStats, ".x") == colnames(LocalStats))
  fvnames(LocalStats, ".y") <- colnames(LocalStats)[1]
  WeightedStatJ <- LocalStats[, c(1, whichr)]
  # make it acceptable to FitFun
  attr(WeightedStatJ, "fname") <- attr(HomStat, "fname")
  
  # Extract matrix of local statistics
  LocalMatrix <- as.matrix(as.data.frame(LocalStats)[ , 1:nX])

  # detect NA or Inf entries, record their positions, replace by zeroes
  LocalMatrixValid <- is.finite(LocalMatrix)
  LocalMatrix[!LocalMatrixValid] <- 0
    
  # save r values
  rvals <- with(WeightedStatJ, .x)
   
  # Now compute local fits
  xx <- X$x
  yy <- X$y

  if(verbose) 
    cat(paste("Fitting", nX, "models...\n"))

  # set up space for results
  pars <- matrix(NA_real_, nrow=nX, ncol=3)
  colnames(pars) <- c(names(hompar), "mu")
  
  for(j in seq_len(nX)) {
    if(verbose)
      progressreport(j, nX)
    # compute vector of weights relative to point j
    localwt <- dnorm(xx - xx[j], sd=sigma) * dnorm(yy - yy[j], sd=sigma)
    localwt <- matrix(localwt, ncol=1)
    # compute weighted average of statistics
    numer <- LocalMatrix %*% localwt
    denom <- LocalMatrixValid %*% localwt
    weightedstat.j <- as.vector(numer/denom)
    WeightedStatJ[,1] <- weightedstat.j
    WeightedStatistics[, j] <- weightedstat.j
    # trim the recommended range of r values
    attr(WeightedStatJ, "alim") <- c(0, max(rvals[denom > 0]))
    # apply minimum contrast
    mcfitj <- do.call.matched(FitFun,
                              resolve.defaults(
                                               list(X=WeightedStatJ),
                                               argh,
                                               list(startpar=hompar)),
                              extrargs=names(formals(optim)))
    # extract fitted parameters
    parj <- mcfitj$par
    lambdaj <- lambda[j]
    if(isPCP) {
      kappaj <- parj[1]
      muj <- lambdaj/kappaj
    } else {
      sigma2j <- parj[1]
      muj <- lambdaj - sigma2j/2
    }
    pars[j, ] <- c(parj, muj)
  }
  if(verbose) cat("Done.\n")
  # assemble results
  ok <- complete.cases(pars)
  smoo <- Smooth(X[ok] %mark% pars[ok,], sigma=tau)
  tau <- attr(smoo, "sigma")
  # pack up
  result <- list(homfit=kfit, pars=pars, ok=ok, smoo=smoo,
                 localstats=LocalStats, weightedstats=WeightedStatistics,
                 sigma=sigma, tau=tau)
  class(result) <- c("locmincon", class(result))
  result <- timed(result, starttime=starttime)
  return(result)
}

plot.locmincon <- function(x, ...,
                           how=c("exact", "smoothed"),
                           which = NULL,
                           sigma=NULL,
                           do.points=TRUE) {
  xname <- deparse(substitute(x))
  how <- match.arg(how)
  X <- x$homfit$X
  switch(how,
           smoothed = {
             if(is.null(sigma)) {
               objects <- x$smoo
             } else {
               pars <- x$pars
               ok   <- x$ok 
               objects <- Smooth(X[ok] %mark% pars[ok,], sigma=sigma)
             }
           },
           exact = {
             objects <- lapply(as.list(as.data.frame(x$pars)),
                               setmarks,
                               x=X)
             objects <- as.listof(objects)
           }
           )
  if(!is.null(which))
    objects <- objects[which]
  pe <- if(do.points) function(...) { plot(X, pch=".", add=TRUE) } else NULL
  do.call("plot", resolve.defaults(list(x=objects),
                                   list(...),
                                   list(main=xname, panel.end=pe)))
}

print.locmincon <- function(x, ...) {
  cat("Local minimum contrast fit\n")
  cat("....................................\n")
  cat("Template model:\n")
  print(x$homfit)
  cat("....................................\n")
  cat(paste("\nLocalisation bandwidth sigma:", x$sigma, "\n"))
  cat(paste("\nInterpolation bandwidth tau:", x$tau, "\n"))
  cat("\nSmoothed coefficient estimates:\n ")
  print(x$smoo)
  return(invisible(NULL))
}

as.ppp.locmincon <- function(X, ...) { (X$homfit$X) %mark% (X$pars) }

Smooth.locmincon <- function(X, tau=NULL, ...) {
  # if sigma is missing, return the precomputed smooth images
  if(is.null(tau))
    return(X$smoo)
  # recompute with new sigma
  XX   <- X$homfit$X
  pars <- X$pars
  ok   <- X$ok
  Y <- XX[ok] %mark% pars[ok,]
  ans <- do.call("Smooth",
                 resolve.defaults(list(X=Y, sigma=tau),
                                  list(...)))
  return(ans)
}

with.locmincon <- function(data, ...) {
  result <- with(as.data.frame(data$pars), ...)
  X <- as.ppp(data)
  if((is.matrix(result) && nrow(result) == npoints(X)) ||
     (is.vector(result) && is.numeric(result) && length(result) == npoints(X)))
    result <- ssf(X, result)
  return(result)
}

     
psib.locmincon <- function(object) {
  clus <- object$homfit$clusters
  info <- spatstatClusterModelInfo(clus)
  if(!info$isPCP)
    stop("The model is not a cluster process")
  X <- data.ppm(as.ppm(object$homfit))
  xtra <- extraClusterModelInfo(clus)
  pmap <- xtra$psib
  if(is.null(pmap))
    stop(paste("Not implemented for", sQuote(clus)))
  p <- applymaps(pmap, object$pars)
  P <- ssf(X, p)
  return(P)
}
