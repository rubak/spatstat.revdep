#
# loccit.R
#
# Local second-order composite likelihood
# and local Palm likelihood
# for cluster/Cox processes.
#
#  $Revision: 2.17 $ $Date: 2019/04/13 15:37:08 $
#

loccit <- local({

  # define objective functions for optimization
  
  PalmObjPar <- function(par, objargs) {
    with(objargs,
#         sumWJloglamJ +
         sum(wJ * log(paco(dIJ, par)))
         - Gscale * unlist(stieltjes(paco, Gv, par=par)),
         enclos=objargs$envir)
  }
  PalmObjTheta <- function(par, objargs) {
    with(objargs,
#         sumWJloglamJ +
         sum(wJ * log(paco(dIJ, par)))
         - Gscale * unlist(stieltjes(paco, Gv, theta=par)),
         enclos=objargs$envir)
  }

  # kernel function 
  wtfun <- function(x,y, xv, yv, sigma) {
    dnorm(x - xv, sd=sigma) * dnorm(y - yv, sd=sigma)
  }

  # list of formal arguments of loccit that are
  # not relevant to the 'template model'
  loccitparams <- c("sigma", "f", "clustargs", "control",
                    "rmax", "verbose")

  loccit <-
    function(X, trend = ~1, 
             clusters = c("Thomas","MatClust","Cauchy","VarGamma","LGCP"),
             covariates = NULL,
             ...,
             diagnostics=FALSE,
             taylor = FALSE,
             sigma=NULL, f=1/4,
             clustargs=list(),
             control=list(),
             rmax,
             covfunargs=NULL,
             use.gam=FALSE,
             nd=NULL, eps=NULL,
             niter = 3,
             fftopt = list(),
             verbose=TRUE
             ) {
  starttime <- proc.time()
  Xname <- short.deparse(substitute(X))
  cl <- match.call()
  stopifnot(is.ppp(X))
  if(is.marked(X))
    stop("Sorry, cannot handle marked point patterns")

  clusters <- match.arg(clusters)
  diagnostics <- identical(diagnostics, TRUE)
  if(diagnostics && !taylor)
    stop("Diagnostics are only implemented for the Taylor approximation")
  
  W <- as.owin(X)
  nX <- npoints(X)
  
  # make a fake call representing the template model
  templatecall <- cl[!(names(cl) %in% loccitparams)]
  templatecall[[1]] <- as.name("kppm")
  
  # determine rmax and bandwidth
  if(missing(rmax) || is.null(rmax))
    rmax <- rmax.rule("K", W, intensity(X))
  if(is.null(sigma)) {
    sigma <- bw.frac(X, f=f)
    if(verbose)
      cat(paste("sigma = ", sigma, "\n"))
  }
  
  # get pair correlation function (etc) for model
  info <- spatstatClusterModelInfo(clusters)
  pcfun      <- info$pcf
  funaux     <- info$funaux
  selfstart  <- info$selfstart
  isPCP      <- info$isPCP
  parhandler <- info$parhandler
  modelname  <- info$modelname
  # process additional parameters of cluster model
  clargs <- if(is.function(parhandler)) do.call(parhandler, clustargs) else NULL
  pcfunargs <- append(clargs, list(funaux=funaux))
  # determine starting parameter values and parameter dimension
  selfstartpar <- selfstart(X)
  npar <- length(unlist(selfstartpar))

  # get additional info
  extra <- extraClusterModelInfo(clusters)
  if(use.transform <- !is.null(extra)) {
    # use transformed parameters
    par2theta <- extra$par2theta
    theta2par <- extra$theta2par
    pcfun      <- extra$pcftheta
  }

  # ................................................................
  #  Step 0: homogeneous model
  # ................................................................
  if(verbose) cat("Fitting homogeneous cluster model ... ")
  homclusfit <- kppm(X,
                     trend=trend,
                     clusters=clusters,
                     covariates=covariates,
                     covfunargs=covfunargs,
                     use.gam=use.gam,
                     nd=nd, eps=eps,
                     method="palm")
  homcluspar <- homclusfit$par
  if(use.transform) 
    homclustheta <- applymaps(par2theta, homcluspar)
  fit0 <- as.ppm(homclusfit)
  Q <- quad.ppm(fit0)
  U <- union.quad(Q)
  Z <- is.data(Q)
  ntrend <- length(coef(fit0))
  if(verbose) cat("Done.\n")

  
  # ................................................................
  #  Step 1: estimate intensity
  # ................................................................
  # ................................................................
  #  fit local Poisson model
  # ................................................................
  if(verbose) cat("Fitting local Poisson model...")
  # compute weights 
  Swt <- strausscounts(U, X, rmax, EqualPairs=equalpairs.quad(Q))
  # fit local Poisson model
  locpoisfit <- locppm(Q,
                       trend=trend,
                       covariates=covariates, 
                       covfunargs=covfunargs,
                       locations="fine",
                       use.gam=use.gam,
                       weightfactor=Swt,
                       sigma=sigma, use.fft=taylor, verbose=FALSE)
  # fitted intensity at quadrature points
  locpoislambdaU <- as.vector(fitted(locpoisfit))
  # fitted local coefficients at quadrature points
  locpoiscoef <- marks(coef(locpoisfit))
  # predict intensity at all locations
  locpoislambda <-
    Smooth(ssf(U, locppmPredict(as.ppm(locpoisfit), locpoiscoef)))
  #
  if(verbose) cat("Done.\n")
  #
  # ..........................................................
  #    Step 2. Fit cluster parameters
  # ..........................................................
  
  # Fix intensity
  lambdaU <- locpoislambdaU
  lambda <- locpoislambda
  #
  
  if(taylor) {
    # ..........................................................
    #          first-order Taylor approximation
    # ..........................................................
    if(verbose) cat("Computing Taylor approximation... ")
#    what <- c("param", "lambda")
    what <- "param"
    if(diagnostics)
      what <- c(what, "influence", "plik", "score", "grad", "delta")
    FT <- loccitFFT(homclusfit, sigma, rmax,
                    base.trendcoef=locpoiscoef,
                    base.lambda=locpoislambdaU,
                    base.lambdaim=locpoislambda,
                    what=what, verbose=verbose,
                    calcopt=fftopt,
                    ...)
    if(verbose) cat("Done.\nExtracting values... ")
    # local cluster parameter estimates
    estpar <- sample.imagelist(FT$parameters, U)[, ntrend+1:npar]
    # fitted intensity
    # DO NOT UPDATE!
    # lambdaU <- sample.imagelist(FT$intensity, U)
    #
    if(diagnostics) {
      score   <- sample.imagelist(FT$score, U)
      hessian <- sample.imagelist(FT$gradient, U)
      delta <- sample.imagelist(FT$delta, U)
      # influence of each data point
      influ <- FT$influence
      # log Palm likelihood
      logplik <- FT$logplik
    }
    if(verbose) cat("Done.\n")
  } else {
    # ..........................................................
    #         Full fitting algorithm
    # ..........................................................
    Xx <- X$x
    Xy <- X$y
    Ux <- U$x
    Uy <- U$y
    nU <- npoints(U)
    ok <- getglmsubset(as.ppm(homclusfit))
    nV <- sum(ok)
    # starting parameter vector
    startpar <- if(use.transform) homclustheta else homcluspar
    # allocate space for results
    estpar <- matrix(NA_real_, nU, npar)
    colnames(estpar) <- names(startpar)
    # identify pairs of points that contribute
    cl <- closepairs(X, rmax)
    I <- cl$i
    J <- cl$j
    dIJ <- cl$d
    xI <- cl$xi
    yI <- cl$yi
    # 
    # create local function with additional parameters in its environment
    paco.par <- function(d, par) {
      do.call(pcfun, append(list(par=par, rvals=d), pcfunargs))
    }
    paco.theta <- function(d, theta) {
      do.call(pcfun, append(list(theta=theta, rvals=d), pcfunargs))
    }
    paco <- if(use.transform) paco.theta else paco.par
    
    # define objective function obj(par, objargs)
    obj <- if(use.transform) PalmObjTheta else PalmObjPar
    
    # collect data that doesn't change
    objargs0 <- list(paco=paco, rmax=rmax, dIJ=dIJ, 
                     envir=environment(paco))
    # set up control arguments for 'optim'
    ctrl <- resolve.defaults(list(fnscale=-1), control, list(trace=0))
    # 
    # ..............................................................
    #             fit local cluster model at each quadrature point
    # ..............................................................
    #
    # intensity at data points
    lambdaX <- lambdaU[Z]
    # total intensity
    Lam <- integral.im(lambda)
    warnlist <- NULL
    if(verbose)
      cat(paste0("Maximising local composite likelihood on ", nV,
                 if(nV == nU) NULL else paste("out of ", nU),
                 " quadrature points..."))
    jok <- which(ok)
    # ...................
    for(k in seq_len(nV)) {
      if(verbose) progressreport(k, nV)
      j <- jok[k]
      xv <- Ux[j]
      yv <- Uy[j]
      wX <- wtfun(Xx, Xy, xv, yv, sigma)
      wI <- wX[I]
      wJ <- wX[J]
      wtf <- as.im(wtfun, W=lambda, xv=xv, yv=yv, sigma=sigma)
      wlam <- eval.im(wtf * lambda)
      Mv <- integral.im(wlam)
      # compute cdf of distance between 
      #   random point in W with density proportional to 'wlam'
      # and
      #   random point in X with equal probability
      #
      Gv <- distcdf(W, X, dW=wlam, dV=1)
      Gscale <- nX * Mv
      # trim Gv to [0, rmax] 
      Gv <- Gv[with(Gv, .x) <= rmax,]
      # pack up necessary information
      objargs <- append(objargs0,
                        list(wJ=wJ,
 #                                    sumWJloglamI = sum(wJ * loglambdaJ),
                             Gv=Gv,
                             Gscale=Gscale))
      # optimize it
      opt <- optim(startpar, obj, objargs=objargs, control=ctrl)
      # save warnings; stop if error
      warnlist <- accumulateStatusList(optimStatus(opt), warnlist) 
      # save fitted parameters
      estpar[j, ] <- unlist(opt$par)
    }
    printStatusList(warnlist)
    # transform back to original parameters
    if(use.transform)
      estpar <- applymaps(theta2par, estpar)
  }
  
  # infer parameter 'mu'
  if(isPCP) {
    # Poisson cluster process: extract parent intensity kappa
    kappaU <- estpar[, "kappa"]
    # mu = mean cluster size
    mu <- lambdaU/kappaU
  } else {
    # LGCP: extract variance parameter sigma2
    sigma2U <- estpar[, "sigma2"]
    # mu = mean of log intensity 
    mu <- log(lambdaU) - sigma2U/2 
  }
  modelpar <- cbind(estpar, mu=as.vector(mu))
  # pack up
  result <-
    list(homclusfit = homclusfit,
         locpoisfit = locpoisfit,
         lambda     = ssf(U, lambdaU),
         modelpar   = ssf(U, modelpar),
         coefs      = coef(locpoisfit),
         score      = if(taylor && diagnostics) ssf(U, score) else NULL,
         hessian    = if(taylor && diagnostics) ssf(U, hessian) else NULL,
         delta      = if(taylor && diagnostics) ssf(U, delta) else NULL,
         influence  = if(taylor && diagnostics) influ else NULL,
         logplik    = if(taylor && diagnostics) logplik else NULL,
         sigma      = sigma,
         templatestring=format(templatecall))
  class(result) <- c("loccit", class(result))
  result <- timed(result, starttime=starttime)
  return(result)
}

  loccit
})


plot.loccit <- function(x, ...,
                        what=c("modelpar", "coefs", "lambda"),
                        how=c("smoothed", "exact"),
                        which=NULL, pre=NULL, post=NULL) {
  xname <- deparse(substitute(x))
  if(!any(what %in% c("score", "hessian")))
     what <- match.arg(what)
  how <- match.arg(how)
  Z <- x[[what]]
  if(!is.null(which))
    Z <- Z[, which]
  if(how == "smoothed") {
    Z <- applymaps(pre, Z)
    Z <- Smooth(Z, ...)
    Z <- applymaps(post, Z)
  }
  do.call("plot", resolve.defaults(list(x=Z),
                                   list(...),
                                   list(main=xname)))
}

print.loccit <- function(x, ...) {
  cat("Local Palm likelihood fit\n")
  cat("Template model:\n")
  cat(paste("\t", x$templatestring, "\n\n"))
  cat(paste("Smoothing parameter sigma: ", x$sigma, "\n\n"))
  cat("Coefficient estimates:\n")
  print(x$coefs, brief=TRUE)
  return(invisible(NULL))
}

predict.loccit <- function(object, ...)  {
  predict(object$locpoisfit, ...)
}

fitted.loccit <- function(object, ..., new.coef=NULL) {
  trap.extra.arguments(...)
  lam <- predict(object, new.coef=new.coef)
  return(marks(lam))
}

with.loccit <- function(data, ...) {
  with(data$coefs, ...)
}

bw.loccit <- function(...,
                      use.fft=TRUE,
                      srange = NULL, ns=9, sigma = NULL,
                      fftopt=list(),
                      verbose=TRUE) {
  starttime <- proc.time()
  parenv <- sys.parent()

  if(!use.fft) stop("Sorry, not yet implemented when use.fft=FALSE")
  
  # fit homogeneous model
  homclusfit <- eval(substitute(kppm(..., method="palm", forcefit=TRUE)),
                     envir=parenv)
  ncluspar <- length(homclusfit$par)

  # extract quadrature info
  hompoisfit <- as.ppm(homclusfit)
  ntrendcoef <- length(coef(hompoisfit))
  X <- data.ppm(hompoisfit)
  Q <- quad.ppm(hompoisfit)
  U <- union.quad(Q)
  wQ <- w.quad(Q)
  Z <- is.data(Q)
  nX <- npoints(X)
  nU <- npoints(U)
  Xindex <- seq_len(nU)[Z]

  p <- ntrendcoef + ncluspar
  
  # determine values of smoothing parameter to be assessed
  if(!missing(sigma) && !is.null(sigma)) {
    stopifnot(is.numeric(sigma))
    ns <- length(sigma)
  } else {
    if(is.null(srange)) {
      srange <- range(bw.diggle(X) * c(1/2, 4),
                      bw.frac(X, f=1/3))
    } else check.range(srange)
    sigma <- exp(seq(log(srange[1]), log(srange[2]), length=ns))
  } 

  # fit using each value of sigma
  if(verbose) cat(paste("Assessing", ns, "values of sigma... "))
  dof <- edof <- logL <- numeric(ns)
  for(k in 1:ns) {
    if(verbose) progressreport(k, ns)
    sigk <- sigma[k]
    fitk <- eval(substitute(loccit(..., sigma=sigk,
                                   taylor=TRUE,
                                   diagnostics=TRUE,
                                   fftopt=fftopt,
                                   verbose=FALSE)),
                 envir=parenv)
    # log Palm likelihood of fit
    logL[k] <- fitk[[ "logplik" ]]
    # influence calculation
    influ <- fitk[[ "influence" ]]
    score.data <- influ$score.data
    delta.score <- influ$delta.score
    hQ <- influ$HP
    # Hessian at data points
#    hQ <- as.matrix(fitk[[ "hessian" ]])
    hX <- hQ[Z,,drop=FALSE]
    # invert the Hessians
    vX <- invert.flat.matrix(hX, p)
    # compute 'leverage' approximation 
    #    log(lambda(x_i)) - log(lambda_{-i}(x_i))
    #    ~ Z(x_i) V(x_i) Y(x_i)^T
    leve <- bilinear3.flat.matrices(score.data, vX, delta.score,
                                    c(1,p), c(p, p), c(1, p))
    dof[k] <- sum(leve)
    # 'expected leverage'..
    eleve <- bilinear3.flat.matrices(score.data, vX, score.data,
                                    c(1,p), c(p, p), c(1, p))
    kernel0 <- 1/(2 * pi * sigk^2)
    edof[k] <- kernel0 * sum(eleve)
  }
  # compute cross-validation criterion
  gcv <- logL - dof
  # optimize
  result <- bw.optim(gcv, sigma, iopt=which.max(gcv), cvname="cv",
                     dof=dof, logL=logL, edof=edof)
  timed(result, starttime=starttime)
}


psib.loccit <- function(object) {
  clus <- object$homclusfit$clusters
  info <- spatstatClusterModelInfo(clus)
  if(!info$isPCP)
    stop("The model is not a cluster process")
  xtra <- extraClusterModelInfo(clus)
  pmap <- xtra$psib
  if(is.null(pmap))
    stop(paste("Not implemented for", sQuote(clus)))
  modpar <- object$modelpar
  p <- applymaps(pmap, marks(modpar))
  P <- ssf(unmark(modpar), p)
  return(P)
}

accumulateStatusList <- function(x, stats=NULL, stoponerror=TRUE) {
  if(inherits(x, "error") && stoponerror)
    stop(conditionMessage(x))
  if(is.null(stats))
    stats <- list(values=list(), frequencies=integer(0))
  if(!inherits(x, c("error", "warning", "message")))
    return(stats)
  same <- unlist(lapply(stats$values, identical, y=x))
  if(any(same)) {
    i <- min(which(same))
    stats$frequencies[i] <- stats$frequencies[i] + 1
  } else {
    stats$values <- append(stats$values, list(x))
    stats$frequencies <- c(stats$frequencies, 1)
  }
  return(stats)
}

