#
#  locppm.R
#
# Local (pseudo)likelihood for point processes
#
#  $Revision: 1.183 $ $Date: 2016/10/11 10:55:06 $
#

locppm <- function(..., sigma=NULL, f = 1/4,
                   vcalc=c("none", "t", "hessian", "hom", "lik", "full"),
                   locations=c("split", "fine", "coarse"),
                   ngrid=NULL, grideps=NULL,
                   verbose=TRUE,
                   use.fft=FALSE, fft.algorithm="closepairs") {

  starttime <- proc.time()

  missloc <- missing(locations)
  locations <- match.arg(locations)
  vcalc <- match.arg(vcalc)

  if(missloc) {
    # Default is locations = "split"
    # except in the following cases
    if(vcalc == "none") locations <- "fine"
    if(vcalc %in% c("hom", "lik"))  locations <- "coarse"
  }

  stopifnot(is.logical(use.fft))
  
  # fit homogeneous model
  if(verbose) cat("Fitting homogeneous model... ")
  parenv <- sys.parent()
  homfit <- eval(substitute(ppm(..., forcefit=TRUE)), envir=parenv)
  if(is.multitype(homfit))
    stop("Sorry, cannot handle marked point processes yet")
  if(verbose) cat("Done.\n")
  templatecall <- format(substitute(ppm(...)))

  ispois <- is.poisson(homfit)

  X <- data.ppm(homfit)
  Q <- quad.ppm(homfit)
  U <- union.quad(Q)
  wU <- w.quad(Q)

  # determine smoothing parameter
  if(is.null(sigma)) {
    sigma <- bw.frac(X, f=f)
    if(verbose)
      cat(paste("sigma = ", sigma, "\n"))
  }

  # default grid dimensions
  need.grid <- (locations %in% c("coarse", "split"))
  if(need.grid && is.null(ngrid) && is.null(grideps))
    ngrid <- 10
    
  # Selected calculations
  opt.none <- locppmOptions(FALSE)
  if(!use.fft) {
    opt.coef <- locppmOptions(cg=TRUE)
    opt.t    <- locppmOptions(tg=TRUE)
    opt.hess <- locppmOptions(Vg=TRUE, Tg=TRUE)
    opt.hom  <- locppmOptions(sh=TRUE, fh=TRUE, gh=TRUE, Xh=!ispois)
    opt.lik  <- locppmOptions(sh=TRUE, fh=TRUE, gh=TRUE, Xh=!ispois, Lg=TRUE)
    opt.full <- locppmOptions(v0=FALSE, other=TRUE, other1=FALSE)
  } else {
    # FFT calculations
    opt.coef <- locppmOptions(cg1=TRUE)
    opt.t    <- locppmOptions(tg1=TRUE)
    opt.hess <- locppmOptions(gg1=TRUE) 
    opt.hom  <- locppmOptions(sh1=TRUE, fh1=TRUE, gh1=TRUE) #, Xh1=!ispois)
    opt.lik  <- opt.hom
    opt.full <- locppmOptions(FALSE, other1=TRUE)
  }
  
  switch(vcalc,
         none = {
           # estimate coefficients only
           opt.fit  <- opt.coef
           opt.var  <- opt.none
         },
         t = {
           # estimate coefficients and calculate t-statistics
           opt.fit <- opt.coef
           opt.var <- opt.t
         },
         hessian = {
           # estimate coefficients and Poisson/Poincare variance only
           opt.fit <- opt.coef
           opt.var <- opt.hess & !opt.coef
         },
         full = {
           # full variance estimation
           opt.fit <- opt.coef
           opt.var <- opt.full & !opt.coef
         },
         hom = {
           # For use by 'homtest'.
           # Don't evaluate fitted coefficients.
           # Calculate variance of local score under homogeneous model
           opt.fit <- opt.none
           opt.var <- opt.hom
         },
         lik = {
           # For use by 'homtest'.
           # Calculate local likelihood ratio test statistic,
           # plus variance of local score under homogeneous model
           opt.fit <- opt.none
           opt.var <- opt.lik
         }
         )

  phase1time <- NULL
  
  # Execute
  switch(locations,
         fine = {
           # quadrature points
           opt <- opt.fit | opt.var
           lpe <- locppmEngine(homfit, sigma, U,
                               weights=wU, opt=opt,
                               scopename="quadrature points",
                               verbose=verbose,
                               fft.algorithm=fft.algorithm)
         },
         coarse = {
           # grid points
           opt <- opt.fit | opt.var
           coarse.to.fine <- gridproxy(U, dimyx=ngrid, eps=grideps)
           G <- U[coarse.to.fine]
           wG <- attr(coarse.to.fine, "weights")
           lpe <- locppmEngine(homfit, sigma, G, weights=wG, opt=opt,
                               scopename="grid points",
                               verbose=verbose, fft.algorithm=fft.algorithm)
         },
         split = {
           # fit on quadrature points, variance estimation on grid points
           lpe.fit <- locppmEngine(homfit, sigma, U, weights=wU, opt=opt.fit,
                                   scopename="quadrature points",
                                   verbose=verbose, fft.algorithm=fft.algorithm)
           # record time taken to complete first phase 
           phase1time <- proc.time() - starttime
           #
           coarse.to.fine <- gridproxy(U, dimyx=ngrid, eps=grideps)
           G <- U[coarse.to.fine]
           wG <- attr(coarse.to.fine, "weights")
           if(verbose) cat("Estimating variance..\n")
           lpe.var <- locppmEngine(homfit, sigma, G, weights=wG, opt=opt.var,
                                   scopename="grid points",
                                   scopeindex=coarse.to.fine,
                                   verbose=verbose, fft.algorithm=fft.algorithm)
           lpe <- resolve.defaults(lpe.fit, lpe.var,
                                   list(coarse.to.fine=coarse.to.fine),
                                   .MatchNull=FALSE)
         })
  # pack up
  result <- append(lpe,
                   list(homfit=homfit, ispois=ispois, sigma=sigma,
                        vcalc=vcalc, locations=locations,
                        ngrid=ngrid, grideps=grideps,
                        templatecall=templatecall,
                        phase1time=phase1time))
  class(result) <- c("locppm", class(result))
  result <- timed(result, starttime=starttime)
  return(result)
}

.locppmOptionTable <- local({
  x <-
    list(cg="coefficient estimates",
         vg="variance of local fit",
         tg="t statistics of local fit",
         Vg="Poincare variance (inverse negative Hessian) of local fit",
         Tg="Poincare approximation of t statistics of local fit",
         xg="leave-one-out coefficient estimates",
         vh="null variance of local fit under homogeneous model",
         fh="local Fisher information under homogeneous model",
         gh="gradient of local score evaluated at homogeneous model",
         sh="local score of homogeneous model",
         Xh="covariance of local and global scores under homogeneous model",
         v0="variance of local fit under reduced model", 
         cg1="Taylor approximation of local coefficient estimates",
         vg1="variance of Taylor approximate local coefficients",
         tg1="t statistics of Taylor approximate local fit",
         gg1="gradient (negative Hessian) of Taylor-approximate local fit",
         vh1="null variance of local fit under homogeneous model (by FFT)",
         fh1="local Fisher information under homogeneous model",
         sh1="local score of homogeneous model",
         gh1="gradient of local score evaluated at homogeneous model",
         Xh1="covariance of local and global scores under homogeneous model",
         #
         Lg="local likelihood ratio test statistic for homogeneity")
  
  not.implemented <- "Xh1"
  require.fastRC <- c("sh", "fh", "Xh")

  y <- names(x)
  z <- unname(unlist(x))

  y1 <- substr(y, 1, 1)
  y2 <- substr(y, 2, 2)

  calcmap <- c(c="coef",
               v="var",
               f="fish",
               s="score",
               g="grad",
               t="tstat",
               V="invgrad",
               T="tgrad",
               L="lrts",
               x="xcoef",
               X="cov")

  dimtype <- c(c="vector",
               v="matrix",
               f="matrix",
               s="vector",
               g="matrix",
               t="vector",
               V="matrix",
               T="vector",
               L="scalar",
               x="vector",
               X="matrix")

  data.frame(tags        = y,
             descrip     = z,
             calctype    = factor(y1),
             calcname    = unname(calcmap[y1]),
             dimtype     = unname(dimtype[y1]),
             modeltype   = factor(ifelse(y2 %in% c("h", "0"), y2, "l")),
             usefft      = (nchar(y) == 3),
             implemented = !(y %in% not.implemented),
             requirefast = y %in% require.fastRC,
             stringsAsFactors=FALSE, row.names=y)
})

locppmOptions <- function(other=FALSE, ..., other1=other) {
  if(!is.logical(other) || !is.logical(other1)) stop("Logical values expected")
  # initialise options
  LOT <- .locppmOptionTable
  opt <- ifelse(LOT$usefft, other1, other)
  names(opt) <- LOT$tags
  # set options given explicitly
  argh <- list(...)
  nama <- names(argh)
  hit <- (nama %in% names(opt))
  if(any(hit)) {
    newvalues <- unlist(argh[hit])
    if(!is.logical(newvalues))  stop("Logical values expected")
    opt[nama[hit]] <- newvalues
  }
  # warn about unrecognised options
  if(any(nbg <- !hit))
    warning(paste("Unrecognised",
                  ngettext(sum(nbg), "argument", "arguments"),
                  commasep(sQuote(nama[!hit]))))
  # return
  class(opt) <- c("locppmOptions", class(opt))
  return(opt)
}

print.locppmOptions <- function(x, ...) {
  cat("Options for locppm fit\n")
  y <- unlist(x)
  if(!any(y)) {
    cat("Selected options: None\n")
  } else {
    cat("Selected options: \n")
    explain <- .locppmOptionTable$descrip[y]
    nama <- names(explain)
    for(i in seq(along=explain))
      cat(paste0("\t$", nama[i], ": ", explain[[i]], "\n"))
  }
  return(invisible(NULL))
}

locppmEngine <- function(model, sigma, V, 
                         ...,
                         weights=NULL,
                         verbose=TRUE,
                         scopename="points",
                         scopeindex = NULL,
                         opt = locppmOptions(cg=TRUE, Tg=TRUE),
                         dropterm = NULL,
                         matrices = FALSE,
                         fastRCinloop = TRUE,
                         internals = NULL,
                         fft.algorithm = "density"
                         ) {
  # compute weighted version of point process model at each point of V
  stopifnot(inherits(model, "ppm"))
  stopifnot(is.ppp(V))
  # ensure 'opt' contains all required entries and resolve defaults
  nullopt <- locppmOptions(cg=TRUE, Tg=TRUE, v0=!is.null(dropterm))
  opt <- resolve.defaults(as.list(opt), as.list(nullopt))
  # get information about *available* options
  LOT <- .locppmOptionTable
  tags        <- LOT$tags
  usefft      <- LOT$usefft
  dimtype     <- LOT$dimtype 
  iscoef      <- (LOT$calctype == "c")
  implemented <- LOT$implemented
  requirefast <- LOT$requirefast
  modeltype   <- LOT$modeltype
  calctype    <- LOT$calctype
  calcname    <- LOT$calcname
  # Detect options that are not yet implemented:
  if(any(nbg <- unlist(opt[tags[!implemented]]))) {
    nama <- names(nbg)[nbg]
    nbad <- sum(nbg)
    warning(paste(ngettext(nbad, "Option", "Options"),
                  commasep(sQuote(nama)),
                  ngettext(nbad, "is", "are"),
                  "not yet implemented"))
  }
  if(!fastRCinloop && any(needfast <- unlist(opt[tags[requirefast]]))) {
    tagz <- names(needfast)[needfast]
    desired <- paste(LOT$descrip[tagz], collapse=" or ")
    warning(paste("Slow code (fastRCinloop=FALSE) does not compute",
                  paste0(desired, ";"), "using fast code instead"))
    fastRCinloop <- TRUE
  }
  # Decide on type(s) of calculation
  # Should we fit the local GLM at each location?
  need.localfit <- any(unlist(opt[tags[modeltype == "l" & !usefft]]))
  # Do we need some kind of iteration over locations?
  need.iteration <- any(unlist(opt[tags[!usefft]]))
  # Do we need to update/refit/change the original model in any way?
  need.update <- any(unlist(opt[tags[modeltype != "h"]]))
  #
  # ... setup ........................
  coef.hom <- coef(model)
  nama <- names(coef.hom)
  ncoef <- length(coef.hom)
  ncoef2 <- ncoef^2
  #
  nV <- npoints(V)
  Vx <- V$x
  Vy <- V$y
  # initialise results
  for(tn in tags) assign(tn, NULL)
  # i.e.   cg <- vg <- tg <- ... <- NULL
  if(any(unlist(opt[tags[dimtype == "scalar"]]))) {
    # create template storage for scalars
    sc.blank <- numeric(nV)
  }
  if(need.localfit || any(unlist(opt[tags[dimtype == "vector"]]))) {
    # create template storage for coefficients and t-statistics
    ct.blank <- matrix(NA_real_, nrow=nV, ncol=ncoef)
    colnames(ct.blank) <- nama
  }
  if(any(unlist(opt[tags[dimtype == "matrix"]]))) {
    # create template storage for variance matrices
    v.blank <-  matrix(NA_real_, nrow=nV, ncol=ncoef2)
    colnames(v.blank) <- as.vector(outer(nama, nama, paste, sep="."))
  }
  # actually assign storage
  if(need.localfit) cg <- ct.blank
  for(i in seq_along(tags)) {
    if(opt[[i]]) {
      blanki <- switch(dimtype[i],
                       matrix = v.blank,
                       vector = ct.blank,
                       scalar = sc.blank)
      assign(tags[i], blanki)
    }
  }

  # .....................................
  # start computin'
  # .....................................

  if(is.null(internals)) {
    # precompute internal data for Rubak-Coeurjolly estimates?
    # Yes if we have to compute Variance, Fisher info or t-statistic
    need.internals <-
      tags[calctype %in% c("v", "f", "t", "X") & ((!usefft & fastRCinloop) |
                                             (usefft  & !is.poisson(model)))]
    if(any(unlist(opt[need.internals]))) 
      internals <- getvcinternals(model, verbose=verbose)
  }
  
  env.here <- sys.frame(sys.nframe())
  m <- getglmfit(model)
  if(is.null(m)) stop("model has no glm fit")
  env.model <- environment(terms(m))
  env.update <- environment(formula(m))
  if(need.update) {
    # manipulate environments so that the update will work
    fmla     <- get("fmla",     envir=env.model)
    gcontrol <- get("gcontrol", envir=env.model)
    glmdata  <- get("glmdata", envir=env.model)
    assign("coef.hom", coef.hom, envir=env.update)
    # to pacify package checker
    .mpl.W <- glmdata$.mpl.W
    .mpl.Y <- glmdata$.mpl.Y
    if(opt$v0) {
      if(is.null(dropterm))
        stop("Argument dropterm is missing")
      # construct GLM formula representing null model
      dropfmla <- paste(". ~ . - ", dropterm)
      assign("dropfmla", dropfmla, envir=env.update)
      # evaluate null model
      m0 <- update(m, dropfmla, evaluate=FALSE)
      m0 <- try(eval(m0, enclos=env.model, envir=env.here), silent=TRUE)
      if(inherits(m0, "try-error"))
        stop("Internal error: evaluation of null model failed")
      coef.hom0 <- coef(m0)
      assign("coef.hom0", coef.hom0, envir=env.update)
      # determine injection of coefficients of null into alternative
      nullmap <- match(names(coef.hom0), nama)
      if(any(is.na(nullmap)))
        stop("Internal error: cannot match coefficients of null to alternative")
      zerocoef <- rep(0, ncoef)
      names(zerocoef) <- nama
    }
  }

  if(any(unlist(opt[tags[usefft]]))) {
    # FFT calculations
    if(verbose) 
      cat("Performing FFT calculations...")
    # calculations based on homogeneous intensity
    ffthom <- tags[usefft & (modeltype == "h" | iscoef)]
    # calculations based on Taylor approximation to locally fitted intensity
    fftapp <- tags[usefft & modeltype == "l" & !iscoef]
    if(any(unlist(opt[fftapp])))
      opt$cg1 <- TRUE
    # 
    if(any(opt.hom <- unlist(opt[ffthom]))) {
      # calculations using homogeneous intensity
      tags.do <- names(opt.hom[opt.hom])
      HF <- locppmFFT(model, sigma=sigma, ..., 
                      what=calcname[tags %in% tags.do],
                      internals=internals,
                      algorithm=fft.algorithm,
                      verbose=verbose)
      # Extract values at locations V
      if(opt$cg1) cg1 <- sample.imagelist(HF$coefficients, V)
      if(opt$vh1) vh1 <- sample.imagelist(HF$variance,     V)
      if(opt$fh1) fh1 <- sample.imagelist(HF$fisher,       V)
      if(opt$sh1) sh1 <- sample.imagelist(HF$score,        V)
      if(opt$gh1) gh1 <- sample.imagelist(HF$gradient,     V)
    } 
    if(any(opt.app <- unlist(opt[fftapp]))) {
      # calculations using Taylor approximation to locally fitted intensity
      tags.do <- names(opt.app[opt.app])
      # extract Taylor approximations to coefficients at each quadrature point
      coTay <- sample.imagelist(HF$coefficients, union.quad(quad.ppm(model)))
      # compute approximate fitted intensity
      lamT <- locppmPredict(model, coTay)
      # discard intensity values, computed for different model
      forbid <- c("lambda", "lamdel")
      internals.clean <- internals[!(names(internals) %in% forbid)]
      # plug into FFT calculation
      TF <- locppmFFT(model, sigma=sigma,
                      lambda=lamT, ...,
                      what=calcname[tags %in% tags.do],
                      internals=internals.clean, 
                      algorithm=fft.algorithm, verbose=verbose)
      # Extract values at locations V
      if(opt$vg1) vg1 <- sample.imagelist(TF$variance,   V)
      if(opt$tg1) tg1 <- sample.imagelist(TF$tstatistic, V)
      if(opt$gg1) gg1 <- sample.imagelist(TF$gradient,   V)
    } 
    if(verbose) cat("Done.\n")
  }

  if(need.iteration) {
    # ............. Compute estimates by local fitting .................
    
    # full pattern of quadrature points 
    U <- union.quad(quad.ppm(model))
    nU <- npoints(U)
    Ux <- U$x
    Uy <- U$y

    if(opt$Lg) {
      ## precompute values for likelihood ratio test
      lambda.hom <- fitted(m)
      log.lambda.hom <- log(ifelse(lambda.hom > 0, lambda.hom, 1))
    }

    if(verbose)
      cat(paste("Processing", nV, scopename, "...\n"))

    for(j in 1:nV) {
      if(verbose)
        progressreport(j, nV)
      localwt <- dnorm(Ux - Vx[j], sd=sigma) * dnorm(Uy - Vy[j], sd=sigma)
      if(need.localfit) {
        assign("localwt", localwt, envir=env.update)
        fitj <- update(m, weights = .mpl.W * localwt, start=coef.hom,
                       evaluate=FALSE)
        fitj <- try(eval(fitj, enclos=env.model, envir=env.here), silent=TRUE)
        if(!inherits(fitj, "try-error") && fitj$converged) {
          # fitted coefficients
          cg[j, ] <- coefj <- coef(fitj)
          # Poincare t-statistic for each parameter (using Hessian)
          if(opt$Tg)
            Tg[j,] <- coef(summary(fitj))[,3]
          # local likelihood ratio test statistic for test of homogeneity
          if(opt$Lg) {
            lambda.j <- fitted(fitj)
            log.lambda.j <- log(ifelse(lambda.j > 0, lambda.j, 1))
            locW <- .mpl.W * localwt
            logLhom <- sum(locW * (.mpl.Y * log.lambda.hom - lambda.hom))
            logLj   <- sum(locW * (.mpl.Y * log.lambda.j   - lambda.j))
            Lg[j] <- 2 * (logLj - logLhom)
          }
          # variance of local parameter estimates under local model
          # (for confidence intervals)
          if(opt$Vg) { # Poincare variance of local fit
            # Note that for a GLM with weights,
            # vcov.glm returns the inverse of the negative Hessian
            Vgj <- try(vcov(fitj, dispersion=1), silent=TRUE)
            if(!inherits(Vgj, "try-error") && !is.null(Vgj)) 
              Vg[j,] <- as.vector(Vgj)
          }
          if(opt$vg || opt$tg) {
            # Rubak-Coeurjolly type estimate of variance of local fit
            if(fastRCinloop) {
              vgj <- vcovlocEngine(internals, localwt,
                                   A1dummy=TRUE, new.coef=coefj)
            } else {
              vgj <- there.is.no.try(vcov(model,
                                          matwt=localwt, new.coef=coefj,
                                          A1dummy=TRUE, matrix.action="silent"),
                                     silent=TRUE)
            }
            if(!is.null(vgj)) {
              if(opt$vg) vg[j,] <- as.vector(vgj)
              if(opt$tg) tg[j,] <- coefj/sqrt(diag(vgj))
            }
          }
        }
        if(opt$xg) {
          ## leave-one-out estimates of coefficients
          kk <- if(is.null(scopeindex)) j else scopeindex[j]
          ## V[j] = U[kk]
          if(glmdata[kk, ".mpl.Y"] > 0) { #' U[kk] is a data point
            fakedata <- glmdata
            fakedata[kk, ".mpl.Y"] <- 0 # pretend there's no data point at U[kk]
            assign("fakedata", fakedata, envir=env.model)
            fut.j <- update(m, data=fakedata,
                            weights = .mpl.W * localwt, start=coef.hom,
                            evaluate=FALSE)
            fut.j <- try(eval(fut.j, enclos=env.model, envir=env.here),
                         silent=TRUE)
            if(!inherits(fut.j, "try-error") && fut.j$converged) 
              xg[j, ] <- coef(fut.j)
          } else xg[j,] <- cg[j,]
        }
      }
      # variance of local parameter estimates under homogeneous model
      # (for tests of homogeneity)
      if(any(unlist(opt[c("vh", "fh", "gh", "sh", "Xh")]))) {
        if(fastRCinloop) {
          vhj <- vcovlocEngine(internals, localwt, A1dummy=TRUE,
                               bananas=opt$Xh)
          if(opt$gh) ghj <- attr(vhj, "Grad")
          if(opt$fh) fhj <- attr(vhj, "Fish")
          if(opt$sh) shj <- attr(vhj, "Score")
          if(opt$Xh) Xhj <- attr(vhj, "Cov")
        } else if(!opt$gh) {
          vhj <- try(vcov(model, matwt=localwt, A1dummy=TRUE,
                          matrix.action="silent"),
                     silent=TRUE)
          if(inherits(vhj, "try-error")) vhj <- NULL
          shj <- fhj <- Xhj <- NULL
        } else {
          tmp <- try(vcov(model, matwt=localwt, A1dummy=TRUE,
                          matrix.action="silent",
                          what="all"),
                     silent=TRUE)
          if(inherits(tmp, "try-error")) vhj <- ghj <- shj <- NULL else {
            vhj <- tmp$varcov
            ghj <- tmp$internals$gradient
          }
          shj <- fhj <- Xhj <- NULL
        }
        if(opt$vh && !is.null(vhj) &&
           length(as.vector(vhj)) == ncoef2)
          vh[j,] <- as.vector(vhj)
        if(opt$gh && !is.null(ghj) &&
           length(as.vector(ghj)) == ncoef2)
          gh[j,] <- as.vector(ghj)
        if(opt$fh && !is.null(fhj) &&
           length(as.vector(fhj)) == ncoef2)
          fh[j,] <- as.vector(fhj)
        if(opt$sh && !is.null(shj) &&
           length(as.vector(shj)) == ncoef)
          sh[j,] <- as.vector(shj)
        if(opt$Xh && !is.null(Xhj) &&
           length(as.vector(Xhj)) == ncoef2)
          Xh[j,] <- as.vector(Xhj)
      }
      # variance of local parameter estimates under local NULL model
      # (for local test of H_0: beta_k = 0)
      if(opt$v0) {
          # fit null model
        fit0j <- update(m, dropfmla,
                        weights = .mpl.W * localwt, start=coef.hom0,
                        evaluate=FALSE)
        fit0j <- try(eval(fit0j, enclos=env.model, envir=env.here), silent=TRUE)
        if(!inherits(fit0j, "try-error")) {
          # extract coefficients of null model
          c0j <- coef(fit0j)
          # inject into alternative
          coef0j <- zerocoef
          coef0j[nullmap] <- c0j
          # compute null variance of local parameter estimates
          if(fastRCinloop) {
            v0j <- vcovlocEngine(internals, localwt,
                                 A1dummy=TRUE, new.coef=coef0j)
          } else {
            v0j <- try(vcov(model, matwt=localwt, new.coef=coef0j,
                            A1dummy=TRUE, matrix.action="silent"),
                       silent=TRUE)
          }
          if(!inherits(v0j, "try-error") && !is.null(v0j) &&
             length(as.vector(v0j)) == ncoef2)
            v0[j,] <- as.vector(v0j)
        }
      }
    }
    # ........ end of loop over V .................................
    if(verbose) cat("Done.\n")
  }

  # convert format and/or add attributes
  make.ssf <- !matrices
  add.scope <- !is.null(scopeindex)
  add.scopename <- !is.null(scopename)
  if(make.ssf || add.scope) {
    for(tn in tags) {
      if(!is.null(thing <- get(tn))) {
        if(make.ssf) {
          thing <- ssf(V, thing)
          if(!is.null(weights)) attr(thing, "weights") <- weights
        }
        if(add.scope)
          attr(thing, "scopeindex") <- scopeindex
        if(add.scopename)
          attr(thing, "scopename") <- scopename
        assign(tn, thing)
      }
    }
  }
        
  result <- mget(tags)
  return(result)
}

# Locally-weighted Rubak-Coeurjolly estimate of variance

vcovlocEngine <- function(internals, localwt=NULL,
                          ...,
                          A1dummy=FALSE, new.coef = NULL,
                          bananas=FALSE) {
  # 'internals' has previously been validated.
  # Components include
  #     'mom'       model matrix 
  #     'lambda'    fitted intensity of the homogeneous model
  #     'Z'         data/dummy indicator
  #     'wQ'        quadrature weight
  #     'ok'        indicator of the domain of the pseudolikelihood
  #     'areaW'     area of the domain of the pseudolikelihood
  #     'nX'        number of data points
  #     'ispois'    TRUE if model is Poisson
  #     'hom.coef'  fitted coefficients of homogeneous model
  with(internals, {
    okX <- ok[Z]
    ncoef <- length(hom.coef)
    # Conditional intensity using the locally-fitted coefficients 'new.coef'
    if(!is.null(new.coef)) 
      lambda <- as.vector(lambda * exp(mom %*% (new.coef - hom.coef)))
    # Matrix A1
    momL <- mom * localwt
    A1 <- with(internals,
               if(A1dummy) sumouter(momL[ok, , drop=FALSE],
                                    w=(lambda * wQ)[ok])
               else sumouter(momL[Z & ok, , drop=FALSE]))
    # Gradient (sensitivity) matrix
    Grad <- with(internals,
               if(A1dummy) sumouter(mom[ok, , drop=FALSE],
                                    w=(localwt* lambda * wQ)[ok])
               else sumouter(mom[Z & ok, , drop=FALSE],
                             w = localwt[Z & ok]))
    #
    # Matrices A2 and A3
    if(ispois) {
      A2 <- A3 <- B2 <- B3 <- matrix(0, ncoef, ncoef)
    } else {
      # Require components 'lamdel' and 'momdel'
      #   lamdel[i,j]   = lambda(X[i] | X[-j]) = lambda(X[i] | X[-c(i,j)])
      #   momdel[ ,i,j] = h(X[i] | X[-j])      = h(X[i] | X[-c(i,j)])
      # adjust lamdel for new coefficient
      if(!is.null(new.coef))
        lamdel <- lamdel * exp(tensor::tensor(new.coef - hom.coef, momdel, 1, 1))
      #   pairweight[i,j] = lamdel[i,j]/lambda[i] - 1 
      pairweight <- lamdel / lambda[Z] - 1
      #   momdelL[ , i, j] = localwt[i] * momdel[ , i, j]
      locwtX <- localwt[Z]
      momdelL <- momdel * rep(locwtX, rep(ncoef, nX))
      # now compute sum_{i,j} for i != j
      # pairweight[i,j] * outer(momdelL[,i,j], momdelL[,j,i])
      # for data points that contributed to the pseudolikelihood
      pwXok <- pairweight[okX, okX]
      locwtXok <- locwtX[okX]
      A2 <- sumsymouter(momdelL[, okX, okX], w=pwXok)
      if(bananas)
        B2 <- sumsymouter( momdel[, okX, okX], w=locwtXok * pwXok)
      # locally-weighted model matrix for data only
      momLX <- momL[Z, , drop=FALSE]
      # deltamomL[ ,i,j] = momLX[j,] - momdelL[,j,i]
      deltamomL <- aperm(rep(t(momLX), nX) - momdelL, c(1, 3, 2))
      A3 <- sumsymouter(deltamomL[, okX, okX])
      if(bananas) {
        momX  <- mom[Z, , drop=FALSE]
        deltamom <- aperm(rep(t(momX), nX) - momdel, c(1, 3, 2))
        nXok <- sum(okX)
        locwtXokMat <- matrix(locwtXok, nXok, nXok)
        B3 <- sumsymouter(deltamom[, okX, okX], w=locwtXokMat)
      }
    }
    ## Finally calculate the Rubak-Coeurjolly estimate:
    Sigma <- A1 + A2 + A3
    U <- try(solve(Grad/areaW))
    if(inherits(U, "try-error")) return(NULL)
    mat <- U %*% (Sigma/areaW) %*% U / areaW
    # add information
    attr(mat, "Grad") <- Grad
    attr(mat, "Fish") <- Sigma
    attr(mat, "Score") <- colSums((Z - lambda * wQ) * momL)
    if(bananas) attr(mat, "Cov") <- Grad + B2 + B3
    return(mat)
  })
}

getvcinternals <- function(model, verbose=TRUE) {
  # Compute the internal data needed for Rubak-Coeurjolly variance estimate
  if(verbose) cat("Computing internal variance data ...")
  ispois <- is.poisson(model)
  internals <- vcov(model, what="internals", saveterms=TRUE)
  needed <- c( c("lambda", "mom", "Z", "ok"),
              if(ispois) NULL else c("lamdel", "momdel"))
  hit <- needed %in% names(internals)
  if(!all(hit))
    stop(paste(ngettext(sum(!hit), "component", "components"),
               commasep(sQuote(needed[!hit])),
               "missing from internals"))
  # add more info
  wQ <- w.quad(quad.ppm(model))
  nX <- npoints(data.ppm(model))
  W <- as.owin(model)
  internals <- append(internals,
                      list(wQ=wQ, nX=nX, ispois=ispois, hom.coef=coef(model)))
  if(is.null(internals$areaW)) {
    internals$areaW <- if(model$correction == "border")
      eroded.areas(W, model$rbord) else area.owin(W)
  }
  if(verbose) cat("Done.\n")
  return(internals)
}

# ............. methods for locppm ................................

plot.locppm <- function(x, ...,
                        what="cg",
                        which=NULL) {
  xname <- deparse(substitute(x))
  if(!missing(what)) {
    what <- match.arg(what, row.names(.locppmOptionTable))
  } else {
    what <- FirstExtantEntry(x, c("cg", "cg1"), "Please specify argument what")
  }
  Z <- x[[what]]
  if(is.null(Z))
    stop(paste("Component", sQuote(what), "is unavailable"))
  if(!is.null(which))
    Z <- Z[, which]
  do.call("plot", resolve.defaults(list(x=Z),
                                   list(...),
                                   list(main=xname)))
}

print.locppm <- function(x, ...) {
  homfit <- as.ppm(x)
  splat("Local",
        if(is.poisson(homfit)) "likelihood" else "pseudolikelihood",
        "fit")
  splat("Template model:")
  splat(x$templatecall, "\n\n", indent=5)
  #' quadrature info
  qux <- quad.ppm(homfit)
  rez <- summary(qux)$resolution
  if(!is.null(rez)) {
    ndum <- npoints(qux$dummy)
    splat(ndum, "dummy points with spacing", format(rez))
  } else print(qux)
  #'
  splat("\nSmoothing parameter sigma: ",
        format(signif(x$sigma, 4) %unit% unitname(qux)),
        "\n")
  # Identify which values have been computed
  LOT <- .locppmOptionTable
  possible <- row.names(LOT)
  present <- possible %in% names(x)
  present[present] <- !unlist(lapply(x[possible[present]], is.null))
  if(any(present)) {
    tags    <- LOT$tags[present]
    explain <- LOT$descrip[present]
    scopes  <- unlist(lapply(lapply(x[tags], attributes),
                             getElement, name="scopename"))
    if(any(nbg <- !nzchar(scopes)))
      scopes[nbg] <- "unspecified locations"
    for(sc in unique(scopes)) {
      splat("Values computed at", sc)
      for(i in which(scopes == sc)) 
        splat(paste0(tags[i], ": ", explain[[i]]), indent=5)
    }
  }
  if(!is.null(ng <- x$ngrid))
    splat("Coarse grid:", paste(ensure2vector(ng), collapse=" x "))
  if(!is.null(ep <- x$grideps)) {
    uh <- unitname(x$homfit)
    splat("Coarse grid: spacing", format(ep %unit% uh))
  }
  if(!is.null(p1t <- x$phase1time))
    splat("Time taken for first phase:", codetime(p1t))
  return(invisible(NULL))
}

contour.locppm <- function(x, ...,
                           what="cg",
                           which=NULL) {
  xname <- deparse(substitute(x))
  if(!missing(what)) {
    what <- match.arg(what, row.names(.locppmOptionTable))
  } else {
    what <- FirstExtantEntry(x, c("cg", "cg1"), "Please specify argument what")
  }
  Z <- x[[what]]
  if(is.null(Z))
    stop(paste("Component", sQuote(what), "is unavailable"))
  W <- rescue.rectangle(as.owin(Z))
  if(!is.null(which))
    Z <- Z[, which]
  Z <- Smooth(Z, ...)
  do.call("contour", resolve.defaults(list(x=Z),
                                      list(...),
                                      list(main=xname)))
  invisible(NULL)
}

Smooth.locppm <- function(X, ..., what="cg") {
  stopifnot(inherits(X, "locppm"))
  if(!missing(what)) {
    what <- match.arg(what, row.names(.locppmOptionTable))
  } else {
    what <- FirstExtantEntry(X, c("cg", "cg1"), "Please specify argument what")
  }
  Y <- X[[what]]
  if(is.null(Y)) return(NULL)
  A <- as.ppp(Y)
  ok <- attr(Y, "ok")
  A <- A[ok]
  sigma0 <- if(any(c("sigma", "varcov") %in% names(list(...))))
            NULL else 1.4 * max(nndist(A))
  out <- do.call("Smooth.ppp",
                 resolve.defaults(list(X = A),
                                  list(...),
                                  list(sigma=sigma0)))
  return(out)
}

coef.locppm <- function(object, ...,
                        which=c("local", "homogeneous")) {
  which <- match.arg(which)
  result <- with(object,
                 switch(which,
                        homogeneous = coef(homfit),
                        local       = {
                          if(!is.null(cg)) cg else cg1
                        }))
  return(result)
}

confint.locppm <- function (object, parm, level = 0.95, ...,
                            which=c("local", "homogeneous")) 
{
  stopifnot(inherits(object, "locppm"))
  which <- match.arg(which)
  ispois <- is.poisson(object)
  cf <- coef(object$homfit)
  pnames <- names(cf)
  if (missing(parm)) 
    parm <- pnames
  else if (is.numeric(parm)) 
    parm <- pnames[parm]
  if(which == "homogeneous")
    return(confint(object$homfit, parm, level, ...))
  nparm <- length(parm)
  a <- (1 - level)/2
  a <- c(a, 1 - a)
  pct <- paste(format(100 * a, trim = TRUE,
                      scientific = FALSE, digits = 3), 
               "%")
  fac <- qnorm(a)
  # local var-cov matrices
  v <- object$vg
  if(is.null(v)) {
    v <- object$Vg
    if(is.null(v))
      stop("Fitted model does not contain variances")
    else
      warning("Using Poincare variance (Hessian matrix)")
  }
  vindex <- attr(v, "scopeindex")
  v <- marks(v)
  # extract standard deviations
  diagindex <- diag(matrix(1:(nparm^2), nparm, nparm))
  sd <- sqrt(v[, diagindex, drop=FALSE])
  # create space for result
  ci <- matrix(NA_real_, nrow=nrow(v), ncol = 2 * nparm)
  colnames(ci) <- as.vector(t(outer(parm, pct, paste)))
  # extract coefficients at same locations
  coefs <- object$cg
  cindex <- attr(coefs, "scopeindex")
  if(!identical(cindex, vindex)) {
    # v is computed on a coarser set of points
    coarse.to.fine <- object$coarse.to.fine
    coefs <- coefs[coarse.to.fine, ]
  }
  co <- marks(coefs)
  P  <- unmark(coefs)
  if(nrow(co) != nrow(v))
    stop("Internal error: mismatch in arrays")
  # calculate confidence intervals
  for(i in seq_len(nrow(v))) 
    ci[i,] <- rep(co[i,], rep(2, nparm)) +
              rep(sd[i, ],   rep(2, nparm)) * rep(fac, nparm)
  # pack up
  ssf(P, ci)
}

as.ppm.locppm <- function(object) {
  return(object$homfit)
}

as.interact.locppm <- function(object) {
  as.interact(as.ppm(object))
}

fitted.locppm <- function(object, ...,
                          type=c("cif", "trend", "intensity"),
                          new.coef=NULL) {
  trap.extra.arguments(...)
  type <- match.arg(type)
  lam <- predict(object, type=type, new.coef=new.coef)
  return(marks(lam))
}

predict.locppm <- function(object, ...,
                           type=c("cif", "trend", "intensity"),
                           locations=NULL, new.coef=NULL) {
  # minimal implementation
  trap.extra.arguments(...)
  type <- match.arg(type)
  locations.given <- !is.null(locations)
  # Extract homogeneous/template model
  homfit <- as.ppm(object)
  # Fitted local coefficients (could be NULL)
  coefs <- coef(object, what="local")
  # Prediction locations
  if(!locations.given) {
    # Use sample points where fitted local coefficients were computed
    if(is.null(coefs)) stop("Unable to determine locations for prediction")
    locations <- unmark(coefs)
  } else stopifnot(is.ppp(locations))
  # Local coefficients
  index <- NULL
  if(is.null(new.coef)) {
     # Extract local coefficients 
    if(is.null(coefs)) stop("Object does not include any fitted coefficients")
    coefmat <- marks(coefs)
    if(locations.given) {
      #' extract coefficient at nearest sample point
      map <- nncross(locations, as.ppp(coefs), what="which")
      coefmat <- coefmat[map, , drop=FALSE]
    }
  } else {
    # New values for local coefficients
    p <- length(coef(as.ppm(object)))
    nloc <- npoints(locations)
    if(is.matrix(new.coef)) {
      if(ncol(new.coef) != p)
        stop("Incorrect number of columns in new.coef")
      if(nrow(new.coef) != nloc)
        stop("Incorrect number of rows in new.coef")
      coefmat <- new.coef
    } else {
      if(length(new.coef) == p) {
        # replicate
        coefmat <- matrix(new.coef, nloc, p, byrow=TRUE)
      } else stop("Incorrect length in new.coef")
    }
  }
  # Check that local coefficients were computed at the quadrature points
  if(is.null(new.coef) &&
     identical(attr(coefs, "scopename"), "quadrature points")) {
    precomputed <- NULL # use defaults
  } else {
    mom <- sample.imagelist(model.images(homfit), locations)
    precomputed <- list(hom.coef = coef(homfit),
                        mom = mom)
    switch(type,
           cif = {
             precomputed$lambda <- predict(homfit, locations=locations)
           },
           trend = ,
           intensity = {
             precomputed$trend <- predict(homfit, locations=locations,
                                          type="trend")
           })
  }
  # compute
  values <- locppmPredict(homfit, coefmat, type=type,
                          locations=locations, precomputed=precomputed)
  if(locations.given) return(values)
  return(ssf(locations, values))
}

locppmPredict <- function(homfit, coefs,
                          type = c("cif", "trend", "intensity"),
                          locations = NULL,
                          precomputed=NULL, details=FALSE,
                          index = NULL) {
  stopifnot(is.ppm(homfit))
  stopifnot(is.matrix(coefs))
  type <- match.arg(type)
  # Start by computing fitted conditional intensity of homogeneous model
  # This ensures correct handling of offsets etc
  coef0   <- precomputed$hom.coef %orifnull%  coef(homfit)
  mom     <- precomputed$mom %orifnull% 
             (if(is.null(locations)) model.matrix(homfit) else
              sample.imagelist(model.images(homfit), locations))
  switch(type,
         cif = {
           result0 <- precomputed$lambda %orifnull%
                      (if(is.null(locations)) fitted(homfit) else
                       predict(homfit, locations=locations))
         },
         trend = ,
         intensity = {
           result0 <- precomputed$trend %orifnull%
                      (if(is.null(locations)) fitted(homfit, type="trend") else
                       predict(homfit, type="trend", locations=locations))
         })
  if(!is.null(index)) {
    ## restrict to subset of quadrature points U[index]
    result0 <- result0[index] 
    mom     <- mom[index,,drop=FALSE]
  }
  # compute change in linear predictor
  d.coef <- coefs - matrix(coef0, nrow(coefs), ncol(coefs), byrow=TRUE)
  if(type != "cif") {
    # ignore interaction terms
    Vnames <- homfit$internal$Vnames
    d.coef[, Vnames] <- 0
  }
  # adjust fitted conditional intensities according to changed coefficients
  d.eta <- rowSums(mom * d.coef)
  result <- result0 * exp(d.eta)
  #
  if(type == "intensity")
    stop("Poisson-saddlepoint approximation is not yet implemented")
  #
  ans <- result
  if(details)
    attr(ans, "d.eta") <- d.eta
  return(ans)
}


is.poisson.locppm <- function(x) { is.poisson(x$homfit) }


# local 't' statistic for one coefficient in model
ttestmap <- function(object, term, ...,
                     method = c("exact", "hessian", "taylor"),
                     grid = FALSE,
                     ngrid=NULL, grideps=NULL,
                     verbose=TRUE) {
  starttime <- proc.time()
  method <- match.arg(method)
  stopifnot(inherits(object, "locppm"))
  homfit <- object$homfit
  coef.hom <- coef(homfit)
  nama <- names(coef.hom)
  # validate 'term'
  stopifnot(is.character(term))
  gfit <- getglmfit(homfit)
  tlab <- attr(terms(gfit), "term.labels")
  tpos <- match(term, tlab)
  if(is.na(tpos))
    stop(paste(sQuote(term), "is not a term in the model formula"))
  ass <- attr(model.matrix(gfit), "assign")
  relevant <- which(ass == tpos)
  if(length(relevant) == 0)
    stop("Internal error: cannot match the term to its canonical coefficients")
  if(length(relevant) > 1 && method == "exact")
    stop(paste("For Gibbs models, the exact method is not yet implemented",
               "for multidimensional parameters;",
               "the term", sQuote(term),
               "corresponds to", length(relevant), "parameters",
               commasep(sQuote(nama[relevant]))))
  parm <- nama[relevant]
  iparm <- relevant
  #
  # Which t statistic?
  tname <- switch(method,
                  exact   = "tg",
                  hessian = "Tg",
                  taylor  = "tg1")
  # Is it already available?
  if(!grid && !is.null(tvalues <- object[[tname]])) {
    result <- tvalues[,parm]
    result <- timed(result, starttime=starttime)
    return(result)
  }
  # Further computation required
  coefs <- coef(object)
  sigma <- object$sigma
  P <- unmark(coefs)
  wP <- attr(coefs, "weights")
  # determine points for evaluation
  if(!grid) {
    # use all quadrature points
    Puse <- P
    wPuse <- wP
    scopename <- "quadrature points"
    coarse.to.fine <- seq_len(npoints(P))
  } else {
    # use an approximate grid
    if(is.null(ngrid) && is.null(grideps) && object$locations == "split") {
      # use existing 'grid'
      coarse.to.fine <- object$coarse.to.fine
    } else {
      # generate new 'grid'
      coarse.to.fine <- gridproxy(P, dimyx=ngrid, eps=grideps, weights=wP)
    }
    Puse <- P[coarse.to.fine]
    wPuse <- attr(coarse.to.fine, "weights")
    scopename <- "grid points"
  }
  # Compute
  opt <- switch(method,
                exact   = locppmOptions(tg  = TRUE),
                taylor  = locppmOptions(tg1 = TRUE),
                hessian = locppmOptions(Tg = TRUE))
  z <- locppmEngine(homfit, sigma, Puse, weights=wPuse, 
                    opt=opt, scopename=scopename,
                    verbose=verbose)
  tvalues <- z[[tname]]
  result <- tvalues[,parm]
  result <- timed(result, starttime=starttime)
  return(result)
}

# (signed square root of) local score test statistic
# for test of homogeneity.

homteststat <- function(object, ..., verbose=FALSE) {
  # 'object' can be either a locppm or a homtestmap
  if(inherits(object, "locppm")) {
    Tv <- homtestmap(object, ..., verbose=verbose)
  } else if(inherits(object, "homtestmap")) {
    Tv <- update(object, ...)
  } else stop("Unrecognised format for object")
  starttime <- proc.time()
  Sv <- if(ncol(marks(Tv)) > 1) sqmag(Tv) else Tv
  S <- integral(Sv)/area(Window(Sv))
  timetaken <- proc.time() - starttime + attr(Tv, "timetaken")
  S <- timed(S, timetaken=timetaken)
  attr(S, "nsample") <- npoints(Sv)
  return(S)
}

update.homtestmap <-
  function(object, ...,
           what=NULL, test=NULL, ladjust=NULL,
           calibrate=NULL, saveall=FALSE, poolmoments=NULL) {
  starttime <- proc.time()
  trap.extra.arguments(...)
  stopifnot(inherits(object, "homtestmap"))
  ## arguments default to their current values in the object
  info <- attr(object, "info")
  oldwhat <- info$what %orifnull% "components"
  oldtest <- info$test %orifnull% "score"
  oldcalibrate <- info$calibrate %orifnull% "Satterthwaite"
  oldladjust   <- info$ladjust   %orifnull% "none"
  oldpoolmoments <- !identical(info$poolable, FALSE)
  if(is.null(what)) what <- oldwhat 
  if(is.null(test)) test <- oldtest 
  if(is.null(calibrate)) calibrate <- oldcalibrate
  if(is.null(ladjust)) ladjust <- oldladjust 
  if(is.null(poolmoments)) poolmoments <- oldpoolmoments
  what <- match.arg(what, c("components", "statistic", "pvalue"))
  test <- match.arg(test, c("score", "taylor", "likelihood"))
  calibrate <- match.arg(calibrate, c("chisq", "Satterthwaite", "firstmoment"))
  ladjust <- match.arg(ladjust, c("none", "moment", "PSS"))
  poolmoments <- as.logical(poolmoments)
  #' interpret arguments
  do.test <- (what %in% c("statistic", "pvalue"))
  likely <- test %in% c("likelihood", "taylor")
  adjusting <- (do.test && likely && ladjust != "none")
  newinfo <- replace(info,
                     c("what", "test", "calibrate", "ladjust",
                       "poolmoments", "adjusting"),
                     list(what, test, calibrate, ladjust,
                          poolmoments, adjusting))
  # is complete data available?
  if(!is.null(localdata <- attr(object, "localdata"))) {
    # yes - compute whatever is required
    resultvalues <- HomTestMapEngine(localdata, newinfo)
    P <- unmark(object)
    result <- ssf(P, resultvalues)
    class(result) <- c("homtestmap", class(result))
    attr(result, "info") <- newinfo
    if(saveall) attr(result, "localdata") <- localdata
    result <- timed(result, starttime=starttime)
    return(result)
  }
  # check compatibility
  sameladjust    <- (ladjust == oldladjust) 
  samecalibrate <- (calibrate == oldcalibrate) &&
                   (poolmoments == oldpoolmoments || calibrate == "chisq")
  sametestname <- (test == oldtest)
  sametestdetails <- switch(test,
                            score = TRUE,
                            taylor = samecalibrate,
                            likelihood = sameladjust && samecalibrate)
  sametest <- sametestname && sametestdetails

  ## 
  if(!sametest)
    stop("Cannot convert between different tests: no saved data")

  ##
  if(what == oldwhat) return(object)
  
  if(oldwhat == "components" && what == "statistic" && test != "likelihood") {
    result <- sqmag(object)
    class(result) <- c("homtestmap", class(result))
    attr(result, "info") <- newinfo
    result <- timed(result, starttime=starttime)
    return(result)
  }

  if(what == "pvalue" && test == "score" && calibrate == "chisq") {
    switch(oldwhat,
           pvalue = return(object),
           statistic = {
             stat <- marks(object)
           },
           components = {
             stat <- marks(sqmag(object))
           })
    ncoef <- info$p 
    pvals <- pchisq(stat, df=ncoef, lower.tail=FALSE)
    result <- ssf(unmark(object), pvals)
    result <- timed(result, starttime=starttime)
    return(result)
  }
  
  stop("Cannot recover information: no saved data")
}

homtestmap <-  function(object, ...,
                        what = c("components", "statistic", "pvalue"),
                        test = c("score", "taylor", "likelihood"),
                        ladjust = c("none", "moment", "PSS"),
                        calibrate = c("chisq", "Satterthwaite", "firstmoment"),
                        simple = !is.null(theta0), 
                        theta0 = NULL, 
                        poolmoments = NULL,
                        sigma = NULL, 
                        saveall = FALSE, 
                        use.fft = TRUE, 
                        verbose = TRUE) {
  starttime <- proc.time()
  test <- match.arg(test)
  what <- match.arg(what)
  ladjust <- match.arg(ladjust)
  told.use.fft <- !missing(use.fft) && use.fft
  stopifnot(inherits(object, "locppm"))
  if(is.null(poolmoments)) poolmoments <- is.stationary(as.ppm(object))
  trap.extra.arguments(...)
  #' interpret arguments
  do.test <- (what %in% c("statistic", "pvalue"))
  likely <- test %in% c("likelihood", "taylor")
  adjusting <- do.test && likely && (ladjust != "none")
  #' calibration i.e. reference distribution for p-value
  miss.cal <- missing(calibrate)
  should.cal.chisq <- adjusting || (do.test && (test == "score"))
  should.cal.other <- do.test && (test != "score") && (ladjust == "none")
  if(miss.cal && should.cal.chisq) {
    #' score test statistic &
    #' adjusted composite likelihood ratio test statistic
    #' should be referred to chi-squared 
    calibrate <- "chisq"
  } else if(miss.cal && should.cal.other) {
    calibrate <- "Satterthwaite"
  } else {
    calibrate <- match.arg(calibrate)
    if(calibrate != "chisq" && should.cal.chisq) {
      statname <- if(test == "score") "Score test statistic" else
                  "Adjusted composite likelihood ratio test statistic" 
      warning(paste(statname, "should be referred to",
                    "the chi-squared distribution",
                    "by setting calibrate='chisq'"),
              call.=FALSE)
    }
    if(calibrate == "chisq" && should.cal.other) {
      warning(paste("Unadjusted likelihood-type statistic should be calibrated",
                    "by setting calibrate='firstmoment' or 'Satterthwaite'"),
              call.=FALSE)
    }
  }
  ispois <- is.poisson(object)
  if(test == "likelihood" && use.fft) {
    if(told.use.fft)
      stop("Cannot use FFT code for test='likelihood'")
    use.fft <- FALSE
  }
  # Determine whether variance calculation under H0 is required
  calc.pvals <- (what == "pvalue")
  calc.var <- calc.pvals || (test != "taylor")
  # Determine whether to assume a simple or composite null hypothesis
  nulltype <- if(simple) "simple" else "composite"
  if(calc.var && nulltype == "simple" && verbose)
    message("Note: calculation ignores the effect of estimating theta")
  if(!is.null(theta0) && use.fft) {
    use.fft <- FALSE
    if(told.use.fft) 
      stop("Cannot use FFT code when theta0 is given")
  }
  homfit <- object$homfit
  hom.coef <- coef(homfit)
  p <- ncoef <- length(hom.coef)
  if(is.null(sigma))
    sigma <- object$sigma

  vec.names <- names(hom.coef)
  mat.names <- as.vector(outer(vec.names, vec.names,
                               function(a,b){paste(a,b,sep=".")}))
  
  # determine what values are required
  # U = score, H = hessian, F = local Fisher info, L = local LRTS
  # X = cross-covariance 
  # taylor approximation to likelihood ratio requires H
  # adjusted composite likelihood ratio requires F, H
  # variance calculation for simple null requires F
  # variance calculation for composite null requires F, H and (if Gibbs) X
  # p-value for exact likelihood ratio requires F, H and (if Gibbs) X
  needed <- c(U = TRUE,
              F = saveall || calc.var,
              H = saveall || calc.pvals ||
                  (calc.var && nulltype == "composite") ||
                  (test=="taylor") || adjusting,
              L = saveall || (test == "likelihood"),
              X = !ispois && (calc.pvals ||
                              (calc.var && (nulltype == "composite"))))

  # try to use FFT results
  if(use.fft) {
    # check that required data are available
    ftable = c(U="sh1", F="fh1", H="gh1", L="NA", X="Xh1")
    missed <- unlist(lapply(object[ftable[needed]], is.null))
    if(any(missed)) {
      use.fft <- FALSE
      if(told.use.fft) {
        want.tags <- ftable[needed[missed]]
        want.desc <- with(.locppmOptionTable, descrip[tags %in% want.tags])
        gripe <- paste("Cannot use FFT code; object does not contain",
                       paste0(commasep(want.desc), ";"),
                       "refit the model with vcalc='hom' or 'lik'")
        warning(gripe, call.=FALSE)
      }
    }
  }
  # otherwise try to use data in object
  if(!use.fft) {
    ftable <- c(U="sh", F="fh", H="gh", L="Lg", X="Xh")
    missed <- unlist(lapply(object[ftable[needed]], is.null))
    names(missed) <- names(needed)[needed]
    if(needed["L"] && missed["L"])
      stop(paste("Cannot perform test='likelihood':",
                 "object does not contain",
                 "local likelihood ratio test statistic;",
                 "refit the model with vcalc='lik'"))
    use.object <- !any(missed) && is.null(theta0)
  }

  # OK, all parameters are decided
  resultinfo <- list(test=test, what=what, use.fft=use.fft, p=p,
                     sigma=sigma,
                     poolmoments=poolmoments,
                     poolable=is.stationary(as.ppm(object)),
                     calibrate=calibrate,
                     ladjust=ladjust,
                     nulltype = nulltype,
                     adjusting = adjusting,
                     saveall=saveall)

  # ................ START COMPUTING .......................................
  # .........  Compute the local score and other required matrices .........
  Fmat <- Hmat <- Umat <- VUmat <- Lvec <- Xmat <- NULL
  if(use.fft) {
    # use FFT results
    sh1 <- object$sh1
    Umat <- marks(sh1)
    P <- unmark(sh1)
    nP <- npoints(P)
    wP <- attr(sh1, "weights")
    if(needed["F"]) Fmat <- marks(object$fh1)
    if(needed["H"]) Hmat <- marks(object$gh1)
    if(needed["X"]) Xmat <- marks(object$Xh1)
  } else if(use.object) {
    # use values in object
    sh <- object$sh
    Umat <- marks(sh)
    P <- unmark(sh)
    nP <- npoints(P)
    wP <- attr(sh, "weights")
    if(needed["F"]) Fmat <- marks(object$fh)
    if(needed["H"]) Hmat <- marks(object$gh)
    if(needed["L"]) Lvec <- marks(object$Lg)
    if(needed["X"]) Xmat <- marks(object$Xh)
  } else {
    # compute directly (slower...)
    internals <- getvcinternals(homfit, verbose=verbose)
    # extract quadrature points used to fit model
    Q    <- quad.ppm(homfit, drop=FALSE)
    w    <- w.quad(Q)
    Z    <- is.data(Q)
    U    <- union.quad(Q)
    nU   <- npoints(U)
    wU <- w.quad(Q)
    Ux <- U$x
    Uy <- U$y
    ok <- getglmsubset(homfit)
#    suf <- model.matrix(homfit)[ok, , drop=FALSE]
    suf <- model.matrix(homfit)
    # fitted intensity of homogeneous model
    lam <- fitted(homfit, drop=FALSE, new.coef=theta0)
    # increments of residual measure of homogeneous fit
    homresid <- Z - lam * w
    #
    if(ispois) {
      # use all quadrature points
      P <- U
      wP <- wU
      scopename <- "quadrature points"
      coarse.to.fine <- seq_len(nU)
    } else {
      # use coarse grid
      # find components inside the object which are 'ssf'
      izsf <- unlist(lapply(object, inherits, what="ssf"))
      sfs <- object[izsf]
      # count number of points in each ssf object
      nps <- unlist(lapply(sfs, npoints))
      # select smallest number of points
      ibest <- which.min(nps)
      bestname <- names(sfs)[ibest]
      best <- object[[bestname]]
      P <- unmark(best)
      wP <- attr(best, "weights")
      scopename <- attr(best, "scopename") %orifnull% "grid points"
      coarse.to.fine <- object$coarse.to.fine
    }
    nP <- npoints(P)
    Px <- P$x
    Py <- P$y
    # create space
    ncoef <- length(hom.coef)
    Umat <- matrix(, nrow=nP, ncol=ncoef, dimnames=list(NULL, vec.names))
    if(needed["F"]) Fmat <- matrix(, nrow=nP, ncol=ncoef^2,
                                dimnames=list(NULL, mat.names))
    if(needed["H"]) Hmat <- matrix(, nrow=nP, ncol=ncoef^2,
                                dimnames=list(NULL, mat.names))
    if(needed["X"]) Xmat <- matrix(, nrow=nP, ncol=ncoef^2,
                                   dimnames=list(NULL, mat.names))
    #
    if(verbose)
      splat("Processing", nP, scopename)
    for(k in 1:nP) {
      if(verbose) 
        progressreport(k, nP)
      localwt <- dnorm(Ux - Px[k], sd=sigma) * dnorm(Uy - Py[k], sd=sigma)
      if(ispois) {
        # local score at P[k] for homogeneous model
        locscore <- matrix(localwt * homresid, nrow=1) %*% suf
        Umat[k,] <- locscore
        if(needed["F"] || needed["H"]) {
          # local Fisher information at P[k] for homogeneous model
          if(needed["F"]) {
            locfish <- sumouter(localwt * suf, w * lam)
            Fmat[k,] <- as.vector(locfish)
          }
          # local Hessian at P[k] for homogeneous model
          if(needed["H"]) {
            locHess <- sumouter(suf, w * lam * localwt)
            Hmat[k,] <- as.vector(locHess)
          }
        }
      } else {
        vcl <- vcovlocEngine(internals, localwt, A1dummy=TRUE, new.coef=theta0)
        Umat[k,] <- as.vector(attr(vcl, "Score"))
        if(needed["F"]) Fmat[k,] <- as.vector(attr(vcl, "Fish"))
        if(needed["H"]) Hmat[k,] <- as.vector(attr(vcl, "Grad"))
        if(needed["X"]) Xmat[k,] <- as.vector(attr(vcl, "Cov"))
      }
    }
    if(verbose) splat("\t Done.")
  }

  #' ...... Local variance calculation ...........................
  
  if(saveall || calc.var) {
    ## compute variance of local score at global estimate of theta
    switch(nulltype,
           simple = {
             ## theta is fixed.
             ## variance of local score
             VUmat <- Fmat
           },
           composite = {
             ## theta has been estimated.
             ## Calculate variance of
             ##    {local score evaluated at global estimate of theta}
             ##     allowing for estimation of theta 
             ##     assuming template model is true.
             usepois <- ispois || is.null(Xmat)
             if(!ispois && usepois)
               warning("cross-covariance term not calculated!")
             ## Variance of estimate for *template* model
             ##    var0 = H^{-1} I H^{-1}
             ## where I = information, H = gradient
             ## (H = G if poisson)
             if(usepois) {
               homV <- vcov(homfit)
             } else {
               a <- vcov(homfit, what="all")
               homV <- a$varcov
               homH <- a$internals$hessian %orifnull% a$internals$A1
             }
             ## Make it a flat matrix stack
             homVmat <- as.flat.matrix(homV, nP)
             ## Cross-covariance term(s)
             ## C(s) = H(s) var0 H(s)
             Cmat <- quadform2.flat.matrices(homVmat, Hmat,
                                             c(p,p), c(p,p))
             if(ispois || is.null(Xmat)) {
               VUmat <- Fmat - Cmat
             } else {
               ## D(d) = H(s) H^{-1} cov(local score, global score)
               invhomHmat <- as.flat.matrix(solve(homH), nP)
               Dmat <- multiply3.flat.matrices(Hmat, invhomHmat, Xmat,
                                               c(p,p), c(p,p), c(p,p))
               tDmat <- transpose.flat.matrix(Dmat, c(p,p))
               VUmat <- Fmat + Cmat - Dmat - tDmat
               ## 
               if(spatstat.options('developer')) {
                 ei <- eigenvalues.flat.matrix(VUmat, p)
                 bad <- apply(ei < 0, 1, any)
                 if(any(bad))
                   warning(paste0(signif(100 * mean(bad), 3),
                                  "% of matrices are not positive definite\n"))
               }
             }
           })
  }
  
  # ....   Assemble 'sufficient statistics' ...............
  
  localdata <- list(Fmat   = Fmat,
                    Hmat   = Hmat,
                    Umat   = Umat,
                    VUmat  = VUmat,
                    Lvec   = Lvec,
                    Xmat   = Xmat,
                    ncoef  = ncoef,
                    p      = p,
                    nP     = nP,
                    wts    = wP,
                    hom.coef = hom.coef)

  # ....  Compute the desired map (p-value, test statistic etc) .....

  resultvalues <- HomTestMapEngine(localdata, resultinfo)

  # ....  Wrap up ............................................

  result <- ssf(P, resultvalues)
  class(result) <- c("homtestmap", class(result))
  attr(result, "weights") <- wP
  #
  attr(result, "info") <- resultinfo
  if(saveall) attr(result, "localdata") <- localdata
  # 
  result <- timed(result, starttime=starttime)
  return(result)
}

print.homtestmap <- function(x, ...) {
  with(attr(x, "info"), {
    splat("Homogeneity test map (class 'homtestmap')")
    switch(test,
           score = {
             testname <- "Score Test"
             shortstat <- paste(testname, "Statistic")
           },
           taylor = {
             testname <- paste("Taylor approximation to",
                               "Composite Likelihood Ratio Test")
             shortstat <- "Approximate CLRTS"
             if(adjusting) {
               testname <- paste("Adjusted", testname)
               shortstat <- paste("Adjusted", shortstat)
             }
           },
           likelihood = {
             testname <- "Composite Likelihood Ratio Test"
             shortstat <- "CLRTS"
             if(adjusting) {
               testname <- paste("Adjusted", testname)
               shortstat <- paste("Adjusted", shortstat)
             }
           })
    statname <- paste(testname, "Statistic")
    vname <-
      switch(what,
             components = {
               issqrt <- (test != "likelihood")
               isvector <- issqrt && (p > 1) 
               paste0(if(isvector) "Vector components of " else "",
                      if(issqrt) "signed square root of " else "",
                      paste(testname, "Statistic"))
             },
             statistic = paste(testname, "Statistic"),
             pvalue    = paste("p-values for", testname))
    splat("\nFunction values:", vname)
    if(adjusting)
      switch(ladjust,
             moment = splat(shortstat, "adjusted by first moment matching"),
             PSS = splat(shortstat, "adjusted by Pace-Salvan-Sartori"),
             none = {})
    if(what == "pvalue")
      splat(
        "The p-values were obtained by referring the",
        shortstat, "to the",
        switch(calibrate,
               chisq = paste("chi-squared distribution with", p, "d.f."),
               Satterthwaite =
               "gamma distribution (Satterthwaite approx.)",
               firstmoment =
               "rescaled chi-squared distribution (first moment approx.)")
        )
  })
  cat("\n")
  NextMethod("print")
}


HomTestMapEngine <- function(x, info) {
  #' some of these may be NULL 
  Lvec  <- x$Lvec  # local pseudolikelihood
  Umat  <- x$Umat  # local pseudoscore
  VUmat <- x$VUmat # null variance of local pseudoscore
  Hmat  <- x$Hmat  # Hessian of local log pseudolikelihood
  p     <- x$p     # dimension of parameter
  wts   <- x$wts   # weights associated with sample locations
  hom.coef <- x$hom.coef # global estimate of theta under H0
  ncoef <- length(hom.coef)
  #'
  if(spatstat.options('developer')) {
    splat("homtestmap info:")
    print(unlist(info))
  }
  with(info, {
    ## .....  Now calculate test statistic .....................................
    switch(what,
           components = {
             ## compute components of 'square root' of local test statistic
             switch(test,
                    score = {
                      ## score test statistic
                      ## inverse square root of variance of local score
                      ISVmat <- handle.flat.matrix(VUmat, c(p, p), invsqrtmat)
                      ## premultiply local score
                      rootstat <- multiply2.flat.matrices(ISVmat, Umat,
                                                          c(p, p), c(p, 1))
                      colnames(rootstat) <- names(hom.coef)
                    },
                    taylor = {
                      ## Taylor approximation to local likelihood ratio test
                      ## inverse square root of variance of local Hessian
                      InvsqHmat <- handle.flat.matrix(Hmat, c(p, p), invsqrtmat)
                      ## premultiply local score
                      rootstat <- multiply2.flat.matrices(InvsqHmat, Umat,
                                                          c(p, p), c(p, 1))
                      colnames(rootstat) <- names(hom.coef)
                    },
                    likelihood = {
                      rootstat <- Lvec
                    })
           },
           pvalue = ,
           statistic = {
             ## compute local test statistic
             switch(test,
                    score = {
                      ## score test statistic
                      IVmat <- invert.flat.matrix(VUmat, p)
                      stat <- quadform2.flat.matrices(IVmat, Umat,
                                                      c(p,p), c(1,p))
                    },
                    taylor = {
                      ## Taylor approx of local likelihood ratio test statistic
                      IHmat <- invert.flat.matrix(Hmat, p)
                      stat <- quadform2.flat.matrices(IHmat, Umat,
                                                      c(p,p), c(1,p))
                    },
                    likelihood = {
                      stat <- Lvec
                    })
           })
    ## .....  Adjustment of (approximate) CLRTS .......................
    if(adjusting) {
      switch(ladjust,
             moment = {
               ## first moment matching adjustment
               GinVmat <- solve2.flat.matrices(Hmat, VUmat,
                                               c(p,p), c(p,p))
               E <- trace.flat.matrix(GinVmat, p)
               stat <- (ncoef/E) * stat
             },
             PSS = {
               ## Pace-Salvan-Sartori adjustment
               IVmat <- invert.flat.matrix(VUmat, p)
               scorestat <- quadform2.flat.matrices(IVmat, Umat,
                                                    c(p,p), c(1,p))
               IHmat <- invert.flat.matrix(Hmat, p)
               taylorstat <- quadform2.flat.matrices(IHmat, Umat,
                                                     c(p,p), c(1,p))
               stat <- (scorestat/taylorstat) * stat
             },
             none = { }
             )
    }
    ## .....  Null moments for calibration ...................................
    if(calibrate != "chisq" && test %in% c("likelihood", "taylor")) {
      if(poolmoments) {
        Hpool <- average.flat.matrix(Hmat, c(p, p), wts)
        VUpool <- average.flat.matrix(VUmat, c(p, p), wts)
        GinVpool <- solve(Hpool, VUpool)
        E <- sum(diag(GinVpool))
        if(calibrate == "Satterthwaite")
          V <- 2 * sum(diag(GinVpool %*% GinVpool))
      } else {
        GinVmat <- solve2.flat.matrices(Hmat, VUmat,
                                        c(p,p), c(p,p))
        E <- trace.flat.matrix(GinVmat, p)
        if(calibrate == "Satterthwaite") {
          GinV2mat <- multiply2.flat.matrices(GinVmat, GinVmat,
                                              c(p,p), c(p,p))
          V <- 2 * trace.flat.matrix(GinV2mat, p)
        }
      }
    }
    ## .....  Local p-values ...................................
    if(what == "pvalue") {
      switch(calibrate,
             chisq = {
               pvals <- pchisq(stat, df=ncoef, lower.tail=FALSE)
             },
             Satterthwaite = {
               shape <- E^2/V
               scale <- V/E
               pvals <- pgamma(stat, shape=shape, scale=scale,
                               lower.tail=FALSE)
             },
             firstmoment = {
               shape <- ncoef/2  ## i.e. chi^2, df=ncoef
               scale <- 2 * E/ncoef  ## rescaled so mean =  E
               pvals <- pgamma(stat, shape=shape, scale=scale,
                               lower.tail=FALSE)
             })
    }
    ## .....  Finally assemble the result ...................................
    resultvalues <- switch(what,
                           components = rootstat,
                           statistic  = stat,
                           pvalue     = pvals)
    return(resultvalues)
  })
}


sqmag <- function(x) {
  stopifnot(inherits(x, "ssf"))
  val <- marks(x)
  y <- matrix(rowSums(val^2), ncol=1)
  z <- ssf(unmark(x), y)
  attr(z, "weights") <- attr(x, "weights")
  return(z)
}


homtest <- function(X, ...,
                    nsim=19,
                    test = c("residuals", "score", "taylor", "likelihood"),
                    locations = c("coarse", "fine", "split"),
                    ladjust = NULL,
                    use.fft = NULL,
                    simul = NULL,
                    verbose = TRUE,
                    Xname=NULL) {
  starttime <- proc.time()
  if(is.null(Xname))
    Xname <- short.deparse(substitute(X))
  test <- match.arg(test)
  locations <- match.arg(locations)
  if(is.null(use.fft))
    use.fft <- (locations == "fine") && (test != "likelihood")
  stopifnot(is.ppp(X))
  argh <- list(...)
  # Text string representing the locppm call
  cl <- sys.call()
  if(!is.null(clnames <- names(cl))) {
    discard <- clnames %in% c("nsim", "test", "simul",
                              "verbose", "Xname", "ladjust", "calibrate")
    cl <- cl[!discard]
  }
  cl[[1]] <- as.name("locppm")
  cl <- format(cl)
  testname <- c("Monte Carlo test of homogeneity for", cl)
  testbase <- paste("based on", nsim, "simulations")
  #'
  ladjust <- if(is.null(ladjust)) "PSS" else
             match.arg(ladjust, c("none", "moment", "PSS"))
  ladjname <- if(ladjust == "none") "" else paste0(ladjust, "-adjusted ")
  #' 
  switch(test,
         likelihood = {
           methodname <- paste0("using ", ladjname,
                                "composite likelihood ratio test statistic")
           calcname <- if(use.fft) "computed by FFT" else "computed directly"
           doit <- function(Y, argh, locations, use.fft, ladjust) { 
             fit <- do.call(locppm,
                            c(list(Y),
                              argh,
                              list(vcalc="lik",
                                   locations=locations,
                                   use.fft=use.fft,
                                   verbose=FALSE)))
             h <- homteststat(fit, test="likelihood", ladjust=ladjust,
                              use.fft=use.fft, verbose=FALSE)
             return(h)
           }
         },
         taylor = {
           methodname <- paste0("using ", ladjname, "Taylor approximation to ",
                                "composite likelihood ratio test statistic")
           calcname <- if(use.fft) "computed by FFT" else "computed directly"
           doit <- function(Y, argh, locations, use.fft, ladjust) { 
             fit <- do.call(locppm,
                            c(list(Y),
                              argh,
                              list(vcalc="hom",
                                   locations=locations,
                                   use.fft=use.fft,
                                   verbose=FALSE)))
             h <- homteststat(fit, test="taylor", ladjust=ladjust,
                              use.fft=use.fft, verbose=FALSE)
             return(h)
           }
         },
         score = {
           methodname <- "using score test statistic"
           calcname <- if(use.fft) "computed by FFT" else "computed directly"
           doit <- function(Y, argh, locations, use.fft, ladjust) { 
             fit <- do.call(locppm,
                            c(list(Y),
                              argh,
                              list(vcalc="hom",
                                   locations=locations,
                                   use.fft=use.fft,
                                   verbose=FALSE)))
             h <- homteststat(fit, test="score",
                              use.fft=use.fft, verbose=FALSE)
             return(h)
           }
         },
         residuals = {
           methodname <- "using score residuals"
           calcname <- "computed directly"
           doit <- function(Y, argh, locations, use.fft, ladjust) {
             fit <- do.call(ppm, append(list(Y), argh))
             res <- residuals(fit, "score")
             smr <- Smooth(res)
             if(is.im(smr)) smr <- list(smr)
             sqr <- lapply(smr, "^", e2=2)
             sumsqr <- Reduce("+", sqr)
             h <- mean(sumsqr)
             attr(h, "nsample") <- n.quad(quad.ppm(fit))
             return(h)
           }
         })

  if(verbose) cat("Computing observed value...")
  obs <- doit(X, argh, locations=locations, use.fft=use.fft, ladjust=ladjust)
  if(verbose) cat("Done.\n")
  # extract info
  nsample <- attr(obs, "nsample")
  if(!is.null(nsample))
    calcname <- paste(calcname, "on", nsample, "sample points")
  # save for use in return value
  obs <- as.numeric(obs)
  names(obs) <- "S"
  # fit 'homogeneous' model
  homfit <- ppm(X, ...)
  # generate the simulated patterns
  if(is.null(simul)) {
    Xsim <- simulate(homfit, nsim=nsim, progress=verbose)
  } else if(is.expression(simul)) {
    Xsim <- list()
    if(verbose) cat(paste("Generating", nsim,
                          "simulated realizations by evaluating expression..."))
    for(i in 1:nsim)
      Xsim[[i]] <- eval(simul)
    if(!all(unlist(lapply(Xsim, is.ppp))))
      stop("Evaluating the expression did not yield a point pattern")
  } else if(is.list(simul)) {
    Xsim <- simul
    if(!all(unlist(lapply(Xsim, is.ppp))))
      stop("Entries in the list should be point patterns")
  } else stop("Unrecognised format of argument simul")
  # evaluate the statistic
  sim <- numeric(nsim)
  if(verbose) cat("Computing test statistics...")
  for(i in 1:nsim) {
    if(verbose) progressreport(i, nsim)
    sim[i] <- doit(Xsim[[i]], argh, locations=locations,
                   use.fft=use.fft, ladjust=ladjust)
  }
  if(verbose) cat("Done.\n")
  pval <- (1 + sum(sim >= obs))/(nsim+1)
  result <- list(statistic = obs,
                 p.value = pval,
                 alternative = "inhomogeneous",
                 method = c(testname, methodname, calcname, testbase),
                 data.name = Xname,
                 simulated = sim)
  class(result) <- "htest"
  attr(result, "homfit") <- homfit
  result <- timed(result, starttime=starttime)
  return(result)
}
 
bw.locppm <- function(...,
                      method = c("fft", "exact", "taylor"), 
                      srange = NULL, ns=9, sigma = NULL,
                      additive = TRUE,
                      verbose=TRUE) {
  starttime <- proc.time()
  parenv <- sys.parent()

  method <- match.arg(method)
  additive <- additive && (method == "taylor")

  # fit homogeneous model
  if(verbose) cat("Fitting homogeneous model... ")  
  homfit <- eval(substitute(ppm(..., forcefit=TRUE)), envir=parenv)
  if(is.multitype(homfit))
    stop("Sorry, cannot handle marked point processes yet")
  p <- length(coef(homfit))
  ispois <- is.poisson(homfit)

  if(verbose) cat("Computing model properties... ")
  
  # extract quadrature info
  X <- data.ppm(homfit)
  Q <- quad.ppm(homfit)
  U <- union.quad(Q)
  wQ <- w.quad(Q)
  Z <- is.data(Q)
  nX <- npoints(X)
  nU <- npoints(U)
  Xindex <- seq_len(nU)[Z]
  xU <- U$x
  yU <- U$y
  
  # save internal structures for use in prediction
  gfit <- getglmfit(homfit)
  gdat <- getglmdata(homfit)

  mm <- model.matrix(homfit)
  mmX <- mm[Z, , drop=FALSE]

  lambda0  <- fitted(homfit)
  coef0    <- coef(homfit)
  coef0mat <- matrix(coef0, nU, length(coef0), byrow=TRUE)
  
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

  # specify fitting options and pre-compute relevant data
  switch(method,
         fft = {
           ## Taylor approximation coefficients and gradient
           ## on quadrature points
           opt.taylor <- locppmOptions(cg1=TRUE, gg1=TRUE)
           #' precompute some values
           internals <- list(lambda=lambda0,
                             hom.coef=coef0,
                             mom=mm)
         },
         exact = {
           ## fit on quadrature points
           opt.quad <- locppmOptions(cg=TRUE)
           ## leave-one-out calculation on data points
           opt.data <- locppmOptions(xg=TRUE)
           internals <- NULL
         },
         taylor = {
           ## fit on quadrature points
           opt.quad <- locppmOptions(cg=TRUE)
           ## leverage approximation on data points using Hessian
           opt.data <- locppmOptions(Vg=TRUE)
           ## precompute internal data for variance estimates
           internals <- getvcinternals(homfit, verbose=verbose)
         })

  if(!ispois && method != "exact") {
    # pre-compute data for leverage approximation
    if(verbose) cat("Leverage... ")
    internals$coef <- internals$hom.coef
    a <- ppmInfluence(homfit, what=c("increments"), precomputed=internals)
    # change in canonical statistic for X[i] when X[j] is deleted
    ddS <- a$increm$ddS[ Z, Z, , drop=FALSE]
    # change in integral increment at U[i] when X[j] is deleted
    ddSincrements <- wQ * a$increm$ddSintegrand[ , Z, , drop=FALSE]
  }
  if(verbose) cat("Done.\n")

  ## allocate storage
  datasum <- dof <- intlam <- numeric(ns)

  ## fit using each value of sigma
  if(verbose) {
    cat(paste("Assessing", ns, "values of sigma... "))
    pstate <- list()
  }
  for(k in 1:ns) {
    sigk <- sigma[k]
    kernel0 <- 1/(2 * pi * sigk^2)
    switch(method,
           fft = {
             #' fit and Hessian approximated on quadrature points
             lpe <- locppmEngine(homfit, sigk, U, opt=opt.taylor,
                                 scopename="quadrature points",
                                 verbose=FALSE, internals=internals)
             #' fitted [conditional] intensity at each quadrature point
             lambda <- locppmPredict(homfit, marks(lpe$cg1),
                                     precomputed=internals, details=TRUE)
             d.eta <- attr(lambda, "d.eta")
             #' Hessian at data points
             hQ <- lpe[[ "gg1" ]]
             hX <- marks(hQ)[Z,,drop=FALSE]
             #' invert the Hessians
             vX <- invert.flat.matrix(hX, p)
           },
           exact = {
             #' fit on quadrature points
             lpe.quad <- locppmEngine(homfit, sigk, U, opt=opt.quad,
                                      scopename="quadrature points",
                                      verbose=FALSE, internals=internals)
             ## compute fitted conditional intensity at each quadrature point
             lambda <- locppmPredict(homfit, marks(lpe.quad$cg), 
                                     precomputed=internals, details=TRUE)
             #' leave-one-out fit on data points
             lpe.data <- locppmEngine(homfit, sigk, X, opt=opt.data,
                                      scopename="data points",
                                      scopeindex=Xindex,
                                      verbose=FALSE, internals=internals)
             ## leave-one-out estimates at data points
             lambdaXminus <- locppmPredict(homfit, marks(lpe.data$xg),
                                           index=Xindex,
                                           precomputed=internals, details=TRUE)
           },
           taylor = {
             #' fit on quadrature points
             lpe.quad <- locppmEngine(homfit, sigk, U, opt=opt.quad,
                                      scopename="quadrature points",
                                      verbose=FALSE, internals=internals)
             ## compute fitted conditional intensity at each quadrature point
             lambda <- locppmPredict(homfit, marks(lpe.quad$cg), 
                                     precomputed=internals, details=TRUE)
             d.eta <- attr(lambda, "d.eta")
             ## leverage approximation on data points
             lpe.data <- locppmEngine(homfit, sigk, X, opt=opt.data,
                                      scopename="data points",
                                      scopeindex=Xindex,
                                      verbose=FALSE, internals=internals)
             ## leverage approximation uses inverse Hessian at data points
             vX <- lpe.data[[ "Vg" ]]
             vX <- marks(vX)
           })
    
    # compute terms in cross-validation criterion
    datasum[k] <- sum(log(lambda[Z])) 
    intlam[k] <- if(anyNA(lambda)) Inf else sum(lambda * wQ)
    #'
    if(method == "exact") {
      dof[k] <- datasum[k] - sum(log(lambdaXminus))
    } else {
      ## compute 'leverage' approximation 
      if(ispois) {
        ##    log(lambda(x_i)) - log(lambda_{-i}(x_i))
        ##    ~ kernel(0) Z(x_i) V(x_i) Z(x_i)^T
        leve <- kernel0 * quadform2.flat.matrices(vX, mmX, c(p, p), c(1, p))
      } else {
        ##    log(lambda(x_i)) - log(lambda_{-i}(x_i))
        ##    ~ kernel(0) Z(x_i) V(x_i) Y(x_i)^T
        ## compute Y(x_i)
        dScore <- kernel0 * mmX
        lrat <- exp(d.eta) # adjustment factor for lambda
        seqX <- 1:nX
        for(j in seqX) {
          weiU <-
            dnorm(xU, mean=xU[j], sd=sigk) * dnorm(yU, mean=yU[j], sd=sigk)
          weiX <- weiU[seqX]
          A <- apply(weiX * ddS[, j, , drop=FALSE], c(2,3), sum)
          B <- apply(weiU * lrat * ddSincrements[ , j, , drop=FALSE],
                     c(2,3), sum)
          dScore[j,] <- dScore[j,] + A + B
        }
        leve <- bilinear3.flat.matrices(mmX, vX, dScore, c(1,p), c(p,p), c(1,p))
      }
      dof[k] <- if(!additive) sum(leve) else
                if(all(leve < 1)) -sum(log(1-leve)) else Inf
    }
    #' end of loop
    if(verbose) pstate <- progressreport(k, ns, state=pstate)
  }

  # cross-validation criterion
  gcv <- datasum - dof - intlam
  result <- bw.optim(gcv, sigma, iopt=which.max(gcv), cvname="cv",
                     dof=dof, intlam=intlam, datasum=datasum)
  attr(result, "method") <- method
  attr(result, "resolution") <- summary(quad.ppm(homfit))[["resolution"]]
  timed(result, starttime=starttime)
}

