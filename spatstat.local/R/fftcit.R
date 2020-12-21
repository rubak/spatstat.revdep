#
# fftcit.R
#
# Calculations for loccit that are possible using Fast Fourier Transform
#
# $Revision: 2.15 $ $Date: 2013/10/31 07:33:38 $

loccitFFT <- local({
  
  # Dependence between calculations
  #  ("Score", "Hessian" etc represent blocks of code in the main function)
  .Score     <- "Score"
  .Hessian   <- "Hessian"
  .HessMat   <- c(.Hessian, "HessMat")
  .ScoreMat  <- c(.Score, "ScoreMat")
  .InvHess   <- c(.HessMat, "InvHess")
  .Param     <- c(union(.HessMat, .ScoreMat), "Param")
  .Delta     <- c(.Param, "Delta")
  .Lambda    <- c(.Param, "Lambda") 
  .Tgrad     <- c(union(.InvHess, .Param), "Tgrad")
  .PalmLik   <- c(union(.Score, .Lambda), "PalmLik")
  .Influence  <- c(.PalmLik, "Influence")

  # Table mapping user options to calculations required
  Dependence.Table <-
    list(param     = .Param,
         score     = .Score,
         grad      = .Hessian,
         invgrad   = .InvHess,
         delta     = .Delta,
         lambda    = .Lambda,
         tgrad     = .Tgrad,
         influence  = .Influence,
         plik      = .PalmLik)

  opties <- names(Dependence.Table)

  # corresponding output names
  
  FullNames <- 
    c("parameters",
      "score",
      "gradient",
      "invgradient",
      "delta",
      "intensity",
      "tgradient",
      "influence",
      "logplik")
  
  sumby <- function(X, Y) {
    v <- by(data=X, INDICES=Y, FUN=colSums)
    if(any(probs <- unlist(lapply(v, is.null)))) {
      zeroes <- rep(0, ncol(X))
      for(i in which(probs)) v[[i]] <- zeroes
    }
    a <- simplify2array(v)
    b <- if(length(dim(a)) > 1) t(a) else matrix(a, ncol=1)
    b
  }

  # 2D gaussian kernel as a function of interpoint distance
  gauss <- function(d, sigma=1) { exp(-d^2/(2 * sigma^2))/(2 * pi * sigma^2)}

  # convert array to flat matrix
  flatten <- function(x) { matrix(as.vector(x), nrow=dim(x)[1])} 

  # weed out values arising from divide-by-zero
  clampfinite <- function(x, hi=.Machine$double.xmax, lo=-hi) {
    x[is.nan(x)] <- 0
    x[x == Inf] <- hi
    x[x == -Inf] <- lo
    x
  }
  
  loccitFFT <- function(fit0, sigma, rmax, ...,
                        base.trendcoef=NULL,
                        base.cluspar=NULL,
                        base.lambda=fitted(fit0, new.coef=base.trendcoef),
                        base.lambdaim=NULL,
                        clusters=NULL,
                        hompoisfit=as.ppm(fit0), 
                        what = "param",
                        do.trend = TRUE,
                        do.clusters = is.kppm(fit0),
                        calcopt = list(), 
                        verbose=TRUE) {
    starttime <- proc.time()

    stopifnot(is.kppm(fit0) || is.ppm(fit0) || inherits(fit0, "locppm"))
    
    # what to calculate
    hit <- pmatch(what, opties)
    if(any(bad <- is.na(hit)))
      warning(paste("Unrecognised option: what=", commasep(sQuote(what[bad]))))
    what <- opties[hit[!bad]]
    needed <- unique(unlist(Dependence.Table[what]))

    # options for calculations
    calcopt.default <- list(use.D2 = FALSE,
                            fishscore = FALSE,
                            refit.for.palmlik = FALSE,
                            refit.for.influence = FALSE,
                            dropidpairs=FALSE)
    if(any(unknown <- !(names(calcopt) %in% names(calcopt.default))))
      stop(paste("Unrecognised",
                 ngettext(sum(unknown), "option", "options"),
                 commasep(sQuote(names(calcopt)[unknown]))))
    calcopt <- resolve.defaults(calcopt, calcopt.default)
    use.D2              <- calcopt$use.D2
    fishscore           <- calcopt$fishscore
    refit.for.palmlik   <- calcopt$refit.for.palmlik
    refit.for.influence <- calcopt$refit.for.influence
    dropidpairs         <- calcopt$dropidpairs

    # extract information about Poisson trend model
    stopifnot(is.ppm(hompoisfit))
    # Fitted trend parameters
    homtrendcoef <- coef(hompoisfit)
    ntrend <- length(homtrendcoef)
    # basic data
    Q <- quad.ppm(hompoisfit)
    P <- union.quad(Q)
    nP <- npoints(P)
    isdata <- is.data(Q)
    wQ <- w.quad(Q)
    X <- Q$data
    nX <- npoints(X)
    if(do.trend) {
      mom <- model.matrix(hompoisfit)
      momX <- mom[isdata, , drop=FALSE]
    }

    # make a template mask for images
    W <- as.mask(as.owin(P), ...)
    Wmat <- as.matrix(W)
    npixels <- length(Wmat)

    # .................. baseline trend ....................
    # validate base.trendcoef
    if(is.null(base.trendcoef)) {
      # default: homogeneous coefficients
      base.trendcoef <- homtrendcoef
      base.lambda <- fitted(fit0)
    } else {
      # validate dimensions 
      if(is.matrix(base.trendcoef)) {
        check.nmatrix(base.trendcoef, ntrend, things="trend parameters",
                      squarematrix=FALSE, matchto="ncol", naok=TRUE)
        check.nmatrix(base.trendcoef, nP, things="quadrature points",
                      squarematrix=FALSE, naok=TRUE)
        base.lambda <- locppmPredict(hompoisfit, base.trendcoef)
      } else {
        check.nvector(base.trendcoef, ntrend,
                      things="trend parameters", naok=TRUE)
        base.lambda <- fitted(fit0, new.coef=base.trendcoef)
      }
    }
    # coerce to matrix
    if(!is.matrix(base.trendcoef))
      base.trendcoef <- matrix(base.trendcoef,
                               nP, length(base.trendcoef),
                               dimnames=list(NULL, names(base.trendcoef)),
                               byrow=TRUE)

    # Baseline intensity at data points
    base.lambdaX <- base.lambda[isdata]
    
    # .................. baseline clusters ....................
    if(do.clusters) {
      # extract information about cluster model
      if(is.null(clusters)) {
        # default: cluster model is fit0
        if(!is.kppm(fit0)) 
          stop("If fit0 is not a kppm object, clusters must be specified")
        clusters <- fit0$clusters 
      }
      info <- spatstatClusterModelInfo(clusters)
      extra <- extraClusterModelInfo(clusters)
      paco      <- extra$pcf  # parallelised version of pair correlation
      # conversion to/from canonical parameters (may be NULL)
      par2theta <- extra$par2theta
      theta2par <- extra$theta2par
      # derivatives (parallelised)
      Dlogpcf  <- extra$DlogpcfDtheta
      if(is.null(Dlogpcf))
        stop(paste("Derivatives not available for", sQuote(clusters)))
      if(use.D2) {
        D2logpcf <- extra$D2logpcfDtheta
        if(is.null(D2logpcf))
          stop(paste("Second derivatives not available for", sQuote(clusters)))
      }

      # extract cluster parameters of baseline model
      if(is.null(base.cluspar)) {
        # default: baseline is fit0
        if(!is.kppm(fit0)) 
          stop("If fit0 is not a kppm object, base.cluspar must be specified")
        base.cluspar <- fit0$par
      } else {
        # validate 'base.cluspar'
        if(is.matrix(base.cluspar)) 
          check.nmatrix(base.cluspar, nP, things="quadrature points",
                        squarematrix=FALSE, naok=TRUE)
      }
      # coerce to matrix
      if(!is.matrix(base.cluspar)) 
        base.cluspar <- matrix(base.cluspar,
                               nrow=nP, ncol=length(base.cluspar),
                               dimnames=list(NULL, names(base.cluspar)),
                               byrow=TRUE)
      ncluspar <- ncol(base.cluspar) 
      # Transform to canonical parameters
      base.cluscoef <- applymaps(par2theta, base.cluspar)
    }
    
    # all canonical parameters
    base.theta <- cbind(if(do.trend) base.trendcoef else NULL,
                        if(do.clusters) base.cluscoef else NULL)
    p <- ncol(base.theta)
    vec.names <- colnames(base.theta)
    mat.names <- as.vector(outer(vec.names, vec.names, paste, sep="."))

    # close pairs of data points
    clX <- closepairs(X, rmax)
    I <- clX$i
    J <- clX$j
    dIJ <- clX$d
    nIJ <- length(I)
    facI <- factor(I, levels=1:nX)
    facJ <- factor(J, levels=1:nX)

    # close pairs of (data, quadrature) points
    clXP <- crosspairs(X, P, rmax)
    IX <- clXP$i
    JP <- clXP$j
    dXIPJ <- clXP$d
    if(dropidpairs) {
      # exclude identical pairs - assume I and J are comparable
      ne <- (IX != JP)
      IX <- IX[ne]
      JP <- JP[ne]
      dXIPJ <- dXIPJ[ne]
    }
    #
    nXIPJ <- length(IX)
    lamPJ <- base.lambda[JP]
    wPJ <- wQ[JP]
    facIX <- factor(IX, levels=1:nX)
    facJP <- factor(JP, levels=1:nP)
    
    if("Score" %in% needed) {
      # 'data' contribution to the score from X[J]
      momXJ <- if(do.trend) momX[J, , drop=FALSE] else NULL
      if(do.clusters) {
        base.cluscoefJ <- base.cluscoef[J, , drop=FALSE]
        zetaIJ <- Dlogpcf(base.cluscoefJ, dIJ)
      } else zetaIJ <- NULL
      zetaIJ <- cbind(momXJ, zetaIJ)     # = zeta_P(X[J] | X[I])
      score.data <- sumby(zetaIJ, facJ)  # = \sum_i zeta_P(X[j] | X[i])
      # 'background' contribution to the score from P[JP]
      momPJ    <- if(do.trend) mom[JP, , drop=FALSE] else NULL
      if(do.clusters) {
        base.cluscoefPJ <- base.cluscoef[JP, , drop=FALSE]
        zetaXIPJ <- Dlogpcf(base.cluscoefPJ, dXIPJ)
      } else zetaXIPJ <- NULL
      zetaXIPJ <- cbind(momPJ, zetaXIPJ)
      if(do.clusters) {
        base.clusparPJ <- base.cluspar[JP, , drop=FALSE]
        palmXIPJ <- lamPJ * paco(base.clusparPJ, dXIPJ)
      } else palmXIPJ <- lamPJ
      score.back.XIPJ <- zetaXIPJ * palmXIPJ * wPJ
      score.back <- sumby(score.back.XIPJ, facJP)
      # put together
      scoreP <- -score.back
      scoreP[isdata, ] <- scoreP[isdata, ] + score.data
      colnames(scoreP) <- vec.names
      # apply kernel
      U <- density(P, weights=scoreP, sigma=sigma, edge=FALSE, w=W)
    } else U <- NULL

    if("ScoreMat" %in% needed) {
      Umat <- imagelist2matrix(U)
    } else Umat <- NULL
    
    if("Hessian" %in% needed) {
      if(!fishscore) {
        # calculate the Hessian
        if(use.D2 && do.clusters) {
          # 'data' contribution to the Hessian (as flat matrix)
          # involves cluster parameters only.
          kappaIJ <- flatten(D2logpcf(base.cluscoefJ, dIJ))
          # flatten array
          kappaIJ <- matrix(as.vector(kappaIJ), nrow=nIJ)
          # sum over i
          kappa.data <- sumby(kappaIJ, facJ)
          if(do.trend) {
            # indices of the submatrix of the Hessian
            # related to the cluster parameters
            m <- matrix(, p, p)
            clusind <- which(as.vector(row(m) > ntrend & col(m) > ntrend))
            hess.data <- matrix(0, nX, p^2)
            hess.data[, clusind] <- kappa.data
          } else hess.data <- kappa.data
        } else hess.data <- 0
        # 'background' contribution to the Hessian
        outerzetaXIPJ <- outersquare.flat.vector(zetaXIPJ)
        if(use.D2 && do.clusters) {
          kapXIPJ <- flatten(D2logpcf(base.cluscoefPJ, dXIPJ))
          kapXIPJ <- matrix(as.vector(kapXIPJ), nrow=nXIPJ)
          if(do.trend) {
            kappaXIPJ <- matrix(0, nXIPJ, p^2)
            kappaXIPJ[, clusind] <- kapXIPJ
          } else kappaXIPJ <- kapXIPJ
        } else kappaXIPJ <- 0
        hess.back <- sumby((kappaXIPJ + outerzetaXIPJ) * palmXIPJ * wPJ,
                             facJP)
        # put together
        hessP <- hess.back
        hessP[isdata,] <- hessP[isdata, ] - hess.data
        colnames(hessP) <- mat.names
        # apply kernel
        H <- density(P, weights=hessP, sigma=sigma, edge=FALSE, w=W)
      } else {
        # Analogue of Fisher scoring - take E of Hessian
        #
        # close pairs of (quadrature, quadrature) points
        clP <- closepairs(P, rmax)
        IPP <- clP$i
        JPP <- clP$j
        dPIPJ <- clP$d
        nPIPJ <- length(IPP)
        lamPPI <- base.lambda[IPP]
        lamPPJ <- base.lambda[JPP]
        wPPJ <- wQ[JPP]
        facIPP <- factor(IPP, levels=1:npoints(P))
        facJPP <- factor(JPP, levels=1:npoints(P))
        # 
        momPPJ   <- if(do.trend) mom[JPP, , drop=FALSE] else NULL
        if(do.clusters) {
          base.cluscoefPPJ <- base.cluscoef[JPP, , drop=FALSE]
          zetaPIPJ <- Dlogpcf(base.cluscoefPPJ, dPIPJ)
        } else zetaPIPJ <- NULL
        zetaPIPJ <- cbind(momPPJ, zetaPIPJ)
        outerzetaPIPJ <- outersquare.flat.vector(zetaPIPJ)
        # second moment intensity
        lambda2PIPJ <- base.lambda[IPP] * base.lambda[JPP]
        if(do.clusters) {
          base.clusparPPJ <- base.cluspar[JPP, , drop=FALSE]
          lambda2PIPJ <- lambda2PIPJ * paco(base.clusparPPJ, dPIPJ)
        }
        #
        hessP <- sumby(outerzetaPIPJ * lambda2PIPJ * wQ[IPP],
                       facJPP)
        colnames(hessP) <- mat.names
        # apply kernel
        H <- density(P, weights=wQ * hessP, sigma=sigma, edge=FALSE, w=W)
      }
    } else H <- NULL

    if("HessMat" %in% needed) {
      # convert to flat matrix
      Hmat <- imagelist2matrix(H)
    } else Hmat <- NULL
    
    if("InvHess" %in% needed) {
      # Compute inverse of negative local Hessian as matrix
      invHmat <- invert.flat.matrix(Hmat, p)
      # convert back to images
      invH <- matrix2imagelist(invHmat, W)
      names(invH) <- mat.names
    } else invH <- invHmat <- NULL
    #

    if("Param" %in% needed) {
      # evaluate Taylor approximation of locally fitted coefficients
      if(!is.null(invHmat)) {
        # use inverse Hessian: delta(theta) = invH %*% U
        Deltamat <- multiply2.flat.matrices(invHmat, Umat, c(p,p), c(p, 1))
      } else if(!is.null(Hmat)) {
        # use Hessian: delta(theta) = solve(H, U)
        Deltamat <- solve2.flat.matrices(Hmat, Umat, c(p, p), c(p, 1))
      } else stop("Internal error: coefficient calculation needs Hessian")
      # If inversion is ill-conditioned, set delta(theta) = 0
      if(any(fail <- !complete.cases(Deltamat)))
        Deltamat[fail, ] <- 0
      # Baseline coefficients
      Theta0 <- nnmark(P %mark% base.theta, w=W)
      Theta0mat <- imagelist2matrix(Theta0)
      # update by Taylor expansion/Fisher scoring
      Thetamat <- clampfinite(Theta0mat + Deltamat)
      # convert back to images
      Theta <- matrix2imagelist(Thetamat, W)
      names(Theta) <- vec.names
      if(!do.clusters) {
        Param <- Theta
      } else {
        # extract cluster coefficients
        ClusCoef <- if(!do.trend) Theta else Theta[ntrend + 1:ncluspar]
        # convert back to original parameters
        ClusPar <- applymaps(theta2par, ClusCoef)
        # ensure values are finite
        ClusPar <- as.listof(lapply(ClusPar,
                                    function(z) eval.im(clampfinite(z))))
        # put back
        Param <- if(!do.trend) ClusPar else append(Theta[1:ntrend], ClusPar)
      }
    } else Param <- NULL
    
    if("Delta" %in% needed) {
      Delta <- matrix2imagelist(Deltamat, W)
      names(Delta) <- vec.names
    } else Delta <- NULL
    
    if("Tgrad" %in% needed)
      warning("tgradient not yet implemented")
    Tgrad <- NULL

    if("Lambda" %in% needed) {
      if(!do.trend)
        stop("Lambda cannot be computed when do.trend=FALSE")
      # Change in trend coefficients
      DeltaPhimat <- Deltamat[, 1:ntrend, drop=FALSE]
      # Canonical sufficient statistics for trend
      Mom <- model.images(hompoisfit, W=W)
      Mommat <- imagelist2matrix(Mom)
      # Intensity for homogeneous model
      if(is.null(base.lambdaim)) {
        base.lambdaim <- if(inherits(fit0, "locppm") || !is.matrix(base.trendcoef)) {
          predict(fit0, locations=W, new.coef=base.trendcoef)
        } else {
          Smooth(ssf(P, locppmPredict(as.ppm(fit0), base.trendcoef)))
        }
      } else {
        stopifnot(is.im(base.lambdaim))
        stopifnot(compatible(as.im(W), base.lambdaim))
      }
      Lam0vec <- as.vector(as.matrix(base.lambdaim))
      # Intensity of local model
      lam <- Lam0vec * exp(clampfinite(rowSums(DeltaPhimat * Mommat)))
      lam <- clampfinite(lam)
      lam <- matrix(lam, ncol=1)
      Lambda <- matrix2imagelist(lam, W)
      names(Lambda) <- "lambda"
    } else Lambda <- NULL

    if("PalmLik" %in% needed) {
      # Calculate Palm likelihood
      if(refit.for.palmlik) {
        # use updated intensity
        lamP <- sample.imagelist(Lambda, P)[,1]
      } else {
        # use baseline intensity
        lamP <- base.lambdaim[P]
      }
      lamX <- lamP[isdata]
      # compute pair correlation terms
      if(do.clusters) {
        # fitted cluster parameters
        Psi <- if(do.trend) Param[-(1:ntrend)] else Param
        # fitted cluster parameters at quadrature points
        parP <- sample.imagelist(Psi, P)
        parX <- parP[isdata, , drop=FALSE]
        # evaluate pair correlation for (x_i, x_j)
        # using parameters fitted at x_j
        parXJ <- parX[J, , drop=FALSE]
        rhoIJ <- paco(parXJ, dIJ)
        # evaluate pair correlation for (x_i, u)
        # using parameters fitted at u
        parXIPJ <- parP[JP, , drop=FALSE]
        rhoXIPJ <- paco(parXIPJ, dXIPJ)
      } else rhoIJ <- rhoXIPJ <- 1
      # 'data' contribution
      plik.data <- sum(log(lamX)[J]) + sum(log(rhoIJ))
      # 'background' contribution
      plik.back <- sum(lamP[JP] * rhoXIPJ * wPJ)
      # put together
      LogPalmLik <- clampfinite(plik.data - plik.back)
    } else LogPalmLik <- NULL

    if("Influence" %in% needed) {
      # Influence calculation
      if(!refit.for.influence) {
        # Use score and Hessian for reference model, from above
        ZetaIJ <- zetaIJ
        Score.back.XIPJ <- score.back.XIPJ
        HP <- sample.imagelist(H, P)
      } else {
        # Recompute score and Hessian using fitted parameters
        if(!(do.clusters && do.trend))
          stop("Influence calculation requires do.clusters=do.trend=TRUE")
        # extract fitted parameters at quadrature points
        fittheta <- sample.imagelist(Theta, P)   # canonical parameters
        fitparam <- sample.imagelist(Param, P)   # natural parameters
        fittrendcoef <- fittheta[, 1:ntrend, drop=FALSE]
        fitcluscoef  <- fittheta[, -(1:ntrend), drop=FALSE]
        fitcluspar   <- fitparam[, -(1:ntrend), drop=FALSE]
        # 'data' contribution to the score from X[J]
        #   using fitted coefficients at X[J]
        fitcluscoefJ <- fitcluscoef[J, , drop=FALSE]
        ZetaIJ <- clampfinite(Dlogpcf(fitcluscoefJ, dIJ))
        ZetaIJ <- cbind(momXJ, ZetaIJ)
        #
        # 'background' contribution to the score from P[JP]
        #
        fitcluscoefJP <- fitcluscoef[JP, , drop=FALSE]
        ZetaXIPJ <- clampfinite(Dlogpcf(fitcluscoefJP, dXIPJ))
        ZetaXIPJ <- cbind(momPJ, ZetaXIPJ)
        LamPJ <- lamP[JP]
        fitclusparPJ <- fitcluspar[JP, , drop=FALSE]
        PalmXIPJ <- LamPJ * paco(fitclusparPJ, dXIPJ)
        Score.back.XIPJ <- ZetaXIPJ * (PalmXIPJ * wPJ)
        # Also compute (expected) Hessian at data points using fitted parameters
        if(!fishscore) {
          # Hessian
          if(use.D2) {
            fitcluscoefJ <- fitcluscoef[J, , drop=FALSE]
            KappaIJ <- clampfinite(flatten(D2logpcf(fitcluscoefJ, dIJ)))
            # sum over i
            Kappa.data <- sumby(KappaIJ, facJ)
            Hess.data <- matrix(0, nX, p^2)
            Hess.data[, clusind] <- Kappa.data
          } else Hess.data <- 0
          # 'background' contribution to the Hessian
          OuterzetaXIPJ <- outersquare.flat.vector(ZetaXIPJ)
          if(use.D2) {
            fitcluscoefJP <- fitcluscoef[JP, , drop=FALSE]
            KapXIPJ <- D2logpcf(fitcluscoefJP, dXIPJ)
            KappaXIPJ <- matrix(0, nXIPJ, p^2)
            KappaXIPJ[, clusind] <- clampfinite(flatten(KapXIPJ))
          } else KappaXIPJ <- 0
          Hess.back <-
            sumby(clampfinite((KappaXIPJ + OuterzetaXIPJ) * PalmXIPJ * wPJ),
                  facJP)
          # put together
          HessP <- Hess.back
          HessP[isdata,] <- HessP[isdata, ] - Hess.data
          colnames(HessP) <- mat.names
          # apply kernel
          HP <- density(P, weights=clampfinite(HessP),
                        at="points", leaveoneout=FALSE,
                        sigma=sigma, edge=FALSE, w=W)
          HP <- as.matrix(HP)
        } else {
          # Analogue of Fisher scoring - take E of Hessian
          #
          fitcluscoefPPJ <- fitcluscoef[JPP, , drop=FALSE]
          ZetaPIPJ <- clampfinite(Dlogpcf(fitcluscoefPPJ, dPIPJ))
          ZetaPIPJ <- cbind(momPPJ, ZetaPIPJ)
          OuterzetaPIPJ <- outersquare.flat.vector(ZetaPIPJ)
          # second moment intensity
          fitclusparPPJ <- fitcluspar[JPP, , drop=FALSE]
          Lambda2PIPJ <- lamP[IPP] * lamP[JPP] * paco(fitclusparPPJ, dPIPJ)
          #
          HessP <- sumby(OuterzetaPIPJ * Lambda2PIPJ * wQ[IPP],
                         facJPP)
          colnames(hessP) <- mat.names
          # apply kernel
          HP <- density(P, weights=wQ * clampfinite(hessP),
                        at="points", leaveoneout=FALSE, 
                        sigma=sigma, edge=FALSE, w=W)
          HP <- as.matrix(HP)
        }
      }
      # Score.data[j] = \sum_i zeta_P(X[j] | X[i])
      #     Derivative of \alpha(X[j] | X, theta) w.r.t. theta
      Score.data <- sumby(ZetaIJ, facJ)
      colnames(Score.data) <- vec.names
      # delta.score[j] = w(0) score.data[j] 
      #            + \sum_i w(X[j]-X[i]) zeta_P(X[i] | X[j])
      #            - \int w(X[j] - u) zeta_P(u | X[j]) lambda_P(u | X[j]) du
      #      Total effect on score of deleting X[j]
      #
      Delta.score <-
        gauss(0, sigma) * Score.data + 
          sumby(gauss(dIJ, sigma) * ZetaIJ, facI) -
            sumby(gauss(dXIPJ, sigma) * Score.back.XIPJ, facIX)
      colnames(Delta.score) <- vec.names
      # pack up
      Influence <- list(score.data=Score.data,
                        delta.score=Delta.score,
                        HP=HP)
    } else Influence <- NULL

    # finally pack up
    answer <- list(parameters  = Param,
                   score       = U,
                   gradient    = H,
                   invgradient = invH,
                   delta       = Delta, 
                   intensity   = Lambda,
                   tgradient   = Tgrad,
                   influence   = Influence,
                   logplik     = LogPalmLik)

    # return only those quantities that were commanded
    indx <- pmatch(what, opties) # there are no NA's
    answer <- answer[FullNames[indx]]
  
    return(timed(answer, starttime=starttime))
  }

  loccitFFT
})

