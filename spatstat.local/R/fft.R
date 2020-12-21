#
# fft.R
#
# Calculations for locppm that are possible using Fast Fourier Transform
#
# $Revision: 1.40 $ $Date: 2013/09/02 06:44:02 $

locppmFFT <- local({

  # Dependence between calculations
  #  ("Score", "Hessian" etc represent blocks of code in the main function)
  .Score   <- "Score"
  .Hessian <- c("Hdens", "Hessian")
  .Fish    <- c("Hdens", "Fish")
  .InvHess <- c(.Hessian, "InvHess")
  .Var     <- c(union(.InvHess, .Fish), "Var")
  .Coef    <- c(union(.Hessian, .Score), "Coef")
  .Tstat   <- c(union(.Var, .Coef), "Tstat")
  .Tgrad   <- c(union(.InvHess, .Coef), "Tgrad")

  # Table mapping user options to calculations required
  Dependence.Table <-
    list(coef    = .Coef,
         score   = .Score,
         var     = .Var,
         fish    = .Fish,
         tstat   = .Tstat,
         grad    = .Hessian,
         invgrad = .InvHess,
         tgrad   = .Tgrad)

  FullNames <- 
    c("coefficients",
      "score",
      "variance",
      "fisher",
      "tstatistic",
      "gradient",
      "invgradient",
      "tgradient")
  
  opties <- names(Dependence.Table)


  locppmFFT <- function(model, sigma, ...,
                        lambda=fitted(model, new.coef=new.coef),
                        new.coef=NULL,
                        what = "coef",
                        internals = NULL,
                        algorithm=c("density","closepairs", "chop"),
                        verbose=TRUE) {
  starttime <- proc.time()

  stopifnot(is.ppm(model))
  algorithm <- match.arg(algorithm)

  homcoef <- coef(model)
  use.coef <- if(!is.null(new.coef)) new.coef else homcoef

  if(missing(lambda) && !is.null(internals$lambda)) 
    lambda <- internals$lambda
  what <- match.arg(what, opties, several.ok=TRUE)
  give <- as.list(!is.na(match(opties, what)))
  names(give) <- opties

  needed <- unique(unlist(Dependence.Table[what]))

  p <- length(homcoef)
  Q <- quad.ppm(model)
  UQ <- union.quad(Q)
  X <- Q$data
  nX <- npoints(X)
  
  # make a template mask for images
  W <- as.mask(as.owin(UQ), ...)
  Wmat <- as.matrix(W)
  npixels <- length(Wmat)

  vec.names <- names(homcoef)
  mat.names <- as.vector(outer(vec.names, vec.names, paste, sep="."))

  if("Hdens" %in% needed) {
    # Evaluate integrand of Hessian at quadrature points
    if(is.null(m <- internals$mom))
      internals$mom <- m <- model.matrix(model)
    mo <- outersquare.flat.vector(m)
    hQ <- UQ %mark% (lambda * mo)
    # Interpolate onto pixel grid by nearest neighbour rule
    Hdens <- nnmark(hQ, at="pixels", ...)
    if(p == 1) Hdens <- list(Hdens)
    names(Hdens) <- mat.names
  } 
  
  if("Hessian" %in% needed) {
    # Compute negative local Hessian as list of images
    H <- lapply(Hdens, blur, sigma=sigma, bleed=FALSE, normalise=FALSE)
    H <- as.listof(H)
    # Convert to matrix
    Hmat <- imagelist2matrix(H)
  } else H <- NULL

  if("InvHess" %in% needed) {
    # Compute inverse of negative local Hessian as matrix
    invHmat <- invert.flat.matrix(Hmat, p)
    # convert back to images
    invH <- matrix2imagelist(invHmat, W)
    names(invH) <- mat.names
  } else invH <- NULL

  if("Score" %in% needed) {
    # Compute local score residuals as list of images
    R <- residuals(model, type="score")
    U <- Smooth(R, sigma=sigma, ..., edge=FALSE)
    if(is.im(U)) U <- list(U)
    U <- as.listof(U)
    names(U) <- vec.names
    # convert to matrix
    Umat <- imagelist2matrix(U)
  } else U <- NULL

  if("Coef" %in% needed) {
    # evaluate Taylor approximation of locally fitted coefficients
    if(!is.null(invH)) {
      # use inverse Hessian: delta(theta) = invH %*% U
      delta <- multiply2.flat.matrices(invHmat, Umat, c(p,p), c(p, 1))
    } else if(!is.null(H)) {
      # use Hessian: delta(theta) = solve(H, U)
      delta <- solve2.flat.matrices(Hmat, Umat, c(p, p), c(p, 1))
    } else stop("Internal error: coefficient calculation needs Hessian")
    Cmat0 <- matrix(homcoef, nrow=npixels, ncol=p, byrow=TRUE)
    Cmat <- Cmat0 + delta
    # convert back to images
    Coefs <- matrix2imagelist(Cmat, W)
    names(Coefs) <- vec.names
  } else Coefs <- NULL
  
  if("Fish" %in% needed) {
    # compute local Fisher information
    # Poisson/Poincare term: integral weighted by square of smoothing kernel
    A1 <- lapply(Hdens, blur, sigma=sigma/sqrt(2),
                   bleed=FALSE, normalise=FALSE)
    A1 <- lapply(A1, function(z, a) { eval.im(a * z) },
                 a = 1/(4 * pi * sigma^2))
    if(is.poisson(model)) {
      # Poisson process
      Fish <- A1
    } else {
      # Gibbs process
      # Compute cross-dependence terms A2, A3
      Z <- internals$Z %orifnull% is.data(Q)
      ok <- internals$ok %orifnull% getglmsubset(model)
      okX <- ok[Z]
      Xok <- X[okX]
      # W is template mask
      xcolW <- W$xcol
      yrowW <- W$yrow
      xrangeW <- W$xrange
      yrangeW <- W$yrange
      #
      # Sufficient statistic h:
      # momX[i, ] = h(X[i] | X)
      momX <- internals$mom[Z, ,drop=FALSE]
      # momdel[ ,i,j] = h(X[i] | X[-j])
      momdel <- internals$momdel
      # mom.array[ ,i,j] = momX[i, ] = h(X[i] | X)
      mom.array <- array(t(momX), dim=c(p, nX, nX))
      # ddS[ ,i,j] = h(X[i] | X) - h(X[i] | X[-j])
      ddS <- mom.array - momdel
      #
      # Conditional intensity lambda:
      lamX <- lambda[Z]
      # lamdel[i,j] = lambda(X[i] | X[-j])
      lamdel <- internals$lamdel
      if(is.null(lamdel))
        lamdel <- matrix(lamX, nX, nX) * exp(tensor::tensor(-use.coef, ddS, 1, 1))
      # pairwise weight for A2
      #   pairweight[i,j] = lamdel[i,j]/lambda[i] - 1
      pairweight <- lamdel/lamX - 1
      # 
      if(algorithm %in% c("closepairs", "chop")) {
        R <- reach(model)
        if(!is.finite(R)) {
          warning(paste("Reach of model is infinite or NA;",
                        "algorithm", dQuote(algorithm), "not used"))
          algorithm <- "density"
        }
      }
      # start
      A2 <- A3 <- rep(list(0), p^2)
      switch(algorithm,
             density = {
               # run through points X[i] that contributed to pseudolikelihood
               ii <- which(okX)
               nok <- length(ii)
               if(verbose)
                 cat(paste("Running through", nok,
                           "data points out of", nX, "..."))
               for(k in seq_len(nok)) {
                 if(verbose) progressreport(k, nok)
                 i <- ii[k]
                 # form a flat matrix dS_i[j] = outer(ddS[ ,i,j], ddS[,j,i])
                 dSi <- outer2.flat.vectors(t(ddS[,i,okX]), t(ddS[,okX,i]))
                 # form a flat matrix containing
                 #   hout_i[j] = pairweight[i,j] h(X[i]|X[-j]) h(X[j]|X[-i])^T 
                 hout.i <- pairweight[i,okX] *
                   outer2.flat.vectors(t(momdel[,i,okX]), t(momdel[,okX,i]))
                 # associate dS_i[j] to X[j] and convolve with kernel
                 Si <- lapply(as.data.frame(dSi),
                              function(z, XX, sigma, W) {
                                density(XX, sigma=sigma, weights=z,
                                        W=W, edge=FALSE)
                              }, XX=Xok, sigma=sigma, W=W)
                 # associate hout_i[j] to X[j] and convolve with kernel
                 h.i <- lapply(as.data.frame(hout.i),
                               function(z, XX, sigma, W) {
                                 density(XX, sigma=sigma, weights=z,
                                         W=W, edge=FALSE)
                               }, XX=Xok, sigma=sigma, W=W)
                 # kernel at X[i]
                 #   Di <- density(X[i], sigma=sigma, W=W, edge=FALSE)
                 dx <- dnorm(xcolW - X$x[i], sd=sigma)
                 dy <- dnorm(yrowW - X$y[i], sd=sigma)
                 Di <- im(outer(dy, dx, "*"),
                          xcol=xcolW, yrow=yrowW,
                          xrange=xrangeW, yrange=yrangeW)
                 # multiply pointwise
                 DiSi <- lapply(Si, function(z, y) eval.im(z * y), y=Di)
                 Dihi <- lapply(h.i, function(z, y) eval.im(z * y), y=Di)
                 # increment sums
                 A2 <- Map(function(a,b) eval.im(a+b), A2, Dihi)
                 A3 <- Map(function(a,b) eval.im(a+b), A3, DiSi)
               }
               Fish <- Map(function(x,y,z) eval.im(x+y+z), A1, A2, A3)
             },
             closepairs = {
               # identify close pairs of data points
               cl <- closepairs(X, R, what="indices")
               I <- cl$i
               J <- cl$j
               # restrict to pairs of points which both contributed to PL
               okIJ <- okX[I] & okX[J]
               I <- I[okIJ]
               J <- J[okIJ]
               npairs <- length(I)
               # group by index I
               JsplitI <- split(J, I)
               Igroup <- as.integer(names(JsplitI))
               ngroups <- length(JsplitI)
               # compute kernel centred at each data point
               kernels <- vector(mode="list", length=nX)
               for(k in sort(unique(c(I,J))))
                 kernels[[k]] <-
                   as.vector(outer(dnorm(yrowW - X$y[k], sd=sigma),
                                   dnorm(xcolW - X$x[k], sd=sigma)))
               # sum over close pairs
               arrdim <- c(p^2, length(yrowW), length(xcolW))
               A2array <- A3array <- numeric(prod(arrdim))
               #
               if(verbose)
                 cat(paste("Running through", ngroups,
                           "close groups of points..."))
               for(k in seq_len(ngroups)) {
                 if(verbose) progressreport(k, ngroups)
                 i <- Igroup[k]
                 ki <- kernels[[i]]
                 for(j in JsplitI[[k]]) {
                   kj <- kernels[[j]]
                   kikj <- ki * kj
                   A3ij <- outer(ddS[,i,j], ddS[,j,i])
                   A2ij <- pairweight[i,j] * outer(momdel[,i,j],momdel[,j,i])
                   A2array <- A2array + as.vector(outer(as.vector(A2ij), kikj))
                   A3array <- A3array + as.vector(outer(as.vector(A3ij), kikj))
                 }
               }
               A2array <- array(A2array, dim=arrdim)
               A3array <- array(A3array, dim=arrdim)
               # convert to images
               for(k in seq_len(p^2)) {
                 A2[[k]] <- im(A2array[k,,],
                               xcol=xcolW, yrow=yrowW,
                               xrange=xrangeW, yrange=yrangeW)
                 A3[[k]] <- im(A3array[k,,],
                               xcol=xcolW, yrow=yrowW,
                               xrange=xrangeW, yrange=yrangeW)
               }
               # Finally evaluate variance
               Fish <- Map(function(x,y,z) eval.im(x+y+z), A1, A2, A3)
             },
             chop = {
               # identify close pairs of data points
               cl <- closepairs(X, R, what="indices")
               I <- cl$i
               J <- cl$j
               # restrict to pairs of points which both contributed to PL
               okIJ <- okX[I] & okX[J]
               I <- I[okIJ]
               J <- J[okIJ]
               npairs <- length(I)
               # compute Gaussian kernel on larger image raster
               nxc <- length(xcolW)
               nyr <- length(yrowW)
               rx <- W$xrange
               ry <- W$yrange
               xstep <- W$xstep
               ystep <- W$ystep
               xx <- seq(-nxc * xstep, nxc * xstep, length = 2 * nxc + 1)
               yy <- seq(-nyr * ystep, nyr * ystep, length = 2 * nyr + 1)
               Gauss <- outer(dnorm(yy, sd=sigma),
                              dnorm(xx, sd=sigma),
                              "*")
               # sum over close pairs
               arrdim <- c(p^2, length(yrowW), length(xcolW))
               A2array <- A3array <- numeric(prod(arrdim))
               if(verbose)
                 cat(paste("Running through", npairs,
                           "close pairs of points..."))
               for(k in seq_len(npairs)) {
                 if(verbose) progressreport(k, npairs)
                 i <- I[k]
                 j <- J[k]
                 # kernel centred at X[i]
                 rowXi <- grid1index(X$y[i], ry, nyr)
                 colXi <- grid1index(X$x[i], rx, nxc)
                 ki <- Gauss[nyr - rowXi + 1:nyr,
                             nxc - colXi + 1:nxc]
                 ki <- as.vector(ki)
                 # kernel centred at X[j]
                 rowXj <- grid1index(X$y[j], ry, nyr)
                 colXj <- grid1index(X$x[j], rx, nxc)
                 kj <- Gauss[nyr - rowXj + 1:nyr,
                             nxc - colXj + 1:nxc]
                 kj <- as.vector(kj)
                 # multiply
                 kikj <- ki * kj
                 A3ij <- outer(ddS[,i,j], ddS[,j,i], "*")
                 A2ij <- pairweight[i,j] * outer(momdel[,i,j],momdel[,j,i], "*")
                 A2array <- A2array + as.vector(outer(as.vector(A2ij), kikj))
                 A3array <- A3array + as.vector(outer(as.vector(A3ij), kikj))
               }
               A2array <- array(A2array, dim=arrdim)
               A3array <- array(A3array, dim=arrdim)
               # convert to images
               for(k in seq_len(p^2)) {
                 A2[[k]] <- im(A2array[k,,],
                               xcol=xcolW, yrow=yrowW,
                               xrange=xrangeW, yrange=yrangeW)
                 A3[[k]] <- im(A3array[k,,],
                               xcol=xcolW, yrow=yrowW,
                               xrange=xrangeW, yrange=yrangeW)
               }
               # Finally evaluate variance
               Fish <- Map(function(x,y,z) eval.im(x+y+z), A1, A2, A3)
             })
    }
    names(Fish) <- mat.names
    Fmat <- imagelist2matrix(Fish)
  } else Fish <- NULL

  if("Var" %in% needed) {
    # evaluate local variance (assuming Poisson)
    #  v = H^{-1} %*% I %*% H^{-1}
    Vmat <- quadform2.flat.matrices(Fmat, invHmat, c(p,p), c(p,p))
    # convert back to images
    Vars <- matrix2imagelist(Vmat, W)
    names(Vars) <- mat.names
  } else Vars <- NULL

  if("Tstat" %in% needed) {
    # evaluate t statistics of Taylor approximate fit
    # Extract standard errors
    diagindices <- diag(matrix(1:(p^2), p, p))
    Vdiags <- Vmat[, diagindices, drop=FALSE]
    # Standardise 
    Tmat <- Cmat/sqrt(Vdiags)
    # convert back to images
    Tstat <- matrix2imagelist(Tmat, W)
    names(Tstat) <- vec.names
  } else Tstat <- NULL

  if("Tgrad" %in% needed) {
    # evaluate surrogate t statistic based on inverse Hessian
    diagindices <- diag(matrix(1:(p^2), p, p))
    invHdiags <- invHmat[, diagindices, drop=FALSE]
    TGmat <- Cmat/sqrt(invHdiags)
    # convert back to images
    Tgrad <- matrix2imagelist(TGmat, W)
    names(Tgrad) <- vec.names
  } else Tgrad <- NULL

  answer <- list(coefficients = Coefs,
                 score        = U,
                 variance     = Vars,
                 fisher       = Fish,
                 tstatistic   = Tstat,
                 gradient     = H,
                 invgradient  = invH,
                 tgradient    = Tgrad)

  # return only those quantities that were commanded
  indx <- pmatch(what, opties) # there are no NA's
  answer <- answer[FullNames[indx]]
  
  return(timed(answer, starttime=starttime))
}

  locppmFFT
})


  matrix2imagelist <- function(mat, W) {
    # 'mat' is a matrix whose rows correspond to pixels
    # and whose columns are different quantities e.g. components of score.
    # Convert it to a list of images using the template 'W'
    if(!is.matrix(mat)) mat <- matrix(mat, ncol=1)
    Mats <- lapply(as.list(as.data.frame(mat)),
                   matrix,
                   nrow=W$dim[1], ncol=W$dim[2])
    result <- as.listof(lapply(Mats, as.im, W=W))
    return(result)
  }

  imagelist2matrix <- function(x) {
    if(is.im(x)) x <- list(x)
    y <- as.matrix(as.data.frame(lapply(x, 
                                        function(z) as.vector(as.matrix(z)))))
    colnames(y) <- names(x)
    return(y)
  }

