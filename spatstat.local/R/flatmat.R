#
# flatmat.R
#
# Code for performing linear algebra operations "in parallel"
# on each of a list of matrices/vectors, e.g.
#     answer[[i]] <- A[[i]] %*% B[[i]]
# where dim(A[[i]]) is the same for all i,
# and dim(B[[i]]) is the same for all i. 
#
# Instead of storing the arguments as lists of matrices,
# we flatten each matrix A[[i]] into a vector,
# and store them as the rows of a single gigantic matrix,
# which can then be processed in parallel.
# 
# In the spatstat package, for example,
# this code is used for handling 'vector-valued'
# and 'matrix-valued' spatial covariates.
#
#  $Revision: 1.14 $  $Date: 2016/10/02 02:06:32 $
#

check.flat.matrix <- function(A, dimA) {
  # 'A' represents a list of matrices of dimension 'dimA'
  Aname <- short.deparse(substitute(A))
  dimAname <- short.deparse(substitute(dimA))
  if(!is.matrix(A))
    stop(paste(Aname, "should be a matrix"))
  if(!is.vector(dimA) || length(dimA) != 2)
    stop(paste(dimAname, "should be a vector of length 2"))
  if(ncol(A) != prod(dimA))
    stop(paste("Incorrect dimensions supplied for", Aname))
  return(invisible(NULL))
}

as.flat.matrix <- function(X, ncopies) {
  stopifnot(is.matrix(X))
  Y <- matrix(as.vector(X), ncol=prod(dim(X)), nrow=ncopies, byrow=TRUE)
  dn <- as.list(dimnames(X))
  if(is.null(dn[[1]])) dn[[1]] <- 1:nrow(X)
  if(is.null(dn[[2]])) dn[[2]] <- 1:ncol(X)
  colnames(Y) <- outer(dn[[1]], dn[[2]], paste, sep=".")
  return(Y)
}

# general operation on one matrix

handle.flat.matrix <- local({

  handle.flat.matrix <- function(A, dimA, matrixop, ...) {
    ## apply some matrix operation to a stack of flat matrices
    y <- apply(A, 1, handleone, 
               nr=dimA[1], nc=dimA[2], matrixop=matrixop, ...)
    y <- if(is.matrix(y)) t(y) else matrix(y, ncol=1)
    return(y)
  }

  handleone <- function(z, nr, nc, matrixop, ...) {
    x <- matrix(z, nr, nc)
    y <- matrixop(x, ...)
    return(y)
  }

  handle.flat.matrix
})
  
eigenvalues.flat.matrix <- local({

  eigenvalues.flat.matrix <- function(X, p) {
    check.flat.matrix(X, c(p,p))
    handle.flat.matrix(X, c(p,p), eigenvals)
  }

  eigenvals <- function(M) {
    eigen(M, symmetric=TRUE, only.values=TRUE)$values
  }

  eigenvalues.flat.matrix
})

invert.flat.matrix <- function(X, p, special=TRUE) {
  # X is a matrix whose rows can be interpreted as p * p matrices.
  # Invert them.
  check.flat.matrix(X, c(p, p))
  if(special && p == 1) {
    # scalar
    Y <- 1/X
  } else if(special && p == 2) {
    # 2 * 2 matrix, by hand
    aa <- X[,1]
    bb <- X[,3]
    cc <- X[,2]
    dd <- X[,4]
    dets <- aa * dd - bb * cc
    invdet <- ifelse(dets == 0, NA, 1/dets)
    Y <- invdet * cbind(dd, -cc, -bb, aa)
  } else {
    # use general algorithm
    Y <- apply(X, 1, invert.slice, p=p)
    Y <- if(p != 1) t(Y) else matrix(Y, ncol=1)
  }
  colnames(Y) <- colnames(X)
  return(Y)
}

trace.flat.matrix <- function(X, p) {
  # X is a matrix whose rows can be interpreted as p * p matrices.
  # Extract the traces of the matrices.
  check.flat.matrix(X, c(p,p))
  ind <- diag(matrix(1:(p^2), p, p))
  Y <- if(p == 1) as.vector(X) else rowSums(X[,ind])
  return(Y)
}

transpose.flat.matrix <- function(X, dimX) {
  check.flat.matrix(X, dimX)
  indices <- matrix(1:prod(dimX), dimX[1], dimX[2])
  indices <- as.vector(t(indices))
  X[, indices, drop=FALSE]
}

average.flat.matrix <- function(X, dimX, weights=NULL) {

  if(is.null(weights)) {
    Xbar <- colMeans(X, na.rm=TRUE)
  } else {
    check.nvector(weights, nrow(X), things="matrices", oneok=TRUE)
    if(length(weights) == 1) weights <- rep(weights, nrow(X))
    Xbar <- apply(X, 2, weighted.mean, w=weights, na.rm=TRUE)
  }
  Xbar <- matrix(Xbar, dimX[1], dimX[2])
  return(Xbar)
}

# functions named '*.slice' act on an individual row.
# They unpack the matrix, and perform the desired operation.

invert.slice <- function(x, p) {
  mat <- matrix(x, p, p)
  y <- try(solve(mat), silent=TRUE)
  if(!inherits(y, "try-error")) return(y)
  return(matrix(NA, p, p))
}

# compute outer(A[[i]], A[[i]])
outersquare.flat.vector <- function(A) {
  nc <- ncol(A)
  if(nc == 1) return(A^2)
  if(nc == 2) {
    A1A1 <- A[,1]^2
    A1A2 <- A[,1] * A[,2]
    A2A2 <- A[,2]^2
    return(cbind(A1A1, A1A2, A1A2, A2A2))
  }
  ans <- apply(A, 1, function(z) { outer(z, z, "*") })
  return(t(ans))
}

# ................... two matrices ..........................

check2.flat.matrices <- function(A, B, dimA, dimB) {
  # 'A', 'B' represent lists of matrices of dimension 'dimA', 'dimB' resp
  check.flat.matrix(A, dimA)
  check.flat.matrix(B, dimB)
  # check that A and B are lists of the same length
  stopifnot(nrow(A) == nrow(B))
  return(invisible(NULL))
}

handle2.flat.matrices <- function(A, B, dimA, dimB, operation) {
  # apply some operation to a pair of compatible stacks of matrices
  lA <- ncol(A)
  lB <- ncol(B)
  z <- apply(cbind(A,B),
             1,
             operation,
             dimA = dimA, dimB=dimB,
             indA = 1:lA, indB = lA + (1:lB))
  z <- if(is.vector(z)) matrix(z, ncol=1) else t(z)
  return(z)
}

multiply2.flat.matrices <- function(A, B, dimA, dimB) {
  # multiply two stacks of matrices
  check2.flat.matrices(A, B, dimA, dimB)
  z <- if(prod(dimA) == 1 && prod(dimB) == 1) {
    # scalars
    A * B
  } else if(dimA[1] == 1 && dimB[2] == 1) {
    # row vector * column vector = scalar
    matrix(rowSums(A * B), ncol=1)
  } else {
    handle2.flat.matrices(A, B, dimA, dimB, multiply2.slice)
  }
  return(z)
}

multiply2.slice <- function(x, dimA, dimB, indA, indB) {
  A <- matrix(x[indA], dimA[1], dimA[2])
  B <- matrix(x[indB], dimB[1], dimB[2])
  A %*% B
}

solve2.flat.matrices <- function(A, B, dimA, dimB) {
  # solve(A,b) = A^{-1} b
  check2.flat.matrices(A, B, dimA, dimB)
  if(dimA[1] != dimA[2])
    stop("The dimensions of A should be square")
  if(dimA[2] != dimB[1])
    stop("Incompatible matrix dimensions for solve(A,B)")
  z <- if(prod(dimA) == 1) {
    # A is a scalar
    B/A
  } else {
    handle2.flat.matrices(A, B, dimA, dimB, solve2.slice)
  }
  return(z)
}

solve2.slice <- function(x, dimA, dimB, indA, indB) {
  A <- matrix(x[indA], dimA[1], dimA[2])
  B <- matrix(x[indB], dimB[1], dimB[2])
  y <- try(solve(A, B), silent=TRUE)
  if(!inherits(y, "try-error")) return(y)
  return(matrix(NA, dimB[1], dimB[2]))
}

outer2.flat.vectors <- function(A, B) {
  # compute outer(A[[i]], B[[i]])
  stopifnot(identical(dim(A), dim(B)))
  nc <- ncol(A)
  if(nc == 1) return(A * B)
  AB <- apply(cbind(A,B), 1,
               function(z, nc) { outer(z[1:nc], z[nc + (1:nc)]) },
               nc=nc)
  return(t(AB))
}

quadform2.flat.matrices <- function(A, B, dimA, dimB, Bsymmetric=FALSE) {
  # compute the quadratic form B[[i]] %*% A[[i]] %*% t(B[[i]])
  check2.flat.matrices(A, B, dimA, dimB)
  z <- if(prod(dimA) == 1 && prod(dimB) == 1) {
    # scalars
    A * B^2
  } else if(Bsymmetric) {
    handle2.flat.matrices(A, B, dimA, dimB, quadform2symm.slice)
  } else {
    handle2.flat.matrices(A, B, dimA, dimB, quadform2.slice)
  }
  return(z)
}

quadform2.slice <- function(x, dimA, dimB, indA, indB) {
  A <- matrix(x[indA], dimA[1], dimA[2])
  B <- matrix(x[indB], dimB[1], dimB[2])
  return(B %*% A %*% t(B))
}
  
quadform2symm.slice <- function(x, dimA, dimB, indA, indB) {
  A <- matrix(x[indA], dimA[1], dimA[2])
  B <- matrix(x[indB], dimB[1], dimB[2])
  return(B %*% A %*% B)
}
  
# ................... three matrices ..........................

check3.flat.matrices <- function(X, Y, Z, dimX, dimY, dimZ) {
  check.flat.matrix(X, dimX)
  check.flat.matrix(Y, dimY)
  check.flat.matrix(Z, dimZ)
  stopifnot(nrow(X) == nrow(Y))
  stopifnot(nrow(X) == nrow(Z))
  return(invisible(NULL))
}

handle3.flat.matrices <- function(X, Y, Z, dimX, dimY, dimZ, operation) {
  lX <- ncol(X)
  lY <- ncol(Y)
  lZ <- ncol(Z)
  ans <- apply(cbind(X,Y,Z),
               1,
               operation,
               dimX = dimX, dimY=dimY, dimZ=dimZ,
               indX = 1:lX, indY = lX + (1:lY), indZ=lX + lY + (1:lZ))
  ans <- if(is.vector(ans)) matrix(ans, ncol=1) else t(ans)
  return(ans)
}

bilinear3.slice <- function(x, dimX, dimY, dimZ, indX, indY, indZ) {
  X <- matrix(x[indX], dimX[1], dimX[2])
  Y <- matrix(x[indY], dimY[1], dimY[2])
  Z <- matrix(x[indZ], dimZ[1], dimZ[2])
  return(X %*% Y %*% t(Z))
}

bilinear3.flat.matrices <- function(X, Y, Z, dimX, dimY, dimZ) {
  # compute the quadratic form X[[i]] %*% Y[[i]] %*% t(Z[[i]])
  check3.flat.matrices(X, Y, Z, dimX, dimY, dimZ)
  z <- if(prod(dimX) == 1 && prod(dimY) == 1 && prod(dimZ) == 1) {
    # scalars
    X * Y * Z
  } else {
    handle3.flat.matrices(X, Y, Z, dimX, dimY, dimZ, bilinear3.slice)
  }
  return(z)
}

multiply3.slice <- function(x, dimX, dimY, dimZ, indX, indY, indZ) {
  X <- matrix(x[indX], dimX[1], dimX[2])
  Y <- matrix(x[indY], dimY[1], dimY[2])
  Z <- matrix(x[indZ], dimZ[1], dimZ[2])
  return(X %*% Y %*% Z)
}

multiply3.flat.matrices <- function(X, Y, Z, dimX, dimY, dimZ) {
  # compute X[[i]] %*% Y[[i]] %*% Z[[i]]
  check3.flat.matrices(X, Y, Z, dimX, dimY, dimZ)
  z <- if(prod(dimX) == 1 && prod(dimY) == 1 && prod(dimZ) == 1) {
    # scalars
    X * Y * Z
  } else {
    handle3.flat.matrices(X, Y, Z, dimX, dimY, dimZ, multiply3.slice)
  }
  return(z)
}
