\name{spatstat.local-internal}
\title{Internal Functions of the spatstat.local Package}
\alias{print.loccit}
\alias{locppmEngine}
\alias{locppmFFT}
\alias{vcovlocEngine}
\alias{getvcinternals}
\alias{sqmag}
\alias{sqrtmat}
\alias{invsqrtmat}
\alias{gridproxy}
\alias{locppmOptions}
\alias{print.locppmOptions}
\alias{FirstExtantEntry}
\alias{sample.imagelist}
\alias{locppmPredict}
\alias{check.flat.matrix}
\alias{as.flat.matrix}
\alias{handle.flat.matrix}
\alias{invert.flat.matrix}
\alias{trace.flat.matrix}
\alias{transpose.flat.matrix}
\alias{eigenvalues.flat.matrix}
\alias{average.flat.matrix}
\alias{invert.slice}
\alias{outersquare.flat.vector}
\alias{check2.flat.matrices}
\alias{handle2.flat.matrices}
\alias{multiply2.flat.matrices}
\alias{multiply2.slice}
\alias{solve2.flat.matrices}
\alias{solve2.slice}
\alias{outer2.flat.vectors}
\alias{quadform2.flat.matrices}
\alias{quadform2.slice}
\alias{quadform2symm.slice}
\alias{check3.flat.matrices}
\alias{handle3.flat.matrices}
\alias{bilinear3.slice}
\alias{bilinear3.flat.matrices}
\alias{multiply3.slice}
\alias{multiply3.flat.matrices}
\alias{matrix2imagelist}
\alias{imagelist2matrix}
\alias{loccitFFT}
\alias{extraClusterModelInfo}
\alias{applymaps}
\alias{accumulateStatusList}
\alias{HomTestMapEngine}
\alias{print.homtestmap}
\description{
  Internal functions of the \pkg{spatstat.local} package.
}
\usage{
 sqmag(x)
 invsqrtmat(M)
 sqrtmat(M)
 gridproxy(P, \dots, dimyx = NULL, eps = NULL, xy = NULL, weights=NULL)
 \method{print}{loccit}(x, \dots)
 locppmEngine(model, sigma, V, \dots, weights,
               verbose, scopename, scopeindex,
               opt, dropterm, matrices, fastRCinloop, internals,
               fft.algorithm)
 locppmFFT(model, sigma, \dots, lambda, new.coef, what, internals,
  algorithm, verbose)
 vcovlocEngine(internals, localwt=NULL, \dots,
               A1dummy=FALSE, new.coef=NULL, bananas=FALSE)
 getvcinternals(model, verbose)
 locppmOptions(other, \dots, other1)
 \method{print}{locppmOptions}(x, \dots)
 FirstExtantEntry(xlist, tags, whinge)
 sample.imagelist(X,V)
 locppmPredict(homfit, coefs, type, locations, precomputed, details, index)
 check.flat.matrix(A, dimA)
 as.flat.matrix(X, ncopies)
 handle.flat.matrix(A, dimA, matrixop, \dots)
 invert.flat.matrix(X,p,special)
 trace.flat.matrix(X,p)
 transpose.flat.matrix(X,dimX)
 eigenvalues.flat.matrix(X,p)
 average.flat.matrix(X,dimX,weights)
 invert.slice(x,p)
 outersquare.flat.vector(A)
 check2.flat.matrices(A, B, dimA, dimB)
 handle2.flat.matrices(A, B, dimA, dimB, operation)
 multiply2.flat.matrices(A, B, dimA, dimB)
 multiply2.slice(x, dimA, dimB, indA, indB)
 solve2.flat.matrices(A, B, dimA, dimB)
 solve2.slice(x, dimA, dimB, indA, indB)
 outer2.flat.vectors(A,B)
 quadform2.flat.matrices(A, B, dimA, dimB, Bsymmetric)
 quadform2.slice(x, dimA, dimB, indA, indB)
 quadform2symm.slice(x, dimA, dimB, indA, indB)
 check3.flat.matrices(X, Y, Z, dimX, dimY, dimZ)
 handle3.flat.matrices(X, Y, Z, dimX, dimY, dimZ, operation)
 bilinear3.slice(x, dimX, dimY, dimZ, indX, indY, indZ)
 bilinear3.flat.matrices(X, Y, Z, dimX, dimY, dimZ)
 multiply3.slice(x, dimX, dimY, dimZ, indX, indY, indZ)
 multiply3.flat.matrices(X, Y, Z, dimX, dimY, dimZ)
 matrix2imagelist(mat, W)
 imagelist2matrix(x)
 loccitFFT(fit0, sigma, rmax, \dots,
           base.trendcoef, base.cluspar, base.lambda, base.lambdaim,
           clusters, hompoisfit,
           what, do.trend, do.clusters, calcopt, verbose)
 extraClusterModelInfo(name, warn)
 applymaps(maplist, x)
 accumulateStatusList(x, stats, stoponerror)
 HomTestMapEngine(x, info)
 \method{print}{homtestmap}(x, \dots)
}
\details{
  These internal functions should not be called directly by the user.
  They may change at any time.
}
\keyword{internal}
