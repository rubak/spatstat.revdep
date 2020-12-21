\name{globalKinhom-package}
\alias{globalKinhom-package}
\alias{globalKinhom}
\docType{package}
\encoding{UTF-8}
\title{
\packageTitle{globalKinhom}
}
\description{
\packageDescription{globalKinhom}
}
\details{

The DESCRIPTION file:
\packageDESCRIPTION{globalKinhom}
\packageIndices{globalKinhom}

This package accompanies Shaw et al (2020). It provides \dQuote{global}
estimators for the non-parametric \eqn{K}- and pair correlation functions, which
summarize the second order interactions of second-order intensity-reweighted
stationary (SOIRS) point processes. These estimators provide an alternative to
those proposed by Baddeley et al (2000) for SOIRS point processes, which we
refer to as \dQuote{local} estimators. The local estimators are implemented in
the spatstat package as \code{\link{pcfinhom}} and \code{\link{Kinhom}}, with
\code{\link{pcfcross.inhom}} and \code{\link{Kcross.inhom}} for the
corresponding cross-pcf and cross-\eqn{K}-function.

Where possible, the interfaces are made to match those used by the \link{spatstat} package.
}
\author{
\packageAuthor{globalKinhom}

Maintainer: \packageMaintainer{globalKinhom}
}
\references{
    T Shaw, J Møller, R Waagepetersen. 2020. \dQuote{Globally Intensity-Reweighted Estimators for
    \eqn{K}- and pair correlation functions}. arXiv:2004.00527 [stat.ME].

    A Baddeley, J Møller, R Waagepetersen. 2000. \dQuote{Non- and Semi-Parametric Estimation of
    Interaction in Inhomogeneous Point Patterns}. Statistica Neerlandica 54, 329-350.
}
\keyword{ package }
\seealso{ \link{spatstat}, \code{\link{Kglobal}}, \code{link{pcfglobal}}
}
