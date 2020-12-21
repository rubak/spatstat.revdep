\name{exposurepps}
\alias{exposurepps}

\title{Example datasets for package spatialsegregation}
\description{Example datasets with 9 different scattering+exposure combinations.}
\usage{
   data(exposurepps)
}
\format{
A list with nine elements of class \code{ppp}.
}
\details{
A list of 9 point patterns with different degrees of intra-species clustering and 
inter-species mingling.

The patterns are synthetically produced using a combination of the functions of this
package as an energy function of a Gibbs model. 
}


\author{
Tuomas Rajala
tuomas.rajala@iki.fi
}

\examples{
 data(exposurepps)
 par(mfrow=c(3,3), mar=c(2,2,2,2))
 for(i in 1:9)plot(exposurepps[[i]])

# upper row has strong inter-type mingling effect
# lower row has strong inter-type repulsion or segregation
# left column has strong intra-type clustering
# right column has strong intra-type repulsion.

 par(mfrow=c(3,3), mar=c(3,3,4,3))
 for(i in 1:9)plot(  isarF(exposurepps[[i]]), cbind(ISARmean,theo)~par)

}

\keyword{ spatial }
