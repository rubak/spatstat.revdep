\name{vf2gcloc}
\alias{vf2gcloc}
\title{Calculates the location of the GC soma corresponding to vf locations}
\description{Calculates the corresponding location of the GC soma for a list of vf locations. See [1,2]}
\usage{vf2gcloc( xy )}
\arguments{
  \item{xy}{data frame with the xy VF locations in degrees of visual angle}
}

\value{returns GC soma position in degrees of visual angle}

\examples{
  vf2gcloc( data.frame( x = 1, y = 1 ) )
}

\references{
[1] D. C. Hood, A. S. Raza, D. M. C. G. V., J. G. Odel, V. C. Greenstein,
J. M. Liebmann, and R. Ritch. Initial arcuate defects within the central 10
degrees in glaucoma. Investigative Ophthalmology and Visual Science,
52(2):940-946, 2011.

[2] A. S. Raza, J. Cho, D. M. C. G. V., H. Wang, X. Zhang, R. H. Kardon,
J. M. Liebmann, R. Ritch, and D. C. Hood. Retinal ganglion cell layer
thickness and local visual field sensitivity in glaucoma. Archives of
Ophthalmology, 129(12):1529-1536, 2011.
}

\author{Ivan Marin-Franch}