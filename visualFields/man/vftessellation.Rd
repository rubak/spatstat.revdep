\name{vftessellation}
\alias{vftessellation}
\title{Voronoi tessellation of vf spatial testing locations}
\description{calculates the Voronoi tessellation of vf spatial testing locations}
\usage{vftessellation( locmap, dist = 3 )}
\arguments{
  \item{locmap}{a vf-object}
  \item{dist}{This function calculates the outer convex hull of the Voronoi diagram. This parameter specifies how much in degrees we want to extend the outer convex hull. Default is 3 degrees}
}
\value{
returns the Voronoi diagrams and the outer convex hull of the visual field area tested.
}

\examples{
  vftess <- vftessellation( saplocmap$p24d2 )
}

\author{Ivan Marin-Franch}