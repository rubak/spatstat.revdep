\name{vfgrayscale}
\alias{vfgrayscale}
\title{gray scales for sensitivities}
\description{maps sensitivity values to grayscales}
\usage{vfgrayscale( sens, age, pattern, algorithm )}
\arguments{
  \item{sens}{an array with sensitivity values}
  \item{age}{age of the subject}
  \item{algorithm}{algorithm used in the visual test}
  \item{pattern}{pattern used in the visual test}
}
\details{This function maps sensitivity values to grayscales. The grayscale depends on age of the person, test pattern and algorithm}

\value{returns the RGB values for the gray scale}

\seealso{\code{\link{vfcolormap}}}

\author{Ivan Marin-Franch}