\name{plot.stim}
\alias{plot.stim}
\title{Plotting 'stim' objects}
\usage{
\method{plot}{stim}(x, fix.range = FALSE, sleep = 0.2, override.par = TRUE, ...)
}
\arguments{
\item{x}{An object of class \code{\link{stim}}.}
\item{fix.range}{Logical value indicating whether use the same color scale limits for each plot in the sequence. Ignored if the user supplies a pre-defined \code{\link[spatstat.geom]{colourmap}} to the \code{col} argument, which is matched to \code{...} and passed to \code{\link[spatstat.geom]{plot.im}}.}
\item{sleep}{Single positive numeric value giving the amount of time (in
seconds) to \code{\link[base]{Sys.sleep}} before drawing the next image in
the animation.}
\item{override.par}{Logical value indicating whether to override the
existing graphics device parameters prior to plotting, resetting
\code{mfrow} and \code{mar}. See `Details' for when you might want to
disable this.}
\item{...}{Additional graphical parameters to be passed to
\code{\link[spatstat.geom]{plot.im}} (see `Details').}
}

\value{
Plots to the active graphics device.
}

\description{
\code{plot} method for class \code{\link{stim}}.
}

\details{
Actual visualisation is deferred to
\code{\link[spatstat.geom]{plot.im}}, for which there are a variety of
customisations available the user can access via \code{...}.

The \code{\link{stim}} object is plotted as an animation, one pixel image
after another, separated by \code{sleep} seconds. If instead you intend the
individual images to be plotted in an array of images, you should first set
up your plot device layout, and ensure \code{override.par = FALSE} so that
the function does not reset these device parameters itself. In such an
instance, one might also want to set \code{sleep = 0}.
}
\examples{
# See help(stgmix) and help(stkey) for examples
}
\author{
T.M. Davies
}

