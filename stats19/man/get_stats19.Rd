% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/get.R
\name{get_stats19}
\alias{get_stats19}
\title{Download, read and format STATS19 data in one function.}
\usage{
get_stats19(
  year = NULL,
  type = "accidents",
  data_dir = get_data_directory(),
  file_name = NULL,
  format = TRUE,
  ask = FALSE,
  silent = FALSE,
  output_format = "tibble",
  ...
)
}
\arguments{
\item{year}{Single year for which file is to be downloaded.}

\item{type}{One of 'Accidents', 'Casualties', 'Vehicles'; defaults to 'Accidents'.
Or any variation of to search the file names with such as "acc" or "accid".}

\item{data_dir}{Parent directory for all downloaded files. Defaults to \code{tempdir()}.}

\item{file_name}{The file name (DfT named) to download.}

\item{format}{Switch to return raw read from file, default is \code{TRUE}.}

\item{ask}{Should you be asked whether or not to download the files? \code{TRUE} by default.}

\item{silent}{Boolean. If \code{FALSE} (default value), display useful progress
messages on the screen.}

\item{output_format}{A string that specifies the desired output format. The
default value is \code{"tibble"}. Other possible values are \code{"data.frame"}, \code{"sf"}
and \code{"ppp"}, that, respectively, returns objects of class \code{\link{data.frame}},
\code{\link[sf:sf]{sf::sf}} and \code{\link[spatstat.geom:ppp]{spatstat.geom::ppp}}. Any other string is ignored and a tibble
output is returned. See details and examples.}

\item{...}{Other arguments that should be passed to \code{\link[=format_sf]{format_sf()}} or
\code{\link[=format_ppp]{format_ppp()}} functions. Read and run the examples.}
}
\description{
Download, read and format STATS19 data in one function.
}
\section{Details}{

This function utilizes \code{dl_stats19} and \verb{read_*} functions and returns a
tibble (default), a data.frame, a sf object or a ppp object (according to the
\code{output_format} parameter).
The file downloaded would be for a specific year (e.g. 2017) or multiple
years (e.g. c(2017, 2018)).
See examples.

As this function uses \code{dl_stats19} function, it can download many MB of data,
so ensure you have a sufficient disk space.

If \code{output_format = "data.frame"} or \code{output_format = "sf"} or \code{output_format = "ppp"} then the output data is transformed into a data.frame, sf or ppp
object using the \code{\link[=as.data.frame]{as.data.frame()}} or \code{\link[=format_sf]{format_sf()}} or \code{\link[=format_ppp]{format_ppp()}}
functions, respectively.
See examples.
}

\examples{
\donttest{
# default tibble output
x = get_stats19(2019)
class(x)
x = get_stats19(2017, silent = TRUE)

# data.frame output
x = get_stats19(2019, silent = TRUE, output_format = "data.frame")
class(x)

# multiple years
get_stats19(c(2017, 2018), silent = TRUE)

# sf output
x_sf = get_stats19(2017, silent = TRUE, output_format = "sf")

# sf output with lonlat coordinates
x_sf = get_stats19(2017, silent = TRUE, output_format = "sf", lonlat = TRUE)
sf::st_crs(x_sf)

# multiple years
get_stats19(c(2017, 2018), silent = TRUE, output_format = "sf")

if (requireNamespace("spatstat.core", quietly = TRUE)) {
# ppp output
x_ppp = get_stats19(2017, silent = TRUE, output_format = "ppp")
spatstat.geom::plot.ppp(x_ppp, use.marks = FALSE)

# Multiple years
get_stats19(c(2017, 2018), silent = TRUE, output_format = "ppp")

# We can use the window parameter of format_ppp function to filter only the
# events occurred in a specific area. For example we can create a new bbox
# of 5km around the city center of Leeds

leeds_window = spatstat.geom::owin(
xrange = c(425046.1, 435046.1),
yrange = c(428577.2, 438577.2)
)

leeds_ppp = get_stats19(2017, silent = TRUE, output_format = "ppp", window = leeds_window)
spatstat.geom::plot.ppp(leeds_ppp, use.marks = FALSE, clipwin = leeds_window)

# or even more fancy examples where we subset all the events occurred in a
# pre-defined polygon area

# The following example requires osmdata package
# greater_london_sf_polygon = osmdata::getbb(
# "Greater London, UK",
# format_out = "sf_polygon"
# )
# spatstat works only with planar coordinates
# greater_london_sf_polygon = sf::st_transform(greater_london_sf_polygon, 27700)
# then we extract the coordinates and create the window object.
# greater_london_polygon = sf::st_coordinates(greater_london_sf_polygon)[, c(1, 2)]
# greater_london_window = spatstat.geom::owin(poly = greater_london_polygon)

# greater_london_ppp = get_stats19(2017, output_format = "ppp", window = greater_london_window)
# spatstat.geom::plot.ppp(greater_london_ppp, use.marks = FALSE, clipwin = greater_london_window)
}
}
}
\seealso{
\code{\link[=dl_stats19]{dl_stats19()}}

\code{\link[=read_accidents]{read_accidents()}}
}
