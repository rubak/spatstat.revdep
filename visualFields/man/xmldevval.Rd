\name{xmldevval}
\alias{xmldevval}
\title{extracts total-deviation values, pattern-deviation values, total-deviation probability values, and pattern-deviation probability values}
\description{extracts total-deviation values, pattern-deviation values, total-deviation probability values, and pattern-deviation probability values}
\usage{xmldevval( xmllines, patternMap, typeData = c( "td" ),
           group = c( 4, 3, 2, 1, 0 ),
           cutoffs = c( 0.5, 1, 2, 5, 95 ) )}
\arguments{
  \item{xmllines}{lines from loaded XML files}
  \item{patternMap}{pattern of stimulus locations. Default is \code{saplocmap$p24d2}}
  \item{typeData}{Type of data to load; visual field (\code{vf}), total deviations (\code{td}), pattern deviations (\code{pd}), global indices (\code{gi}), visual-field index (\code{vfi}), total-deviation p-values (\code{tdp}), pattern-deviation p-values (\code{pdp}), global indices probability maps (\code{gip}), visual-field-index probability map (\code{vfip}). Default is \code{vf}}
  \item{group}{for probability maps: the probability group coding}
  \item{cutoffs}{for probability maps: the corresponding pvalue for each group code}
}

\value{total-deviation values, pattern-deviation values, total-deviation probability values, and pattern-deviation probability values}

\seealso{ \code{\link{loadvfxml}}, \code{\link{xmlvfval}} }

\author{Ivan Marin-Franch}