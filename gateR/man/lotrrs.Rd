% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/lotrrs.R
\name{lotrrs}
\alias{lotrrs}
\title{A single gate for two conditions}
\usage{
lotrrs(
  dat,
  alpha = 0.05,
  p_correct = "none",
  nbc = NULL,
  doplot = FALSE,
  rcols = c("#FF0000", "#cccccc", "#0000FF"),
  lower_lrr = NULL,
  upper_lrr = NULL,
  c1n = NULL,
  c2n = NULL,
  win = NULL,
  verbose = FALSE,
  ...
)
}
\arguments{
\item{dat}{Input data frame flow cytometry data with five (5) features (columns): 1) ID, 2) Condition A ID, 3) Condition B ID, 4) Marker A as x-coordinate, 5) Marker B as y-coordinate.}

\item{alpha}{Numeric. The two-tailed alpha level for significance threshold (default is 0.05).}

\item{p_correct}{Character string specifying whether to apply a correction for multiple comparisons including a Bonferroni correction \code{p_correct = "uncorrelated"} or a correlated Bonferroni correction \code{p_correct = "correlated"}. If \code{p_correct = "none"} then no correction is applied.}

\item{nbc}{Optional. An integer for the number of bins when \code{p_correct = "correlated"}. Similar to \code{nbclass} argument in \code{\link[pgirmess]{correlog}}. The default is the average number of gridded knots in one-dimension (i.e., x-axis).}

\item{doplot}{Logical. If \code{TRUE}, the output includes basic data visualizations.}

\item{rcols}{Character string of length three (3) specifying the colors for: 1) group A (numerator), 2) neither, and 3) group B (denominator) designations. The defaults are \code{c("#FF0000", "#cccccc", "#0000FF")} or \code{c("red", "grey80", "blue")}.}

\item{lower_lrr}{Optional, numeric. Lower cut-off value for the log relative risk value in the color key (typically a negative value). The default is no limit and the color key will include the minimum value of the log relative risk surface.}

\item{upper_lrr}{Optional, numeric. Upper cut-off value for the log relative risk value in the color key (typically a positive value). The default is no limit and the color key will include the maximum value of the log relative risk surface.}

\item{c1n}{Optional, character. The name of the level for the numerator of condition A. The default is null and the first level is treated as the numerator.}

\item{c2n}{Optional, character. The name of the level for the numerator of condition B. The default is null and the first level is treated as the numerator.}

\item{win}{Optional. Object of class \code{owin} for a custom two-dimensional window within which to estimate the surfaces. The default is NULL and calculates a convex hull around the data.}

\item{verbose}{Logical. If \code{TRUE} will print function progress during execution. If \code{FALSE} (the default), will not print.}

\item{...}{Arguments passed to \code{\link[sparr]{risk}} to select bandwidth, edge correction, and resolution.}
}
\value{
An object of class 'list' where each element is a object of class 'rrs' created by the \code{\link[sparr]{risk}} function with two additional components:

\describe{
\item{\code{rr}}{An object of class 'im' with the relative risk surface.}
\item{\code{f}}{An object of class 'im' with the spatial density of the numerator.}
\item{\code{g}}{An object of class 'im' with the spatial density of the denominator.}
\item{\code{P}}{An object of class 'im' with the asymptotic p-value surface.}
\item{\code{lrr}}{An object of class 'im' with the log relative risk surface.}
\item{\code{alpha}}{A numeric value for the alpha level used within the gate.}
}
}
\description{
Estimates a ratio of relative risk surfaces and computes the asymptotic p-value surface for a single gate with two conditions. Includes features for basic visualization. This function is used internally within the \code{\link{gating}} function to extract the points within the significant areas. This function can also be used as a standalone function.
}
\details{
This function estimates a ratio of relative risk surfaces and computes the asymptotic p-value surface for a single gate with two conditions using three successive \code{\link[sparr]{risk}} functions. A relative risk surface is estimated for Condition A at each level of Condition B and then a ratio of the two relative risk surfaces is computed.

\deqn{RR_{Condition B1} = \frac{Condition A2 of B1}{Condition A1 of B1}}
\deqn{RR_{Condition B2} = \frac{Condition A2 of B2}{Condition A1 of B2}}
\deqn{ln(rRR) = ln\left (\frac{RR_{Condition B2}}{CRR_{Condition B2}}\right )}

The p-value surface of the ratio of relative risk surfaces is estimated assuming asymptotic normality of the ratio value at each gridded knot. The bandwidth is fixed across all layers. Basic visualization is available if \code{doplot = TRUE}.

Provides functionality for a correction for multiple testing.  If \code{p_correct = "uncorrelated"}, then a conventional Bonferroni correction is calculated by dividing the \code{alpha} level by the number of gridded knots across the estimated surface. The default in the \code{\link[sparr]{risk}} function is a resolution of 128 x 128 or n = 16,384 knots and a custom resolution can be specified using the \code{resolution} argument within the \code{\link[sparr]{risk}} function. If \code{p_correct = "correlated"} (NOTE: May take a considerable amount of computation resources and time), then a Bonferroni correction that takes into account the spatial correlation of the surface is calculated within the internal \code{pval_correct} function. The \code{alpha} level is divided by the minimum number of knots that are not spatially correlated. The minimum number of knots that are not spatially correlated is computed by counting the knots that are a distance apart that exceeds the minimum distance of non-significant spatial correlation based on a correlogram using the \code{\link[pgirmess]{correlog}} function. If \code{p_correct = "none"}, then the function does not account for multiple testing and uses the uncorrected \code{alpha} level. See the internal \code{pval_correct} function documentation for more details.

The two condition variables (Condition A and Condition B) within \code{dat} must be of class 'factor' with two levels. The first level in each variable is considered the numerator (i.e., "case") value and the second level is considered the denominator (i.e., "control") value. The levels can also be specified using the \code{c1n} and \code{c2n} parameters.
}
\examples{
test_lotrrs <- lotrrs(dat = randCyto)

}
