#' R function for testing the difference in distance of two point feature datasets to
#' a target feature dataset
#'
#' The function allows to perform a permutation-based t-test to test the difference in distance
#' of two point feature datasets to a target feature dataset.
#' The latter can consist of either points, a lines, or polygons.
#'
#' Under the hood, the function relies on the perm.t.test() function out of this same package.
#' First, for each feature of both patterns, the distance to the nearest
#' target feature is calculated; for each set of features, the distances
#' are eventually averaged; the observed difference between the two averages
#' is stored. Then, the individual observed nearest distances are randomly assigned
#' to either group; the re-assignment is performed B times (999 by default) and each time
#' the difference between the two averages is calculated.
#' The distribution of these permuted average differenes represents the distribution of that
#' statistic under the Null Hypothesis of no difference in distance to the target feature.
#' One-sided and two-sided p-values are reported.\cr
#'
#' @param feat1 Point pattern to be tested (of point type; 'SpatialPointsDataFrame' class).
#' @param feat2 Second point pattern to be tested (of point type; 'SpatialPointsDataFrame' class).
#' @param to.feat Target feature (point, polyline, or polygon type; 'SpatialPointsDataFrame',
#'   'SpatialLinesDataFrame', 'SpatialPolygonsDataFrame' class).
#' @param feat1.lab Label to be used in the returned chart to indicate the 'feat1' (default: smpl 1).
#' @param feat2.lab Label to be used in the returned chart to indicate the 'feat2' (default: smpl 2).
#' @param B Desired number of permutations (set at 999 by default).
#'
#' @return The frequency histogram returned by the function displays the distribution of the
#' permuted mean difference between the two samples; a solid dot indicates the observed mean
#' difference, while an hollow dot represents the mean of the permuted differences.
#' Two dashed blue lines indicates the 0.025 and 0.975 percentile of the permuted
#' distribution. A rug plot at the bottom histgram indicates the individual permuted mean
#' differences. At the bottom of the chart, some information are displayed.
#' In particular, the observed mean difference and the permuted p-values are reported.
#' In the last row, the result of the regular (parametric) t-test (both assuming and
#' not assuming equal variances) is reported to allow users to compare the outcome of
#' these different versions of the test.
#'
#' @keywords distDiffTest
#'
#' @export
#'
#' @importFrom rgeos gDistance
#'
#' @examples
#' #test the difference in distance of two sets of points to the nearest geological fault
#' distDiffTest(feat1=springs, feat2=points, to.feat=faults, B=299)
#'
#' @seealso \code{\link{perm.t.test}}
#'
distDiffTest <- function(feat1, feat2, to.feat, feat1.lab=NULL, feat2.lab=NULL, B=999){

  #calculates all the pair-wise distances between feat1 and to-feature
  feat1.dist <- rgeos::gDistance(feat1, to.feat, byid=TRUE)

  #calculates all the pair-wise distances between feat2 and to-feature
  feat2.dist <- rgeos::gDistance(feat2, to.feat, byid=TRUE)

  #for each feat1 (i.e., column-wise), get the minimum distance to the to-feature
  obs.min.dist.feat1 <- apply(feat1.dist, 2, min)

  #for each feat2 (i.e., column-wise), get the minimum distance to the to-feature
  obs.min.dist.feat2 <- apply(feat2.dist, 2, min)

  #create a dataframe merging the two sets of minimum distances
  dtf <- data.frame(c(obs.min.dist.feat1, obs.min.dist.feat2))

  #create an empty column in the previous dataframe
  dtf$grp <- ""

  #populate the second column of the dataframe with a factor indicating feat1
  #the column with the 2 factors will be used within the perm.t.test() function
  dtf$grp[1:length(obs.min.dist.feat1)] <- "A"

  #populate the second column of the dataframe with a factor indicating feat2
  #the column with the 2 factors will be used within the perm.t.test() function
  dtf$grp[length(obs.min.dist.feat1)+1:length(obs.min.dist.feat2)] <- "B"

  #use the dataframe created above to perform the perm.t.test() function
  perm.t.test(dtf, format="long", sample1.lab=feat1.lab, sample2.lab=feat2.lab, B=B)
}
