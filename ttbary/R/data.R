#' Simulated Point Pattern Lists
#'
#' Lists of simulated point patterns for illustrating the computation of barycenters.
#'
#' @name pplist-data
#' @rdname pplist-data
#' @aliases pplist_samecard
#' @aliases pplist_diffcard
#' 
#' @format Objects of class \code{pplist}, which are essentially lists of \code{ppp}-objects.
#' 
#' @details
#' \code{pplist_samecard} contains 80 point patterns of 100 points each. The patterns
#' were independently generated from a distribution that creates quite distinctive clusters.
#' 
#' \code{pplist_diffcard} contains 50 point patterns with cardinalities ranging from 17 to 42. 
#' The patterns were independently generated from a distribution that creates overlapping clusters.
#' 
#' @author Raoul Müller  \email{raoul.mueller@uni-goettingen.de}\cr
#'         Dominic Schuhmacher \email{schuhmacher@math.uni-goettingen.de}
#'         
#' @examples
#' # plot the first eight patterns of each data set
#' plot(superimpose(pplist_samecard[1:8]), legend=FALSE, cex=0.4, cols=1:8)
#' plot(superimpose(pplist_diffcard[1:8]), legend=FALSE, cex=0.4, cols=1:8)
#'               
#' @keywords datasets
NULL

#' @rdname pplist-data
#' @docType data
#' @keywords datasets
#' @export
"pplist_samecard"

#' @rdname pplist-data
#' @docType data
#' @keywords datasets
#' @export
"pplist_diffcard"

