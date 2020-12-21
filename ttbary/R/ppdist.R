#' Compute Distance Between Two Point Patterns
#'
#' Based on an arbitrary matrix of "distances" between the points of two point patterns
#' \eqn{\xi}{xi} and \eqn{\eta}{eta}, this function computes versions of the transport-transform
#' distance between \eqn{\xi}{xi} and \eqn{\eta}{eta}.
#'
#' @param dmat a matrix specifying in its \eqn{(i,j)}-th entry the distance from the
#'        i-th point of \eqn{\xi}{xi} to the j-th point of \eqn{\eta}{eta}.
#' @param penalty a positive number. The penalty for adding/deleting points. 
#' @param type either \code{"tt"}/\code{"TT"} for the transport-transform metric
#'        or \code{"rtt"}/\code{"RTT"} for the relative transport-transform metric.
#' @param ret_matching logical. Shall the optimal point matching be returned?
#' @param p a number \eqn{>0}. The matching is chosen such that the \code{p}-th
#'        order sum (\eqn{\ell_p}{l_p}-norm) is minimized.       
#' @param precision a small positive integer value. The precisions of the computations, which
#'        are currently performed in integers. After correcting for the penalty, \code{dmat^p}
#'        is divided by its largest entry, multiplied by \code{10^precision} and rounded to
#'        compute the optimal matching. The default value \code{NULL} chooses maximal
#'        integer precision possible, which is \code{precision = 9} on almost all systems.
#'     
#' @details The transport-transform (TT) distance gives the minimal total cost for \dQuote{morphing}
#'          the pattern \eqn{\xi}{xi} into the pattern \eqn{\eta}{eta} by way of shifting points (at costs
#'          specified in \code{dmat}) and adding or deleting points (each at cost \code{penalty}).
#'          The total cost is determined as 
#'          \deqn{\biggl( \sum_{j=1}^n c_j^p \biggr)^{1/p},}{(sum_{j=1}^n c_j^p)^{1/p},}
#'          where \eqn{c_j} denotes the cost for the \eqn{j}th individual operation and \eqn{n} is
#'          the cardinality of the larger point pattern.  
#'
#'          The relative transport-transform (RTT) metric is exactly the same, but the sum in the
#'          total cost is divided by the larger cardinality:
#'          \deqn{\biggl( \frac{1}{n}\sum_{j=1}^n c_j^p \biggr)^{1/p}.}{(1/n * sum_{j=1}^n c_j^p)^{1/p}.}
#' 
#'          The TT- and RTT-metrics form an umbrella concept that    
#'          includes the OSPA and Spike Time metrics frequently used in the literature.
#'          See Müller, Schuhmacher and Mateu (2019) for details.
#'          
#' @return The corresponding distance between the point patterns if \code{ret_matching}
#'         is \code{FALSE}. 
#'         
#'         Otherwise a list with components \code{dist} containing
#'         this distance and two vectors \code{target1, target2} of integers, where
#'         \code{target}\eqn{i} specifies the indices of the points in the other pattern
#'         that the points of the \eqn{i}-th pattern are matched to and \code{NA} every
#'         time a point is deleted.
#'         
#'         There may be a minus in front of an index, where
#'         \code{-j} indicates that the corresponding pairing with point \code{j}
#'         would be over a distance of more than \eqn{2^{1/p} \cdot 
#'         \code{penalty}}{2^{1/p} * penalty}. This is 
#'         equivalent to saying that the corresponding point of the first pattern
#'         is deleted and the \eqn{j}-th point of the second pattern is added. 
#'         
#'         Note that having more than one minus implies that the matching is non-unique.
#'    
#' @references Raoul Müller, Dominic Schuhmacher and Jorge Mateu (2019).\cr
#'             Metrics and Barycenters for Point Pattern Data.\cr
#'             Preprint \href{https://arxiv.org/abs/1909.07266}{arXiv:1909.07266}
#'                   
#' @author Dominic Schuhmacher \email{schuhmacher@math.uni-goettingen.de} 
#'
#' @examples
#'   # small example
#'   # -------------
#'   set.seed(181230)
#'   xi <- spatstat::rpoispp(20)
#'   eta <- spatstat::rpoispp(20)
#'   dmat <- spatstat::crossdist(xi,eta)
#'   res <- ppdist(dmat, penalty=1,  type="rtt", ret_matching=TRUE, p=1)
#'   plotmatch(xi, eta, dmat, res, penalty=1, p=1)
#'   res$dist
#' 
#'   # for comparison: ospa-distance computation from spatstat:
#'   res_ospa <- spatstat::pppdist(xi,eta,"spa")
#'   res_ospa$distance  # exactly the same as above because nothing gets cut off 
#' 
#' 
#'   # same example, but with a much smaller penalty for adding/deleting points
#'   # --------------------------------------------------------------- 
#'   res <- ppdist(dmat, penalty=0.1,  type="rtt", ret_matching=TRUE, p=1)
#'   plotmatch(xi, eta, dmat, res, penalty=0.1, p=1)
#'     # dashed lines indicate points that are deleted and re-added at new position
#'     # grey segments on dashed lines indicate the cost of deletion plus re-addition
#'   res$dist
#'   
#'   # for comparison: ospa-distance computation from spatstat
#'   # (if things do get cut off, we have to ensure that the cutoff distances
#'   # are the same, thus cutoff = 2^(1/p) * penalty):
#'   res_ospa <- spatstat::pppdist(xi,eta,"spa",cutoff=0.2)
#'   res_ospa$distance  # NOT the same as above
#'   res_ospa$distance - abs(xi$n-eta$n) * 0.1 / max(xi$n,eta$n)  # the same as above
#'   
#'   
#'   # a larger example
#'   # --------------------------------------------------------------- 
#'   set.seed(190203)
#'   xi <- spatstat::rpoispp(2000)
#'   eta <- spatstat::rpoispp(2000)
#'   dmat <- spatstat::crossdist(xi,eta)
#'   res <- ppdist(dmat, penalty = 0.1,  type = "rtt", ret_matching = TRUE, p = 1)
#'   res$dist
#'   # takes about 2-3 seconds
#'
#' @aliases ttdist
#'
#' @export

ppdist <- function(dmat, penalty = 1, type = c("tt","rtt","TT","RTT"), 
                   ret_matching = FALSE, p = 1, precision = NULL) {
  type <- match.arg(type)
  relative <- (type == "rtt" | type == "RTT")   

  return(ttdist(dmat, penalty = penalty, relative = relative, 
                ret_matching = ret_matching, p = p, precision = precision))
}


ttdist <- function(dmat, penalty = 1, relative = FALSE, 
                   ret_matching = FALSE, p = 1, precision = NULL) {
  stopifnot(is.matrix(dmat))
  
  if (is.null(precision)) {
    precision <- trunc(log10(.Machine$integer.max))
  }
  
  n1 <- dim(dmat)[1]
  n2 <- dim(dmat)[2]
  n <- max(n1,n2)
    
  if (n1 == 0 || n2 == 0) {
    return(zeropoints(n1,n2,(n^(1/p))*penalty,ret_matching))
  }
    
  dfix <- matrix(penalty,n,n)
  dfix[1:n1,1:n2] <- dmat
  # dfix contains d-distances and C for distance to aleph points
  d <- dfix2 <- pmin(dfix,2^(1/p)*penalty)
  # dfix2 contains d'-distances
  # d is a working copy of dfix2 that will be rescaled.
  d <- d/max(d)
  d <- round((d^p) * (10^precision))
    
  epsvec <- prepare_epsvec(d,n)
  neps <- length(epsvec)
  stopifnot(neps >= 1)
  d <- max(d)-d
  # auctionbf uses a "desire matrix"
    
  res <- auctionbf2cpp(d,n,pers_to_obj = rep(-1,n),
                       obj_to_pers = rep(-1,n),
                       price       = rep(0,n),
                       profit      = rep(0,n),
                       neps,
                       epsvec)
  
  target1pre <- res$pers_to_obj+1
  
  distvec <- dfix2[cbind(1:n,target1pre)]
  dist <- Lpsum(distvec, p)
  if (relative) {
    dist <- dist/n^{1/p}
  }
    
  if (!ret_matching) {
    return(dist)  # return tt/rtt distance
  } else {
    target1 <- target1pre[1:n1]
    rawdistvec <- dfix[cbind(1:n1,target1[1:n1])]
    ind1 <- rawdistvec > 2^(1/p)*penalty
    target1[ind1] <- -target1[ind1]
    target2 <- res$obj_to_pers[1:n2]+1
    rawdistvec <- dfix[cbind(target2,1:n2)]
    ind2 <- rawdistvec > 2^(1/p)*penalty
    target2[ind2] <- -target2[ind2]
    if (n1 < n2) {
      target2[target2 > n1] <- NA
    } else if (n1 > n2) {
      target1[target1 > n2] <- NA 
    }
      
    return(list(dist=dist, target1=target1, target2=target2))
  }  
}  

#' Compute Sum of q-th Powers of Distances Between a Point Pattern and a List of Point Patterns
#'
#' Determine the Euclidean distance based TT-p-distances (or RTT-p-distances) between
#' a single point pattern \code{zeta} and each point pattern in a list \code{pplist}. Then
#' compute the sum of \eqn{q}-th powers of these distances.
#'
#' @param zeta an object of class \code{\link[=ppp]{ppp}}.
#' @param pplist an object of class \code{\link[=solist]{ppplist}} or an 
#'        object that can be coerced to this class, such as a list of \code{\link[=ppp]{ppp}}
#'        objects.  
#' @param penalty a positive number. The penalty for adding/deleting points. 
#' @param type either \code{"tt"}/\code{"TT"} for the transport-transform metric
#'        or \code{"rtt"}/\code{"RTT"} for the relative transport-transform metric.
#' @param p a number \eqn{>0}. Matchings between \code{zeta} and the patterns in  
#'        \code{pplist} are chosen such that the \code{p}-th order sums (\eqn{\ell_p}{l_p}-norms)
#'        of the Euclidean distances are minimized.
#' @param q a number \eqn{>0}. 
#'        
#' @return A nonnegative number, the \code{q}-th order sum of the TT-p- or RTT-p-distances 
#'         between \code{zeta} and each pattern in \code{pplist}. This number has an attribute
#'         \code{distances} that contains the individual distances.
#'         
#' @details The main purpose of this function is to evaluate the relative performance
#'         of approximate \eqn{q}-th order barycenters of point patterns. A true 
#'         \eqn{q}-th order barycenter of the point patterns \eqn{\xi_1,\ldots,\xi_k}{xi_1, ..., xi_k}
#'         with respect to the TT-p metric \eqn{\tau_p}{tau_p} minimizes 
#'         \deqn{\sum_{j=1}^k \tau_p(\xi_j, \zeta)^q}{sum_{j=1}^k tau_p(xi_j, zeta)^q}
#'         in \eqn{\zeta}{zeta}.      
#'         
#'         The most common choices are \code{p = q = 1} and \code{p = q = 2}. Other
#'         choices have not been tested.
#'                
#' @author Dominic Schuhmacher \email{schuhmacher@math.uni-goettingen.de} 
#'
#' @examples
#'   # See the examples for kmeansbary
#'   
#' @seealso \code{\link{ppdist}} for computation of TT-p- and RTT-p-metrics,\cr
#'          \code{\link{kmeansbary}} for finding a local minimum of the above sum for \code{p = q = 2} 
#'   
#' @export
sumppdist <- function(zeta, pplist, penalty = 1, type = c("tt","rtt","TT","RTT"),
                       p = 1, q = 1) {
  if (!("ppplist" %in% class(pplist))) {
    pplist <- as.ppplist(pplist)  # allows e.g. also to pass a single point pattern as pplist
  }
  
  if (p != q) {
    warning("Setting p=q might give more desirable results.")
  }
  
  n <- length(pplist)
  disttozeta <- rep(0,n)
  for (i in 1:n) {
    dmat <- crossdist(zeta, pplist[[i]])
    disttozeta[i] <- ppdist(dmat, penalty = penalty, type = type, ret_matching = FALSE, p = p)
  }
  
  res <- Lpsumraw(disttozeta, q)
  attr(res, "distances") <- disttozeta
  # print(disttozeta^p)  # not sure why we printed this with ^p
  return(res)
}

#' Compute Distance Between Two Point Patterns on a Network
#'
#' Based on an arbitrary matrix of "distances" on a network, this function computes versions
#' of the transport-transform distance between two point patterns \eqn{\xi}{xi} and \eqn{\eta}{eta}
#' on this network.
#'
#' @param dmat a matrix specifying in its \eqn{(i,j)}-th entry the shortest-path distance from the
#'        i-th point of \eqn{\xi}{xi} to the j-th point of  \eqn{\eta}{eta}
#'        OR the distance matrix of a whole network. In the latter case arguments \eqn{\xi}{xi} and
#'        \eqn{\eta}{eta} have to be specified.
#' @param xi a vector specifying the vertex-indices of \eqn{\xi}{xi}, only needed if \code{dmat} is the 
#'           distance matrix of a whole network.
#' @param eta a vector specifying the vertex-indices of \eqn{\eta}{eta}, only needed if \code{dmat} is the 
#'           distance matrix of a whole network.
#' @param penalty a positive number. The penalty for adding/deleting points. 
#' @param type either \code{"tt"}/\code{"TT"} for the transport-transform metric
#'        or \code{"rtt"}/\code{"RTT"} for the relative transport-transform metric.
#' @param ret_matching Logical. Shall the optimal point matching be returned?
#' @param p a number \eqn{>0}. The matching is chosen such that the \code{p}-th
#'        order sum (\eqn{\ell_p}{l_p}-norm) is minimized.       
#' @param precision a small positive integer value. The precision of the computations, which
#'        are currently performed in integers. After correcting for the penalty, \code{dmat^p}
#'        is divided by its largest entry, multiplied by \code{10^precision} and rounded to
#'        compute the optimal matching. The default value \code{NULL} chooses maximal
#'        integer precision possible, which is \code{precision = 9} on almost all systems.
#'     
#' @details This function provides a more convenient way for computing (relative)
#'          transport-transform distances on networks if the points of the patterns are given in
#'          terms of indices of network vertices. If \code{dmat} contains only the distances
#'          between the points of \eqn{\xi}{xi} and \eqn{\eta}{eta}, this function
#'          does the same as \code{\link{ppdist}}.
#'          
#' @return The corresponding distance between the point patterns if \code{ret_matching}
#'         is \code{FALSE}. 
#'         
#'         Otherwise a list with components \code{dist} containing
#'         this distance and two vectors \code{target1, target2} of integers, where
#'         \code{target}\eqn{i} specifies the indices of the points in the other pattern
#'         that the points of the \eqn{i}-th pattern are matched to and \code{NA} every
#'         time a point is deleted.
#'
#'         There may be a minus in front of an index, where
#'         \code{-j} indicates that the corresponding pairing with point \code{j}
#'         would be over a distance of more than \eqn{2^{1/p} \cdot 
#'         \code{penalty}}{2^{1/p} * penalty}. This is 
#'         equivalent to saying that the corresponding point of the first pattern
#'         is deleted and the \eqn{j}-th point of the second pattern is added. 
#'         
#'         Note that having more than one minus implies that the matching is non-unique.
#'         
#' @author Raoul Müller  \email{raoul.mueller@uni-goettingen.de}\cr
#'         Dominic Schuhmacher \email{schuhmacher@math.uni-goettingen.de}
#'
#' @examples 
#'   set.seed(123456)
#'   nvert <- 100 #number of vertices in the network
#'   lambda <- 0.5 #expected number of points per unit length
#'   ln <- delaunayNetwork(runifpoint(nvert)) #create an artificial network
#'   ppnetwork <- rpoislpp(lambda, ln, nsim = 2)
#'     #simulate two point patterns on the network
#' 
#'   plot(ppnetwork[[1]]$domain, cex=0.5, main="")
#'   plot(as.ppp(ppnetwork[[1]]),vpch=1,col=2,add=TRUE)
#'   plot(as.ppp(ppnetwork[[2]]),vpch=1,col=4,add=TRUE)
#'
#'   res <- netsplit(ln, ppnetwork)
#'     #incorporate data point patterns into the network
#'     #calculating all pairwise distances between vertices
#'     #and creating matrix of vertex-indices of data point patterns
#'   
#'   xi <- res$ppmatrix[1:npoints(ppnetwork[[1]]), 1]
#'   eta <- res$ppmatrix[1:npoints(ppnetwork[[2]]), 2]
#'   res2 <- ppdistnet(res$network$dpath, xi = xi, eta = eta,
#'                     penalty = 1, type = "tt", ret_matching = TRUE, p = 1)
#'   res2
#'
#' @aliases ttdistnet 
#' 
#' @seealso \code{\link{ppdist}}
#'
#' @export
#' 
ppdistnet <- function(dmat, xi=NULL, eta=NULL, penalty = 1, type = c("tt","rtt","TT","RTT"),
                      ret_matching = FALSE, p = 1, precision = NULL) {
  type <- match.arg(type)
  relative <- (type == "rtt" | type == "RTT")
  
  if (is.null(xi) || is.null(eta)) {
    #just work with distance matrix dmat
    return(ttdist(dmat, penalty = penalty, relative = relative, 
                  ret_matching = ret_matching, p = p, precision = precision))
  }else{
    if (max(eta) > length(dmat[1,]) || max(xi) > length(dmat[,1])) {
      n1 <- length(dmat[,1]) #number of rows
      n2 <- length(dmat[1,]) #number of columns
      n3 <- max(c(xi,eta))
      
      dmat <- rbind(dmat, matrix(penalty^(1/p),nrow = n3-n1,ncol = n2))
      dmat <- cbind(dmat, rbind(matrix(penalty^(1/p),nrow = n1,ncol = n3-n2),matrix(0,nrow = n3-n1,ncol = n3-n2)))
    }
    dmat <- dmat[xi,eta]
    res <- ttdist(dmat, penalty = penalty, relative = relative, ret_matching = ret_matching, p = p, precision = precision)
    dist <- res$dist
    target1 <- eta[res$target1]
    target2 <- xi[res$target2]
    return(list(dist = dist, target1 = target1, target2 = target2))
  }
}

#' Compute Sum of q-th Powers of Distances Between a Point Pattern and a Collection of Point Patterns on a Network
#'
#' Based on the shortest-path metric in a network, determine the TT-p-distances (or RTT-p-distances) 
#' between a single point pattern \code{zeta} and a collection of point patterns. Then
#' compute the sum of \eqn{q}-th powers of these distances. The point patterns are
#' specified by vectors of indices referring to the vertices in the network. 
#'
#' @param dmat the distance matrix of a network containing all shortest-path distances between its vertices.
#' @param zeta a vector specifying the vertex-indices of zeta.
#' @param ppmatrix a matrix specifying in its columns the vertex-indices of the point patterns in the collection.
#'        A virtual index that is one greater than the maximum vertex-index in the network
#'        can be used to fill up columns so that they all have the same length.
#' @param penalty a positive number. The penalty for adding/deleting points. 
#' @param type either \code{"tt"}/\code{"TT"} for the transport-transform metric
#'        or \code{"rtt"}/\code{"RTT"} for the relative transport-transform metric.
#' @param p a number \eqn{>0}. Matchings between \code{zeta} and the patterns in   
#'        \code{ppmatrix} are chosen such that the \code{p}-th order sums (\eqn{\ell_p}{l_p}-norms)
#'        of the shortest-path distances are minimized.
#' @param q a number \eqn{>0}. 
#'         
#' @return A nonnegative number, the \code{q}-th order sum of the TT-p- or RTT-p-distances 
#'         between the patterns represented by \code{zeta} and \code{ppmatrix}. This number has an attribute
#'         \code{distances} that contains the individual distances.
#'
#' @details The main purpose of this function is to evaluate the relative performance
#'          of approximate \eqn{q}-th order barycenters of point patterns. A true 
#'          \eqn{q}-th order barycenter of the point patterns \eqn{\xi_1,\ldots,\xi_k}{xi_1, ..., xi_k}
#'          with respect to the TT-p metric \eqn{\tau_p}{tau_p} minimizes 
#'          \deqn{\sum_{j=1}^k \tau_p(\xi_j, \zeta)^q}{sum_{j=1}^k tau_p(xi_j, zeta)^q}
#'          in \eqn{\zeta}{zeta}.      
#'         
#'         The most common choices are \code{p = q = 1} and \code{p = q = 2}. Other
#'         choices have not been tested.
#' 
#' @author Raoul Müller  \email{raoul.mueller@uni-goettingen.de}\cr
#'         Dominic Schuhmacher \email{schuhmacher@math.uni-goettingen.de}
#'
#' @examples
#'   # See examples for kmeansbarynet
#'   
#' @seealso \code{\link{kmeansbarynet}}, \code{\link{sumppdist}}
#'   
#' @export
sumppdistnet <- function(dmat, zeta, ppmatrix, penalty = 1, type = c("tt","rtt","TT","RTT"),
                      p = 1, q = 1) {

  if (p != q) {
    warning("Setting p=q might give more desirable results.")
  }
  if (max(ppmatrix) > length(dmat[1,])) {
    n1 <- length(dmat[,1]) #number of rows
    n2 <- length(dmat[1,]) #number of columns
    n3 <- max(ppmatrix)
    
    dmat <- rbind(dmat, matrix(penalty^(1/p),nrow = n3-n1,ncol = n2))
    dmat <- cbind(dmat, rbind(matrix(penalty^(1/p),nrow = n1,ncol = n3-n2),matrix(0,nrow = n3-n1,ncol = n3-n2)))
  }
  
  n <- length(ppmatrix[1,])
  disttozeta <- rep(0,n)
  for (i in 1:n) {
    dmattemp <- dmat[zeta,ppmatrix[,i]]
    disttozeta[i] <- ppdistnet(dmattemp, penalty = penalty, type = type, ret_matching = FALSE, p = p)
  }
  
  res <- Lpsumraw(disttozeta, q)
  attr(res, "distances") <- disttozeta
  # print(disttozeta^p)  # not sure why we printed this with ^p
  return(res)
}

# -----------------------------------------------

#' Plot Optimal Matching between Two Point Patterns
#'
#' After calling \code{\link{ppdist}} with argument \code{ret_matching = TRUE}
#' in a situation where it makes sense to assign to the points of the patterns \eqn{\xi}{xi}
#' and \eqn{\eta}{eta} coordinates in \eqn{R^2}{R^2}, this function may be used to display
#' the result graphically.
#' 
#' @param xi,eta objects of class \code{\link[=ppp]{ppp}}.
#' @param dmat a matrix specifying in its \eqn{(i,j)}-th entry the distance from the
#'        i-th point of \eqn{\xi}{xi} to the j-th point of \eqn{\eta}{eta}.
#' @param res the object returned by the call to \code{\link{ppdist}} with \code{ret_matching = TRUE}.
#' @param penalty a positive number. The penalty for adding/deleting points. 
#' @param p a number \eqn{>0}. The order of the TT- or RTT-distance computed.
#' @param cols,pchs,cexs vectors of length 2 specifying the corresponding graphic 
#'        parameters col, pch and cex for plotting the two point patterns.     
#' @param ... further graphic parameters passed to the code that draws the line segments
#'        between the points.
#'     
#' @details The default use-case is to plot a matching obtained with \code{\link{ppdist}}.
#'        In that case \code{dmat}, \code{penalty} and \code{p} should be the same
#'        as in the call to \code{ppdist}. These objects are used to display additional
#'        information about the matching.
#'          
#' @return Used for the side effect of plotting.
#'                   
#' @author Dominic Schuhmacher \email{schuhmacher@math.uni-goettingen.de} 
#'
#' @examples
#'   # See examples for ppdist
#'
#' @seealso \code{\link{ppdist}}
#'         
#' @export
#' 
plotmatch <- function(xi, eta, dmat, res, penalty, p=1, cols=c(2,4), pchs=c(1,1), cexs=c(1,1), ...) {
  stopifnot(is.matrix(dmat))
  n1 <- dim(dmat)[1]
  n2 <- dim(dmat)[2]
  
  # if only col, pch, cex is specified we get strange 
  # behaviour due to partial matching, therefore:
  if (length(cols) == 1) {cols <- c(cols,cols)}
  if (length(pchs) == 1) {pchs <- c(pchs,pchs)}
  if (length(cexs) == 1) {cexs <- c(cexs,cexs)}
  
  #oldpar <- par(mai=c(0,0,0.2,0))
  plot(xi, main="", cols=cols[1], pch=pchs[1], cex=cexs[1])
  plot(eta, add=TRUE, cols=cols[2], pch=pchs[2], cex=cexs[2])
  
  if (n1 <= n2) {
    tt1 <- abs(res$target1)
    segments(xi$x,xi$y,eta$x[tt1],eta$y[tt1], lty=2-sign(res$target1), ...)
    ## full line for + targets, dashed line for - targets
    ind <- which(res$target1 < 0)
    jnd <- tt1[ind]
  } else {
    tt2 <- abs(res$target2)
    segments(xi$x[tt2],xi$y[tt2],eta$x,eta$y, lty=2-sign(res$target2), ...)
    ## full line for + targets, dotted line for - targets
    jnd <- which(res$target2 < 0)
    ind <- tt2[jnd]
  }
  dis <- dmat[cbind(ind,jnd)]
  vecs <- (cbind(eta$x[jnd],eta$y[jnd]) - cbind(xi$x[ind],xi$y[ind]))/dis
  centr <- (cbind(eta$x[jnd],eta$y[jnd]) + cbind(xi$x[ind],xi$y[ind]))/2
  pen <- 2^(1/p) * penalty / 2
  a <- centr - pen*vecs
  b <- centr + pen*vecs
  arrows(a[,1],a[,2],b[,1],b[,2],length=0.05,angle=90,code=3,lwd=1.5,col=grey(0.7))
  #par(oldpar)
  return(invisible())
}




# for treating the special cases
# n1: no. of points in first pattern
# n2: no. of points in second pattern
# rvalue: which value to return if only one n1,n2 is zero 
# ret_matching: shall matching be returned
zeropoints <- function(n1,n2,rvalue,ret_matching) {
  if (!ret_matching) {
    if (n1 == 0 && n2 == 0) {
      return(0)
    } else {
      return(rvalue)
    }
  } else {
    if (n1 == 0 && n2 == 0) {
      return(list(dist=0, target1=numeric(0), target2=numeric(0)))
    } else if (n1 == 0) {
      return(list(dist=rvalue, target1=numeric(0), target2=rep(NA,n2)))
    } else {
      return(list(dist=rvalue, target1=rep(NA,n1), target2=numeric(0)))
    }
  }
}


# prepares epsilon vector for epsilon-scaling
# this should definitely go into the C or C++ code
prepare_epsvec <- function(d,n) {
  dupper <- max(d)/10	
  lasteps <- 1/(n+1)
  epsfac <- 10
  epsvec <- lasteps
  ## Bertsekas: from dupper/2 to 1/(n+1) divide repeatedly by a constant
  while (lasteps < dupper) {
    lasteps <- lasteps*epsfac
    epsvec <- c(epsvec,lasteps)
  }
  epsvec <- rev(epsvec)[-1]
  return(epsvec)
}


Lpmean <- function(x, p) {
  if (p == 1) {
    return(mean(x))
  } else {
    f <- max(x)
    return(ifelse(f==0, 0, f * mean((x/f)^p)^(1/p)))
  }
}


Lpsum <- function(x, p) {
  if (p == 1) {
    return(sum(x))
  } else {
    f <- max(x)
    return(ifelse(f==0, 0, f * sum((x/f)^p)^(1/p)))
  }
}


Lpsumraw <- function(x, p) {
  if (p == 1) {
    return(sum(x))
  } else {
    f <- max(x)
    return(ifelse(f==0, 0, f^p * sum((x/f)^p)))
  }
}
