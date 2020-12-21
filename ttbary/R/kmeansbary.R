#' Compute Pseudo-Barycenter of a List of Point Patterns
#'
#' Starting from an initial candidate point pattern \code{zeta}, use a k-means-like
#' algorithm to compute a local minimum in the barycenter problem based on the TT-2 metric
#' for a list \code{pplist} of planar point patterns. 
#'
#' @param zeta a point pattern. Object of class \code{ppp} or a list with components \code{x} and \code{y}.
#' @param pplist a list of point patterns. Object of class \code{ppplist} or any list where each elements
#'           has components  \code{x} and \code{y}.
#' @param penalty the penalty for adding/deleting points when computing TT-2 distances.
#' @param add_del for how many iterations shall the algorithm add points to / delete points from zeta
#'           if this is favorable? Defaults to Inf.
#' @param surplus by how many points is the barycenter point pattern allowed to be larger than
#'           the largest input point pattern (among pplist and zeta) if add_del > 0.
#'           A larger number increases the computational load.
#' @param N the maximum number of iterations.
#' @param eps the algorithm stops if the relative improvement of the objective function between two
#'          iterations is less than eps.
#' @param verbose the verbosity level. One of 0, 1, 2, 3, where 0 means silent and 3 means full details.
#'
#'
#' @details  Given \eqn{k} planar point patterns \eqn{\xi_1, \ldots, \xi_k}{xi_1, ..., xi_k} (stored in
#'           \code{pplist}), this function finds a local minimizer \eqn{\zeta^*}{zeta*} of
#'           \deqn{\sum_{j=1}^k \tau_2(\xi_j, \zeta)^2,}{sum_{j=1}^k tau_2(xi_j, zeta)^2,}
#'           where \eqn{\tau_2}{tau_2} denotes the TT-2 metric based on the Euclidean metric between points.
#'           
#'           Starting from an initial candidate point pattern \code{zeta}, the algorithm alternates
#'           between assigning a point from each pattern \eqn{\xi_j}{xi_j} to each point of the candidate
#'           and computing new candidate patterns by shifting, adding and deleting points.
#'           A detailed description of the algorithm is given in Müller, Schuhmacher and Mateu (2019).
#'           
#'           For first-time users it is recommended to keep the default values and set \code{penalty}
#'           to a noticeable fraction of the diameter of the observation window, e.g. between
#'           0.1 and 0.25 times this diameter.
#'
#' @return A list with components:
#'         \item{cost}{the sum of squared TT-2 distances between the computed pseudo-barycenter and the point patterns.}
#'         \item{barycenter}{the pseudo-barycenter as a \code{ppp}-object.}
#'         \item{iterations}{the number of iterations required until convergence.}
#'
#' @references Raoul Müller, Dominic Schuhmacher and Jorge Mateu (2019).\cr
#'             Metrics and Barycenters for Point Pattern Data.\cr
#'             Preprint \href{https://arxiv.org/abs/1909.07266}{arXiv:1909.07266}
#'
#' @author Raoul Müller  \email{raoul.mueller@uni-goettingen.de}\cr
#'         Dominic Schuhmacher \email{schuhmacher@math.uni-goettingen.de}
#'
#' @examples
#' data(pplist_samecard)
#' plot(superimpose(pplist_samecard), cex=0.7, legend=FALSE,
#'      xlim=c(-0.2,1.2), ylim=c(-0.1,1.1), main="", use.marks=FALSE) #plotting the data
#'
#' set.seed(12345)
#' zeta <- ppp(runif(100), runif(100))
#' plot(zeta, add=TRUE, col="beige", lwd=2, pch=16) #plotting the start-zeta over the data
#'
#' res <- kmeansbary(zeta, pplist_samecard, penalty=0.1, add_del=Inf)
#' plot(res$barycenter, add=TRUE, col="blue", pch=16) #adding the computed barycenter in blue
#'
#' res$cost
#' #[1] 30.30387
#' sumppdist(res$barycenter, pplist_samecard, penalty=0.1, type="tt", p=2, q=2)
#' #[1] 30.30387
#' #attr(,"distances")
#' #[1] 0.5991515 0.6133397 0.6040680 0.6020058 0.5648000 0.6415018 0.6385394 0.5784291 0.5985299
#' #[10] 0.6313200 0.7186154 ...
#' 
#' @seealso \code{\link{kmeansbaryeps}} for a variant with epsilon relaxation that is typically faster
#'
#' @export
#'
kmeansbary <- function(zeta, pplist, penalty, add_del = Inf,
                       surplus = 0, N = 200L, eps = 0.005, verbose = 0) {
  stopifnot(class(zeta) == "ppp")
  stopifnot(any(class(pplist) == "list"))
  stopifnot(all(sapply(pplist, function(x) {class(x) == "ppp"})))
  npplist <-sapply(pplist, npoints)
  nzeta <- npoints(zeta)
  n <- max(c(nzeta, npplist))
  k <- length(pplist)
  stopifnot(k > 1)

  zetax <- c(zeta$x, rep(NA, n-nzeta+surplus))
  zetay <- c(zeta$y, rep(NA, n-nzeta+surplus))

  ppmatx <- c(pplist[[1]]$x, rep(NA, n-npplist[1]+surplus))
  ppmaty <- c(pplist[[1]]$y, rep(NA, n-npplist[1]+surplus))
  for (i in 2:k) {
    ppmatx <- cbind(ppmatx, c(pplist[[i]]$x, rep(NA, n-npplist[i]+surplus)))
    ppmaty <- cbind(ppmaty, c(pplist[[i]]$y, rep(NA, n-npplist[i]+surplus)))
  }

  if (add_del > .Machine$integer.max) add_del <- .Machine$integer.max
    # this is in fact only to catch the Inf value (since it is an int in C++ not clear if
    # it would work with R_PosInf)
  res <- .Call(`_ttbary_kMeansBary`, zetax, zetay, ppmatx, ppmaty, penalty, add_del, N, eps, verbose)
  win <- bounding.box.xy(c(superimpose(pplist)$x,res$baryx),c(superimpose(pplist)$y,res$baryy))
  #
  #xwin <- c(min(win$xrange[1],0),max(win$xrange[2],1))
  #ywin <- c(min(win$yrange[1],0),max(win$yrange[2],1))
  #win <- owin(xwin,ywin)
  #
  baryx <- na.omit(res$barycenterx)
  baryy <- na.omit(res$barycentery)
  if (!isTRUE(all.equal(attr(baryx, "na.action"), attr(baryy, "na.action")))) {
    warning("NA coordinates of barycenter did not agree")
    message("NA in x-coords: ", paste(which(is.na(baryx)), collapse=" "))
    message("NA in y-coords: ", paste(which(is.na(baryy)), collapse=" "))
  }
  return(list(cost=res$cost, barycenter=ppp(baryx, baryy, window=win), iterations=res$iterations))
}

#' Compute Pseudo-Barycenter of a List of Point Patterns (with epsilon-relaxation)
#'
#' Starting from an initial candidate point pattern \code{zeta}, use a k-means-like
#' algorithm to compute a local minimum in the barycenter problem based on the TT-2 metric
#' for a list \code{pplist} of planar point patterns.
#'
#' @param epsvec a vector containing the values for the relaxed assignment. Last entry should be < 1/n, where n is the largest cardinality among the point 
#'              patterns. Otherwise the algorithm has no guarantee of terminating in a local minimum!
#'              If epsvec[1] is too small, the computational load may be large.
#'              If in doubt, choose c(10^8,10^7,10^6,...,10/(n+1),1/(n+1)).
#' @param zeta a point pattern. Object of class \code{ppp} or a list with components \code{x} and \code{y}.
#' @param pplist a list of point patterns. Object of class \code{ppplist} or any list where each elements
#'           has components  \code{x} and \code{y}.
#' @param penalty the penalty for adding/deleting points when computing TT-2 distances.
#' @param add_del for how many iterations shall the algorithm add points to / delete points from zeta
#'           if this is favorable? Defaults to Inf.
#' @param surplus By how many points is the barycenter point pattern allowed to be larger than
#'           the largest input point pattern (among pplist and zeta) if add_del > 0.
#'           A larger number increases the computational load.
#' @param relaxVec a vector of four integers controlling the epsilon-relaxation of the assignments.
#'           See details below.
#' @param N the maximum number of iterations.
#' @param eps the algorithm stops if the relative improvement of the objective function between two iterations is less
#'          than eps.
#' @param verbose the verbosity level. One of 0, 1, 2, 3, where 0 means silent and 3 means full details.
#'
#'
#' @details  Given \eqn{k} planar point patterns \eqn{\xi_1, \ldots, \xi_k}{xi_1, ..., xi_k} (stored in
#'           \code{pplist}), this function finds a local minimizer \eqn{\zeta^*}{zeta*} of
#'           \deqn{\sum_{j=1}^k \tau_2(\xi_j, \zeta)^2,}{sum_{j=1}^k tau_2(xi_j, zeta)^2,}
#'           where \eqn{\tau_2}{tau_2} denotes the TT-2 metric based on the Euclidean metric between points.
#'           
#'           Starting from an initial candidate point pattern \code{zeta}, the algorithm alternates
#'           between assigning a point from each pattern \eqn{\xi_j}{xi_j} to each point of the candidate
#'           and computing new candidate patterns by shifting, adding and deleting points.
#'           A detailed description of the algorithm is given in Müller, Schuhmacher and Mateu (2019).
#'           
#'           For first-time users it is recommended to keep the default values and set \code{penalty}
#'           to a noticeable fraction of the diameter of the observation window, e.g. between
#'           0.1 and 0.25 times this diameter.
#'           
#'           The argument \code{relaxVec} must be a vector of four integers c(a,b,c,d) > c(0,0,0,0).
#'           For the first \code{a} iterations step by step one entry of \code{epsvec} is additionally considered in the assignment, starting with
#'           only the first entry in the first iteration. In this \code{a} iterations the algorithm can stop if it has improved by less than \code{eps} between iterations. 
#'           After \code{a} iterations all entries of \code{epsvec} before \code{epsvec[b]} are ignored and everytime
#'           the algorithm does not improve, the next \code{d} entries of epsvec are additionally considered in the following iterations. When the last 
#'           entry of \code{epsvec} is considered in the assignments, the entries of epsvec before \code{epsvec[c]} are ignored.
#'           \code{relaxVec} defaults to c(20,1,1,1) meaning that in every one of the first 20 iterations one additional entry of epsvec
#'           is considered until the algorithm converges. This allows the algorithm to converge before the full epsvec was considered! For further
#'           details see example.
#'           
#'           \strong{Warning:} The argument \code{relaxVec} offers many different options for controlling the epsilon-relaxation of the assignments
#'           in order to save computation time. But choosing bad parameters may heavily increase the computational load!
#'           If in doubt, go with c(length(epsvec),1,1,1) (see examples).
#'           
#' @return A list with components:
#'         \item{cost}{the sum of squared TT-2 distances between the computed pseudo-barycenter and the point patterns.}
#'         \item{barycenter}{the pseudo-barycenter as a \code{ppp}-object.}
#'         \item{iterations}{the number of iterations required until convergence.}
#'
#' @references Raoul Müller, Dominic Schuhmacher and Jorge Mateu (2019).\cr
#'             Metrics and Barycenters for Point Pattern Data.\cr
#'             Preprint \href{https://arxiv.org/abs/1909.07266}{arXiv:1909.07266}
#'
#' @author Raoul Müller  \email{raoul.mueller@uni-goettingen.de}\cr
#'         Dominic Schuhmacher \email{schuhmacher@math.uni-goettingen.de}
#'
#' @examples
#' data(pplist_samecard)
#' plot(superimpose(pplist_samecard), cex=0.7, legend=FALSE,
#'      xlim=c(-0.2,1.2), ylim=c(-0.1,1.1), main="", use.marks=FALSE) #plotting the data
#'
#' set.seed(12345)
#' zeta <- ppp(runif(100),runif(100))
#' plot(zeta, add=TRUE, col="beige", lwd=2, pch=16) #plotting the start-zeta over the data
#' 
#' epsvec <- c(1e8,1e7,1e6,1e5,1e4,1e3,1e2,10,1,10/101,1/101)
#' 
#' relaxVec1 <- c(length(epsvec),1,1,1) 
#' #One additional entry of epsvec is considered in each iteration;
#' #algorithm can stop before full epsvec was used.
#' #Runs fast with little to no drawback in the quality of the computed solution.
#' #Time advantage more visible for large datasets.
#' 
#' relaxVec2 <- c(1,1,1,length(epsvec))
#' #In the first iteration only epsvec[1] is used, after that every assignment is exact.
#' #Not as fast as the previous version but usually no drawbacks at all in the computed solution.
#' #Time advantage more visible for large datasets.
#' 
#' relaxVec3 <- c(3,2,3,2)
#' #in the first 3 iterations epsvec[1],epsvec[1:2],epsvec[1:3] are used in the assignments,
#' #after that epsvec[2:x] is used, where x starts at 3 (=maximum(a,b)) and increases
#' #by 2 everytime the algorithm does not improve. When x >= length(epsvec) all assignments
#' #are done with epsvec[3:length(epsvec)].
#'
#' res1 <- kmeansbaryeps(epsvec, zeta, pplist_samecard, penalty=0.1, add_del=5, relaxVec = relaxVec1)
#' res2 <- kmeansbaryeps(epsvec, zeta, pplist_samecard, penalty=0.1, add_del=5, relaxVec = relaxVec2)
#' res3 <- kmeansbaryeps(epsvec, zeta, pplist_samecard, penalty=0.1, add_del=5, relaxVec = relaxVec3)
#' plot(res1$barycenter, add=TRUE, col="blue", pch=16) #adding the computed barycenter in blue
#'
#' @seealso \code{\link{kmeansbary}} for a similar function that works without epsilon relaxation
#'
#' @export
#'
kmeansbaryeps <- function(epsvec, zeta, pplist, penalty, add_del = Inf, surplus = 0,
                          relaxVec = c(20,1,1,1), N = 200L, eps = 0.005, verbose = 0) {
  stopifnot(class(zeta) == "ppp")
  stopifnot(any(class(pplist) == "list"))
  stopifnot(all(sapply(pplist, function(x) {class(x) == "ppp"})))
  #stopifnot(length(relaxVec) == 4)
  if (length(relaxVec) != 4 || any(relaxVec < c(1,1,1,1))) {
    relaxVec <- c(20,1,1,1)
  }
  #if (relaxVec[3] < relaxVec[2]) { does not make too much sense, but is not really my problem
  #  relaxVec <- c(20,1,1,20)
  #}
  relaxVec <- relaxVec - c(0,1,1,0)
  npplist <-sapply(pplist, npoints)
  nzeta <- npoints(zeta)
  n <- max(c(nzeta, npplist))
  k <- length(pplist)
  stopifnot(k > 1)
  
  zetax <- c(zeta$x, rep(NA, n-nzeta+surplus))
  zetay <- c(zeta$y, rep(NA, n-nzeta+surplus))
  
  ppmatx <- c(pplist[[1]]$x, rep(NA, n-npplist[1]+surplus))
  ppmaty <- c(pplist[[1]]$y, rep(NA, n-npplist[1]+surplus))
  for (i in 2:k) {
    ppmatx <- cbind(ppmatx, c(pplist[[i]]$x, rep(NA, n-npplist[i]+surplus)))
    ppmaty <- cbind(ppmaty, c(pplist[[i]]$y, rep(NA, n-npplist[i]+surplus)))
  }
  
  if (add_del > .Machine$integer.max) add_del <- .Machine$integer.max
  # this is in fact only to catch the Inf value (since it is an int in C++ not clear if
  # it would work with R_PosInf)
  res <- .Call(`_ttbary_kMeansBaryEps`, epsvec, zetax, zetay, ppmatx, ppmaty, penalty, add_del, relaxVec, N, eps, verbose)
  win <- bounding.box.xy(c(superimpose(pplist)$x,res$baryx),c(superimpose(pplist)$y,res$baryy))
  
  baryx <- na.omit(res$barycenterx)
  baryy <- na.omit(res$barycentery)
  if (!isTRUE(all.equal(attr(baryx, "na.action"), attr(baryy, "na.action")))) {
    warning("NA coordinates of barycenter did not agree")
    message("NA in x-coords: ", paste(which(is.na(baryx)), collapse=" "))
    message("NA in y-coords: ", paste(which(is.na(baryy)), collapse=" "))
  }
  return(list(cost=res$cost, barycenter=ppp(baryx, baryy, window=win), iterations=res$iterations))
}

sampleFromDatapp <- function(n, pplist) {
  win <- bounding.box.xy(superimpose(pplist))
  allx <- unname(unlist(lapply(pplist, function(ll) {ll$x})))
  ally <- unname(unlist(lapply(pplist, function(ll) {ll$y})))
  res <- sampleFromData(n, allx, ally)

  return(ppp(res$zetax, res$zetay, window=win))
}

#' Compute Pseudo-Barycenter of a List of Point Patterns on a Network
#'
#' Starting from an initial candidate point pattern \code{zeta}, use a k-means-like
#' algorithm to compute a local minimum in the barycenter problem based on the TT-1 metric
#' for a collection of point patterns on a network. The data needs to be in a special 
#' form which can be produced with the function \code{\link{netsplit}}.
#'
#' @param dpath a square matrix whose (\code{i},\code{j})th entry is the shortest-path distance between vertex \code{i}
#'              and vertex \code{j}. Vertex means either network vertex or data point.
#' @param zeta a vector containing the vertex-indices of the initial candidate for the barycenter.
#' @param ppmatrix a matrix specifying in its columns the vertex-indices of the different data point patterns. A virtual
#'           index that is one greater than the maximum vertex-index can be used to fill up columns so they all have
#'           the same length (see examples).
#' @param penalty the penalty for adding/deleting points when computing TT-1 distances.
#' @param N the maximum number of iterations.
#' @param eps the algorithm stops if the relative improvement of the objective function between two iterations is less
#'          than eps.
#'
#' @details  Given \eqn{k} planar point patterns \eqn{\xi_1, \ldots, \xi_k}{xi_1, ..., xi_k} (specified by giving the indices
#'           of their points in the \eqn{k} columns of \code{ppmatrix}), this function finds a local minimizer \eqn{\zeta^*}{zeta*} of
#'           \deqn{\sum_{j=1}^k \tau_1(\xi_j, \zeta),}{sum_{j=1}^k tau_1(xi_j, zeta),}
#'           where \eqn{\tau_1}{tau_1} denotes the TT-1 metric based on the shortest-path metric between points in the network.
#'           
#'           Starting from an initial candidate point pattern \code{zeta} (specified by giving the indices
#'           of its points), the algorithm alternates between assigning a point from each pattern \eqn{\xi_j}{xi_j}
#'           to each point of the candidate and computing new candidate patterns by shifting points (addition and deletion
#'           of points is currently not implemented).
#'           A detailed description of the algorithm is given in Müller, Schuhmacher and Mateu (2019).
#'                      
#'           The most convenient way to obtain objects \code{dpath} and \code{ppmatrix} of the right form is by calling
#'           \code{\link{netsplit}} and extracting components \code{network$dpath} and \code{ppmatrix} from the resulting
#'           object (see examples below). 
#'
#' @return A list containing the following components:
#'         \item{cost}{the sum of TT-1 distances between the computed pseudo-barycenter and the point patterns.}
#'         \item{barycenter}{the pseudo-barycenter as a vector of vertex-indices.}
#'         \item{zetalist}{a list containing the alternative vertex-indices for each point of the pseudo-barycenter.}
#'         \item{barycost}{a vector containing the cluster costs for each point of the pseudo-barycenter 
#'         (the alternative indices in \code{zetalist} lead to the same cluster cost).} 
#'         \item{perm}{the permutation matrix for the clusters.}
#'         \item{iterations}{the number of iterations required until convergence.}
#'
#' @references Raoul Müller, Dominic Schuhmacher and Jorge Mateu (2019).\cr
#'             Metrics and Barycenters for Point Pattern Data.\cr
#'             Preprint \href{https://arxiv.org/abs/1909.07266}{arXiv:1909.07266}
#'
#' @author Raoul Müller  \email{raoul.mueller@uni-goettingen.de}\cr
#'         Dominic Schuhmacher \email{schuhmacher@math.uni-goettingen.de}
#'
#' @examples 
#' set.seed(123456)
#' nvert <- 100 #number of vertices in the network
#' npp <- 5 #number of data point patterns
#' npts <- 40 #number of points per data point pattern
#' ln <- delaunayNetwork(runifpoint(nvert)) #create an artificial network
#' ppnetwork <- runiflpp(npts,ln,nsim = npp)
#'   #simulate npp point patterns with npts points each
#' 
#' plot(ppnetwork[[1]]$domain, cex=0.5, main="")
#' for (i in 1:npp) {
#'   plot(as.ppp(ppnetwork[[i]]),vpch=1,col=i,add=TRUE)
#'      #plotting the point patterns in different colors
#' }
#' 
#' res <- netsplit(ln, ppnetwork)
#'   #incorporate data point patterns into the network
#'   #calculating all pairwise distances between vertices
#'   #and creating matrix of vertex-indices of data point patterns
#'   
#' zeta <- sample(res$nvirtual - 1, median(res$dimensions))
#'   #sample random vertex-indices in the network
#'   #taking as cardinality the median of point pattern cardinalities
#' 
#' res2 <- kmeansbarynet(res$network$dpath, zeta, res$ppmatrix, penalty = 0.1)
#' 
#' barycenter <- ppp(res$network$vertices$x[res2$barycenter], res$network$vertices$y[res2$barycenter])
#'   #construct the barycenter pattern based on the index information in res2
#' points(barycenter,cex = 1.2, lwd = 2, pch = 4, col = "magenta")
#'   #add the computed barycenter as magenta crosses
#' 
#' res2$cost
#' #[1] 18.35171
#' sumppdistnet(res$network$dpath, res2$barycenter, res$ppmatrix, penalty=0.1, type="tt", p=1, q=1)
#' #[1] 18.35171
#' #attr(,"distances")
#' #[1] 3.666471 3.774709 3.950079 3.841166 3.119284
#'
#' @seealso \code{\link{kmeansbary}} for a similar function for point patterns in \eqn{R^2}
#'
#' @export
#'
kmeansbarynet <- function(dpath, zeta, ppmatrix, penalty, 
                             N = 200L, eps = 0.005) {
  nzeta <- length(zeta)
  nppmat <- dim(ppmatrix)[1]
  n <- max(nzeta,nppmat)
  k <- dim(ppmatrix)[2]
  bign <- length(dpath[1,])
  
  zeta <- c(zeta, rep(bign+1, n-nzeta))
  ppmat <- rbind(ppmatrix, matrix(bign+1, n-nppmat, k))
  
  zeta <- zeta-1 #indices now match the column/row indices in cpp
  ppmat <- ppmat-1 
  
  dpath[dpath > 2*penalty] <- 2*penalty
  dpath <- rbind(dpath,rep(penalty,bign))
  dpath <- cbind(dpath,c(rep(penalty,bign),0))
  
  res <- .Call(`_ttbary_kMeansBaryNet`, dpath, zeta, ppmat, penalty, add_del = 0, N, eps)
  
  zeta <- res$barycenter+1
  order <- order(zeta) #virtual points at the end
  zeta <- zeta[order]
  perm <- res$perm+1
  perm <- perm[order,]
  nzeta <- length(zeta) - sum(zeta == bign+1) #right now obsolete, but if zetapoints change from real to virtual and back this is necessary
  
  zetacost <- NULL #clustercosts for every barycenterpoint
  for (i in 1:length(zeta)) {
    zetacost[i] <- sum(dpath[zeta[i],perm[i,]])
  }
  
  zetalist <- vector("list",length(nzeta)) #alternative barycenterpoints for every non-virtual point
  for (i in 1:nzeta) {
    colsum <- rowSums(dpath[,perm[i,]])
    zetalist[[i]] <- seq(1,length(dpath[1,]))[colsum == min(colsum)]
  }
  
  
  return(list("cost" = res$cost, "barycenter" = zeta, "zetalist" = zetalist, "barycost" = zetacost, "perm" = perm, "iterations" = res$iterations))
}

#' Incorporate Point Patterns into a Network
#'
#' Given a network and a list of point patterns on this network, create a new network from all the
#' vertices of the original network plus all the points in the patterns, splitting any edges that
#' contain such points into several shorter edges. This function keeps track which vertex-indices 
#' represent each of the data point patterns. The returned object contains all the components
#' needed for a call to \code{\link{kmeansbarynet}}.
#'
#' @param network an object of class \code{linnet} or \code{lpp}. In the latter case the \code{domain}
#'                component is extracted and any points of the \code{lpp} are ignored.
#' @param pplist a list containing (at least) \code{x}- and \code{y}-coordinates of the point patterns,
#'               which will be projected onto the network
#'
#' @details  This function relies heavily on code from the package \code{spatstat} to create the
#'           new network and efficiently compute all pairwise shortest-path distances between
#'           the new vertices.
#'           
#'           If not all point patterns are of the same size, this function fills up the vertex-indices
#'           of the smaller patterns with a virtual index that is one larger than the maximal
#'           index appearing in the new network. This structure is required for 
#'           calling \code{\link{kmeansbarynet}}.
#'
#' @return A list containing the following components:
#'         \item{network}{the new network with all the points added as vertices. Contains also the matrix
#'         of shortest-path distances between all these points.}
#'         \item{ppmatrix}{a matrix containing the new vertex-indices of the data point patterns, one column corresponds to one point pattern.} 
#'         \item{dimensions}{a vector containing the cardinalities of the data point patterns.} 
#'         \item{nvirtual}{the index of the virtual point.}
#'
#' @author Raoul Müller  \email{raoul.mueller@uni-goettingen.de}\cr
#'         Dominic Schuhmacher \email{schuhmacher@math.uni-goettingen.de}
#'
#' @examples # See the example for kmeansbarynet.
#' 
#' @seealso \code{\link{kmeansbarynet}}
#'
#' @export
#'
netsplit <- function(network, pplist){
  
  dims <- NULL
  rawlist <- FALSE
  if (is.null(pplist[[1]]$data$x)) {
    rawlist <- TRUE
  }
  if (rawlist) {
    for (i in 1:length(pplist)){
      dims <- c(dims,length(pplist[[i]]$x))
    }
  }else{
    for (i in 1:length(pplist)){
      dims <- c(dims,length(pplist[[i]]$data$x))
    }
  }
  
  
  virtualvertex <- ifelse("linnet" %in% class(network),network$vertices$n,network$domain$vertices$n) + sum(dims) + 1 #index of the virtual vertex in the resulting network
  
  ppmatrix <- matrix(NA,max(dims),length(pplist))
  if ("linnet" %in% class(network)) {
    newnetwork <- network
  }else{
    newnetwork <- network$domain
  }
  index <- newnetwork$vertices$n+1 #index of first added point
  if (rawlist) {
    for (i in 1:length(pplist)){
      newnetwork <- insertVertices(newnetwork,cbind(pplist[[i]]$x,pplist[[i]]$y)) #add next point pattern to network
      ppmatrix[,i] <- c(seq(index,(index+dims[i]-1)),rep(virtualvertex,max(dims)-dims[i])) #indices of vertices are added to ppmatrix
      index <- index + dims[i]
    }
  }else{
    for (i in 1:length(pplist)){
      newnetwork <- insertVertices(newnetwork,cbind(pplist[[i]]$data$x,pplist[[i]]$data$y)) #add next point pattern to network
      ppmatrix[,i] <- c(seq(index,(index+dims[i]-1)),rep(virtualvertex,max(dims)-dims[i])) #indices of vertices are added to ppmatrix
      index <- index + dims[i]
    }
  }
  
  
  if (is.null(newnetwork$dpath)){
    newnetwork <- as.linnet(newnetwork,sparse=FALSE)
  }
  
  return(list("network" = newnetwork, "ppmatrix" = ppmatrix, "dimensions" = dims ,"nvirtual" = virtualvertex))
}

kmeansbary_mat <- function(zetax, zetay, ppmatx, ppmaty, penalty,
                           add_del = Inf, surplus = 0, N = 200L, eps = 0.005, verbose = 0) {
  nzeta <- length(zetax)
  stopifnot(nzeta == length(zetay))
  nppmatx <- dim(ppmatx)[1]
  stopifnot(nppmatx == dim(ppmaty)[1])
  n <- max(nzeta,nppmatx)
  k <- dim(ppmatx)[2]
  stopifnot(k == dim(ppmaty)[2])
  
  zetax <- c(zetax, rep(NA, n-nzeta+surplus))
  zetay <- c(zetay, rep(NA, n-nzeta+surplus))
  ppmatx <- rbind(ppmatx, matrix(NA, n-nppmatx+surplus, k))
  ppmaty <- rbind(ppmaty, matrix(NA, n-nppmatx+surplus, k))
  
  if (add_del > .Machine$integer.max) add_del <- .Machine$integer.max
  
  .Call(`_ttbary_kMeansBary`, zetax, zetay, ppmatx, ppmaty, penalty, add_del, N, eps, verbose)
}

kmeansbaryeps_mat <- function(epsvec, zetax, zetay, ppmatx, ppmaty, penalty,
                              add_del = Inf, surplus = 0, N = 200L, eps = 0.005, verbose = 0) {
  nzeta <- length(zetax)
  stopifnot(nzeta == length(zetay))
  nppmatx <- dim(ppmatx)[1]
  stopifnot(nppmatx == dim(ppmaty)[1])
  n <- max(nzeta,nppmatx)
  k <- dim(ppmatx)[2]
  stopifnot(k == dim(ppmaty)[2])
  
  zetax <- c(zetax, rep(NA, n-nzeta+surplus))
  zetay <- c(zetay, rep(NA, n-nzeta+surplus))
  ppmatx <- rbind(ppmatx, matrix(NA, n-nppmatx+surplus, k))
  ppmaty <- rbind(ppmaty, matrix(NA, n-nppmatx+surplus, k))
  
  if (add_del > .Machine$integer.max) add_del <- .Machine$integer.max
  
  .Call(`_ttbary_kMeansBaryEps`, epsvec, zetax, zetay, ppmatx, ppmaty, penalty, add_del, c(20,1,1,1), N, eps, verbose)
}
