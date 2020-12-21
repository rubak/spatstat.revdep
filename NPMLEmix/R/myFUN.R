#' @importClassesFrom Matrix CsparseMatrix
#' @import methods
### #' @import Rmosek
#' @import dplyr
#' @import FDRreg
#' @import Rmosek
#' @import R.utils
### #' @import spatstat
### #' @import latexpdf
### #' @import Matrix
#' @import Rcpp
#' @import REBayes
### #' @import CAMAN
#' @import progress
#' @import pbapply
### #' @import Hmisc
#' @import pracma
### #' @import mosaic
#' @importFrom methods as
#' @importFrom grDevices col2rgb rgb
#' @importFrom stats dnorm lm optim quantile rnorm runif sd
#' @importFrom Rmosek mosek

NULL
########### Packages required #############################
#' Likelihood based inference in mixture models and multiple hypotheses testing
#'
#' The NPMLEmix-package fits nonparametric Gaussian location mixture models for
#' z-scores arising out of several hypotheses, while taking into account any
#' available covariate information. It also provides three important functions:
#' marg1(), marg2() (both based on marginal likelihoods), npmleEM() (based on
#' joint data likelihood) for inference in multiple testing.
#'
#' @section Functions:
#' The principal functions in the NPMLEmix-package: marg1(), marg2(), npmleEM().
#'
#' @details \eqn{(Y_1,X_1),\ldots ,(X_n,Y_n)} are i.i.d. samples drawn from the model,
#' \deqn{Y|X=x\sim (1-\pi^*(x))\phi(y)+\pi^*(x)\underbrace{\int_{\theta}\phi(y-\theta)\,dG(\theta)}_{\phi_1(y)}, \qquad X\sim m_X(\cdot)}
#' where \eqn{\pi^*(\cdot)} represents the logistic link function, \eqn{\phi(\cdot)} is
#' the standard Gaussian density and \eqn{G(\cdot)} is some unknown probability measure on
#' the real line. Usually, \eqn{\pi^*(\cdot)} is referred to as the \emph{signal proportion},
#' \eqn{\phi_1(\cdot)} is called the \emph{signal density} and \eqn{G(\cdot)} is called \emph{mixing distribution}.
#' The \eqn{i^{th}} local false discovery rate is then defined as
#' \deqn{lfdr_i=\frac{(1-\pi^*(X_i))\phi_0(Y_i)}{(1-\pi^*(X_i))\phi_0(Y_i)+\pi^*(X_i)\phi_1(Y_i)}}
#' All the principal functions estimate the unknown parameters - \eqn{\pi^*(\cdot)} and
#' \eqn{\phi_1(\cdot)}, and consequently the \eqn{lfdr_i}'s.
#' The optimization algorithms use quasi-Newton routines such as the BFGS (Broyden-Fletcher-Goldfarb-Shanno)
#' algorithm and the separable convex optimization routine available in the Rmosek
#' optimization suite.
#' The principal functions accept a vector of z-scores (\eqn{Y})'s and a covariate matrix X
#' in their list of arguments. Read the documentations for each function to
#' check whether or not to add a column of \eqn{1}'s to X matrix.
#'
#' @references Deb, N., Saha, S., Guntuboyina, A. and Sen, B., 2018. Two-component Mixture Model in the Presence of Covariates. arXiv preprint arXiv:1810.07897.
#' @references Scott, J.G., Kelly, R.C., Smith, M.A., Zhou, P. and Kass, R.E., 2015. False discovery rate regression: an application to neural synchrony detection in primary visual cortex. Journal of the American Statistical Association, 110(510), pp.459-471.
#' @references Efron, B., 2005. Local false discovery rates.
#' @references Koenker, R. and Mizera, I., 2014. Convex optimization in R. Journal of Statistical Software, 60(5), pp.1-23.
#'
#' @docType package
#' @name NPMLEmix-package
NULL

suppressWarnings(require(latexpdf, quietly = TRUE))
# suppressWarnings(require(splines, quietly = TRUE))
suppressWarnings(require(methods, quietly = TRUE))
suppressWarnings(require(Matrix, quietly = TRUE))
suppressWarnings(require(Rmosek, quietly = TRUE))
suppressWarnings(require(mosaic, quietly = TRUE))
suppressWarnings(require(dplyr, quietly = TRUE))
suppressWarnings(require(spatstat, quietly = TRUE))
suppressWarnings(require(pracma, quietly=TRUE))
suppressWarnings(require(FDRreg, quietly = TRUE))
suppressWarnings(require(R.utils, quietly = TRUE))
suppressWarnings(require(Rcpp, quietly = TRUE))
suppressWarnings(require(REBayes, quietly = TRUE))
suppressWarnings(require(CAMAN, quietly = TRUE))
suppressWarnings(require(progress, quietly = TRUE))
suppressWarnings(require(pbapply, quietly = TRUE))
suppressWarnings(require(Hmisc, quietly = TRUE))

##########################################################
# expanding each column of a matrix via splines
##########################################################

# spline_expand = function(xx, df = 3){
#  require(splines)
#  b1 = bs(xx[,1], df = df)
#  b2 = bs(xx[,2], df = df)
#  model.matrix( ~  b1 + b2 - 1)
# }

##########################################################
# Calculating the non-null proportion
##########################################################

evnnprop=function(beta, xx){
  nxx=cbind(1,xx)
  return(apply(nxx, 1, function(x) return(as.numeric(crossprod(beta,x)))))
}

##############################################################
# main function for generating data
# from logistic pix and gaussian mixture f1
# sx denotes the link within logistic
# tdparams determines the mixing density underlying f1
# covariates are drawn from Unif([0,1]^2)
# third order spline expansion of covariates is also returned
##############################################################

#' Simulates data from the aforementioned model
#'
#' This function can be used to simulate observations from the aforementioned model, if
#' \eqn{G(\cdot)} is chosen as a finite Gaussian mixture. It returns the true local
#' false discovery rates which determine the optimal multiple testing procedure.
#' @param n Number of z-scores to be generated.
#' @param x \eqn{n\times}p data matrix. Do not add an additional column of \eqn{1's}.
#' @param sx The vector of coefficients for the logistic function. The first entry will be considered as the intercept term by default. Requires compatibility with x. See \strong{Details}.
#' @param atoms The vector of means for each component of the mixing distribution.
#' @param probs The probability vector for the mixing distribution.
#' @param variances The vector of variances for each component of the mixing distribution. Requires compatibility with atoms and probs. See \strong{Details}.
#' @details Given \eqn{X=x}, a Bernoulli\eqn{(\pi^*(x))} sample is drawn. If the outcome is
#' 1 (0), a z-score is drawn from \eqn{\phi_1(\cdot)} \eqn{(\phi(\cdot))}. All the observations
#' corresponding to a Bernoulli outcome 1 (0) are termed as \emph{non-null observations}
#' (\emph{null observations}).\cr
#' The length of sx should be 1 more than the number of columns of the data matrix x.\cr
#' The vectors - atoms, probs and variances must have the same length.
#' @return The output is a list with the following entries:
#' @return \item{y}{The vector of simulated z-scores.}
#' @return \item{x}{The input data matrix.}
#' @return \item{pix}{The vector of signal proportions.}
#' @return \item{f0y}{The vector of standard Gaussian densities evaluated at simulated z-scores.}
#' @return \item{f1y}{The vector of signal densities evaluated at simulated z-scores.}
#' @return \item{den}{The vector of conditional densities evaluated at simulated z-scores.}
#' @return \item{localfdr}{The vector of local false discovery rates evaluated at simulated z-scores. Note that the local FDR can be interpreted as one minus the posterior probability that a given observation is non-null.}
#' @return \item{ll}{The average conditional log-likelihood.}
#' @return \item{nnind}{The indices corresponding to non-null observations.}
#' @export
#' @examples
#'
#' x=cbind(runif(1000),runif(1000))
#' n=1000
#' atoms=c(-2,0,2)
#' probs=c(0.48,0.04,0.48)
#' variances=c(1,16,1)
#' sx=c(-3,1.5,1.5)
#' ### Generating the data ###
#' st=makedata(n,x,sx,atoms,probs,variances)
#' ### Output the vector of local false discovery rates ###
#' st$localfdr
#'
#' @references Basu, P., Cai, T.T., Das, K. and Sun, W., 2018. Weighted false discovery rate control in large-scale multiple testing. Journal of the American Statistical Association, 113(523), pp.1172-1183.
#' @references Scott, J.G., Kelly, R.C., Smith, M.A., Zhou, P. and Kass, R.E., 2015. False discovery rate regression: an application to neural synchrony detection in primary visual cortex. Journal of the American Statistical Association, 110(510), pp.459-471.

makedata = function(n, x , sx, atoms, probs, variances){

  # constructing f1
  atoms = atoms; probs = probs; variances = variances;
  # checking validity of tdparams
  stopifnot(length(unique(c(length(atoms), length(probs), length(variances)))) == 1)
  # number of components
  nc = length(atoms)
  # f1 function
  f1 = make_density(atoms, probs, variances)

  # here we generate data
  ####################################################################
  # starting data generation process
  xs = x
  # checking further validity
  stopifnot(length(unique(c(ncol(xs), (length(sx)-1)))) == 1)
  stopifnot(length(unique(c(nrow(xs), n))) == 1)
  # generating pi(.)
  pix = sapply(evnnprop(sx,xs), function(x) return(1/(1+exp(-x))))

  # generating observations
  nulls = (runif(n) >= pix);  nind = which(nulls); nnind = setdiff(1:n, nind);
  # populating data
  theta = array(0, dim = n)
  wc = sample(nc, length(nnind), replace = TRUE, prob = probs)
  theta[nnind] = atoms[wc] + rnorm(length(nnind)) * sqrt(variances[wc])
  y = rnorm(n) + theta
  # completed data generation process
  #####################################################################

  f0y = dnorm(y)
  f1y = f1(y)
  den = (1-pix)*f0y + pix*f1y
  localfdr = (1-pix)*f0y/den
  ll = mean(log(den))

  return(list(y = y, x = x, xs = xs,
              pix = pix,
              f0y = f0y, f1y = f1y, localfdr = localfdr, ll = ll, f1 = f1,
              nnind = nnind, den = den))
}

############################################################
# utility to make a density given atoms, probs and vars
############################################################

make_density = function(atoms,probs,variances){
  f = function(xx) sum(sapply(1:length(atoms), function(j)
    probs[j] * dnorm(xx, atoms[j], sqrt(1 + variances[j]))))
  return(Vectorize(f))
}

############################################################
# Rejection set in the multiple hypotheses testing problem
############################################################

#' Finds the rejection set in a multiple testing problem
#'
#' This function accepts a vector of local false discovery rates from a family of hypotheses
#' and a level parameter, to compute the rejection set.
#'
#' @param locfdr The vector of local false discovery rates (actual or estimated) corresponding to a family of hypotheses.
#' @param level The level at which the false discovery rate is to be controlled. Should ideally be a scalar in \eqn{[0,1]}.
#' @details The problem of optimal inference in multiple hypotheses testing has been widely studied in literature. In particular, this function adopts the framework and algorithm proposed in Basu et al. See \strong{References}.
#' @return A vector of 1s and 0s with 1s indicating the hypotheses which are to be rejected.
#' @export
#' @examples
#'
#' x=cbind(runif(1000),runif(1000))
#' n=1000
#' atoms=c(-2,0,2)
#' probs=c(0.48,0.04,0.48)
#' variances=c(1,16,1)
#' sx=c(-3,1.5,1.5)
#' stdata=makedata(n,x,sx,atoms,probs,variances)
#' ### Obtain the rejection set ###
#' reject=reject_set(stdata$lo)
#' @references Basu, P., Cai, T.T., Das, K. and Sun, W., 2018. Weighted false discovery rate control in large-scale multiple testing. Journal of the American Statistical Association, 113(523), pp.1172-1183.
#' @references Deb, N., Saha, S., Guntuboyina, A. and Sen, B., 2018. Two-component Mixture Model in the Presence of Covariates. arXiv preprint arXiv:1810.07897.

reject_set = function(locfdr, level=0.1){
  asclfdr = sort(locfdr)
  cumlfdr=cummean(asclfdr)
  ind=max(which(cumlfdr <= level))
  rejvec=as.numeric(locfdr <= asclfdr[ind])
  return(rejvec)
}

# this function is used by REBayes to construct the mids
histbin = function(x, m, eps = 1e-06, weights = weights) {
  u = seq(min(x) - eps, max(x) + eps, length = m)
  xu = findInterval(x, u)
  txu = tabulate(xu)
  midu = (u[-1] + u[-m])/2
  wnz = (txu > 0)
  if (length(weights)) {
    if (length(weights) == length(x))
      w = as.numeric(tapply(weights, xu, sum))
    else stop("length(weights) not equal length(x)")
  }
  else w = txu[wnz]/sum(txu[wnz])
  list(midu = midu[wnz], w = w)
}

# let's try implementing primal formulation of kw in mosek
# this can be used to solve both the weighted and pi-constrained problems
kwprimal = function(y, pivec = 1, weights = 1, grid_len = 100){
  LIMIT_THREAD = TRUE
  m = grid_len
  n = length(y)

  # if (length(pivec) == 1) pivec = rep(pivec, n)
  # if (length(weights) == 1) pivec = rep(weights, n)

  atoms = seq(from = min(y), to = max(y), length.out = m)
  fmat = dnorm(outer(y, atoms, "-"))

  # setting up problem
  # want to maximize so sense = max
  # the variables which need to be optimized are (p_j, v_i)
  moseko = list(sense = "max")
  # the problem is stated by mosek in the form:
  # sum of non-linear functions separable in arguments + linear part (c^T\pi)
  # for us there is no linear part so c = 0
  moseko$c = rep(0, m + n)
  # monotone constraints in pi
  A = rbind(cbind(fmat, -diag(n)), c(rep(c(1,0), times = c(m,n))))
  # sparsify to increase speed
  moseko$A = as(A, "CsparseMatrix")
  # moseko$A = Matrix(A, sparse = TRUE)
  # matrix of two rows
  # first row (blc) is lower bound on each linear relationship defined by A
  # second row (blu) is upper bound on each linear relationship defined by A
  moseko$bc = rbind(blc = c(rep(0, n), 1), buc = c(rep(0, n), 1))
  # box constraints: individual constraints on each pi_i
  # it seems to be easier to specify these separately
  # similar to above
  # first row (bxl) is lower bound on each pi_i
  # second row (bxu) is upper bound on each pi_i
  moseko$bx = rbind(blx = rep(0, m + n), bux = c(rep(1, m), rep(Inf, n)))
  # this is the interesting part, specifies the non-linear function in the objective
  # the function needs to separable in pi_i,
  # \sum_k f_k phi_k( g_k \pi_{j_k} + h_k)
  # where k is the k-th column of the following matrix
  # the number of columns is the number of non-linear functions the objective is separated into
  # the k-th column is a function of the j_k variable
  # the non-linearity is determined by 'type' (written as phi above)
  # you can multiply phi_k by some constant: f_k
  # also phi_k can be evaluated at a linear function of the variable in question: (g_k, h_k)
  opro = matrix(list(), nrow = 5, ncol = n)
  rownames(opro) = c("type", "j", "f", "g", "h")
  opro[1, ] = rep("LOG", n)
  opro[2, ] = m + 1:n
  opro[3, ] = weights # coefficients outside log, in this case all 1
  opro[4, ] = pivec # coefficient of v_i inside log
  opro[5, ] = (1-pivec) * dnorm(y) #constants inside log
  moseko$scopt = list(opro = opro)
  if(LIMIT_THREAD) {moseko$iparam$num_threads = 1}


  ans = mosek(moseko, opts = list(verbose = 1))

  probs = ans$sol$itr$xx[1:m]
  f1y = ans$sol$itr$xx[m + 1:n]

  return(list(atoms = atoms, probs = probs, f1y = f1y, ll = mean(log(f1y))))
}

# kwprimal using weights and histogram
kwprimal_weights = function(y, weights = 1, num_atoms = 100, hist_flag = TRUE, num_hist = num_atoms){

  n = length(y);
  if (length(weights) == 1) {weights = rep(weights, n)}
  # constructing histogram object
  if(hist_flag){
    histo = histbin(y, m = num_hist, weights = weights)
    # the compressed dataset is midu and w
    # print(names(histo))
    yy = histo$midu; ww = histo$w;
  } else {
    yy = y; ww = weights;
  }

  kwo = kwprimal(yy, weights = ww, grid_len = num_atoms);

  if (hist_flag){
    fmat = dnorm(outer(y, kwo$atoms, "-"))
    kwo$f1y = as.vector(fmat %*% kwo$probs)
    kwo$ll = mean(log(kwo$f1y))
  }

  return(kwo)
}

lgstep = function(y, x, w, b, num_atoms, blambda, level){
  f0y = dnorm(y)
  # take one step of EM on data
  kwo = kwprimal_weights(y, weights = w, num_atoms)
  lp = lregoptim(f0y, kwo$f1y, x, b, lambda = blambda)
  # compute results
  den = (1 - lp$p) * f0y + lp$p * kwo$f1y; ll = mean(log(den)); localfdr = (1 - lp$p) * f0y/(den);
  rejset = reject_set(localfdr, level)
  return(list(p = lp$p, b = lp$b, f1y = kwo$f1y, kwo=kwo,
              localfdr = localfdr, den = den, ll = ll, rejset = rejset))
}

# solve logistic problem given alternate density (AM version)
lregoptim = function(f0y, f1y, x, binit, lambda = 1e-2/length(f0y)){
  # by default this provides a l2-regularization
  # equivalent to gaussian prior N(0,100) on all parameters
  # lambda = c(1e-2/length(f0y), rep(lambda, length(binit)-1))

  # defining objective function
  llobj = function(bb, f0y, f1y, x, lambda){
    pivec = as.vector(1/(1 + exp(-x %*% bb)))
    - mean(log(pivec*f1y + (1-pivec)*f0y)) + sum(lambda*bb^2)
  }

  # defining first gradient
  llfirst = function(bb, f0y, f1y, x, lambda){
    pivec = as.vector(1/(1 + exp(-x %*% bb)))
    wts = (f1y - f0y)/(pivec*f1y + (1-pivec)*f0y) * pivec * (1-pivec)
    - apply(x, 2, function(vec) mean(vec * wts)) + 2*lambda*bb
  }

  optimres = optim(par = binit, fn = llobj, gr = llfirst,
                   f0y = f0y, f1y = f1y, x = x, lambda = lambda,
                   method = 'BFGS')

  b = optimres$par
  p = as.vector(1/(1 + exp(-x %*% b)))

  return(list(b = b, p = p))
}

###############################################################
# Implementing the marg1() method
###############################################################

#' Implements a profile likelihood based algorithm for estimating signal proportion and density
#'
#' This function estimates the signal proportion and the signal density by using
#' the marginal distribution of \eqn{Y}, followed by a profile likelihood based approach.
#' It returns the vector of estimated local false discovery rates and the corresponding
#' rejection set at a prespecified level for the false discovery rate.
#' @param y The observed vector of z-scores.
#' @param x The \eqn{n\times p} data matrix, where \eqn{n} must be equal to thelength of y. If you are interested in the intercept, you must add a column of \eqn{1's} to \eqn{x}.
#' @param blambda The tolerance threshold while implementing a quasi-Newton approach for estimating the signal proportion. Default is set to \eqn{1e-6/}length(y). We recommend not changing it unless absolutely sure.
#' @param level The level at which the false discovery rate is to be controlled. Should be a scalar in \eqn{[0,1]}. Default set to \eqn{0.05}.
#' @details Note that the marginal distribution of \eqn{Y} based on the aforementioned model
#' is same as that in a standard two-groups model (Efron 2008, see \strong{References}). Fixing
#' \eqn{\bar\pi = \mathbf{E}[\pi(X)]}, the signal density \eqn{\phi_1(\cdot)} is
#' estimated using the Rmosek optimization suite. The primary idea is to approximate
#' the mixing distribution \eqn{G(\cdot)} using \eqn{\max\{100,\sqrt{n}\}} many
#' components, each having a suitable Gaussian distribution. The signal proportion
#' is then estimated
#' using the BFGS algorithm. Finally, the algorithm chooses the best value of \eqn{\bar\pi} based
#' on a profile likelihood approach.
#' @return This function returns a list consisting of the following:
#' @return \item{p}{The estimated prior probabilities, i.e., \eqn{\hat\pi(\cdot)} evaluated at the data points.}
#' @return \item{b}{The estimates for the coefficient vector in the logistic function.}
#' @return \item{f1y}{The vector of estimated signal density evaluated at the data points.}
#' @return \item{kwo}{This is a list with four items - i. \emph{atoms}: The vector of means for the Gaussian distributions used to approximate \eqn{G(\cdot)}, ii. \emph{probs}: The vector of probabilities for each Gaussian component used to approximate \eqn{G(\cdot)}, iii. \emph{f1y}: Same as f1y above, iv. \emph{ll}: The average of the logarithmic values of f1y. }
#' @return \item{localfdr}{The vector of estimated local false discovery rates evaluated at the data points.}
#' @return \item{den}{The vector of estimated conditional densities evaluated at the data points.}
#' @return \item{ll}{The log-likelihood evaluated at the estimated optima.}
#' @return \item{rejset}{The vector of \eqn{1}s and \eqn{0}s where \eqn{1} indicates that the corresponding hypothesis is to be rejected.}
#' @return \item{pi0}{The average of the entries of the vector \emph{p}.}
#' @return \item{ll_list}{The vector of profile log-likelihoods corresponding to a pre-determined set of grid points for \eqn{\bar\pi}. The highest element of this vector is the output in \emph{ll}.}
#' @references Deb, N., Saha, S., Guntuboyina, A. and Sen, B., 2018. Two-component Mixture Model in the Presence of Covariates. arXiv preprint arXiv:1810.07897.
#' @references Koenker, R. and Mizera, I., 2014. Convex optimization, shape constraints, compound decisions, and empirical Bayes rules. Journal of the American Statistical Association, 109(506), pp.674-685.
#' @export
#'
#' @references Efron, B., 2008. Microarrays, empirical Bayes and the two-groups model. Statistical science, pp.1-22.

marg1 = function(y, x, blambda = 1e-6/length(y), level = 0.05){
  # x=cbind(1, datax)
  mt = abs(y) - mean(abs(rnorm(1e4)));
  pi0grid = seq(from = 0.01, to = 0.49, by = 0.01);
  verbose=FALSE;
  grid_len = max(100, round(sqrt(length(y))));
  n = length(y); f0y = dnorm(y);
  kwm = kwprimal_weights(y, num_atoms = grid_len);
  fmy = kwm$f1y;

  m1step = function (pt) {
    weights = pmax(pmin(1 - (1-pt)*f0y/fmy, 0.99), 0.01)
    muinit = mean(mt)/pt
    binit = lm(mosaic::logit(pmax(pmin(mt/muinit, 0.9), 0.1)) ~ 0 + x)$coefficients;
    robj = lgstep(y, x, weights, binit, grid_len, blambda, level)
    robj$pi0 = pt
    return(robj)
  }

  if (verbose) {res = pblapply(pi0grid, m1step)}
  else {res = lapply(pi0grid, m1step)}

  ll_list = sapply(res, function (robj) robj$ll)
  bi = which.max(ll_list)
  robj = res[[bi]]
  robj$ll_list = ll_list

  return(robj)
}

# solve logistic problem given alternate density (EM version)
lregem = function(weights, x, binit, lambda = 1e-2/length(weights)){
  # by default this provides a l2-regularization
  # equivalent to gaussian prior N(0,100) on all parameters
  # lambda = c(1e-2/length(f0y), rep(lambda, length(binit)-1))

  # defining objective function
  llobj = function(bb, weights, x, lambda){
    pivec = as.vector(1/(1 + exp(-x %*% bb)))
    - mean(weights * log(pivec) + (1-weights) * log(1-pivec)) + sum(lambda*bb^2)
  }

  # defining first gradient
  llfirst = function(bb, weights, x, lambda){
    pivec = as.vector(1/(1 + exp(-x %*% bb)))
    tt = pivec * (1-pivec) * (weights/pivec - (1-weights)/(1-pivec))
    - apply(x, 2, function(vec) mean(vec * tt)) + 2*lambda*bb
  }

  optimres = optim(par = binit, fn = llobj, gr = llfirst,
                   weights = weights, x = x, lambda = lambda,
                   method = 'BFGS')

  b = optimres$par
  p = as.vector(1/(1 + exp(-x %*% b)))

  return(list(b = b, p = p))
}
# EM
lgem = function(y, x,
                weights, binit,
                grid_len = 3*max(100, round(sqrt(length(y)))), histFlag = TRUE, timed = 60,
                maxit = 600, tol = 1e-6,
                blambda = 1e-6/length(y), level){

  st = proc.time()['elapsed']

  n = length(y); f0y = dnorm(y);

  # initial values, may need to revisit this
  lp = NULL; lp$b = binit; lp$p = as.vector(1/(1 + exp(-x %*% lp$b)));

  # starting iterations
  convFlag = 1; ll = NULL; lp_list = NULL; kw_list = NULL; time_list = NULL; itcount = 1; err_list = NULL;

  # if(verbose) pb = txtProgressBar(max = maxit+1, style = 3)
  while(convFlag){
    kwo = kwprimal_weights(y, weights = weights, num_atoms = grid_len, hist_flag = histFlag)
    lp = lregem(weights, x, lp$b, lambda = blambda)

    weights_new = lp$p*kwo$f1y/((1-lp$p)*f0y + lp$p*kwo$f1y)
    iter_error = rmse(weights_new, weights); err_list = append(err_list ,iter_error);

    lp_list = append(lp_list, list(lp)); kw_list = append(kw_list, list(kwo));
    ll = append(ll, mean(log((1-lp$p)*f0y + lp$p*kwo$f1y)))

    ct = proc.time()['elapsed'] - st
    time_list = append(time_list, ct)
    convFlag = (ct <= timed) & (itcount < maxit) & (iter_error > tol)
    weights = weights_new
    itcount = itcount + 1

  }

  localfdr = 1 - weights
  rejset = reject_set(localfdr, level)
  den = (1-lp$p)*f0y + lp$p*kwo$f1y

  return(list(atoms = kwo$atoms, probs = kwo$probs, f1y = kwo$f1y,
              b = lp$b, p = lp$p,
              f0y = f0y, den = den,
              localfdr = localfdr, rejset = rejset,
              ll = ll,  lp_list = lp_list, kw_list = kw_list,
              err_list = err_list, time_list = time_list,
              runtime = proc.time()['elapsed'] - st))
}


# marginal method 2, using a grid of potential values of pi0
# remember to use a good choice of mt (marginal transform)


###############################################################
# Implementing the marg2() method
###############################################################

#' Implements a non-linear least squares based algorithm for estimating signal proportion and density
#'
#' This function estimates the signal proportion and the signal density by using
#' the conditional mean \eqn{Y|X=x}, followed by a non-linear least squares regression based approach.
#' It returns the vector of estimated local false discovery rates and the corresponding
#' rejection set at a prespecified level for the false discovery rate.
#' @param y The observed vector of z-scores.
#' @param x The \eqn{n\times p} data matrix, where \eqn{n} must be equal to thelength of y. If you are interested in the intercept, you must add a column of \eqn{1's} to \eqn{x}.
#' @param nlslambda The tolerance threshold while implementing a quasi-Newton approach for the non-linear least squares problem. Default is set to \eqn{1e-6/}length(y). We recommend not changing it unless absolutely sure.
#' @param level The level at which the false discovery rate is to be controlled. Should be a scalar in \eqn{[0,1]}. Default set to \eqn{0.05}.
#' @details Note that the conditional mean of \eqn{Y|X} based on the aforementioned model
#' is a non-linear function of the parameters, i.e., the logistic coefficients
#' and the mean of the marginal distribution of \eqn{Y}, \eqn{\mu^* = \mathbf{E}[Y]}.
#' This is a non-convex optimization problem in the parameters and is solved by varying
#' \eqn{\mu^*} over a predetermined grid, and optimizing over the logistic coefficients.
#' This is the estimate of \eqn{\pi^*(\cdot)} from the marg2() method. The estimate of
#' \eqn{\phi_1(\cdot)} is obtained as in the marg1() method by using the Rmosek optimization
#' suite, and the same discrete approximation to the mixing distribution \eqn{G(\cdot)}.
#' @return This function returns a list consisting of the following:
#' @return \item{p}{The estimated prior probabilities, i.e., \eqn{\hat\pi(\cdot)} evaluated at the data points.}
#' @return \item{b}{The estimates for the coefficient vector in the logistic function.}
#' @return \item{f1y}{The vector of estimated signal densities evaluated at the data points.}
#' @return \item{kwo}{This is a list with four items - i. \emph{atoms}: The vector of means for the Gaussian distributions used to approximate \eqn{G(\cdot)}, ii. \emph{probs}: The vector of probabilities for each Gaussian component used to approximate \eqn{G(\cdot)}, iii. \emph{f1y}: Same as f1y above, iv. \emph{ll}: The average of the logarithmic values of f1y. }
#' @return \item{localfdr}{The vector of estimated local false discovery rates evaluated at the data points.}
#' @return \item{den}{The vector of estimated conditional densities evaluated at the data points.}
#' @return \item{ll}{The log-likelihood evaluated at the estimated optima.}
#' @return \item{rejset}{The vector of \eqn{1}s and \eqn{0}s where \eqn{1} indicates that the corresponding hypothesis is to be rejected.}
#' @return \item{pi0}{The average of the entries of the vector \emph{p}.}
#' @return \item{ll_list}{The vector of profile log-likelihoods corresponding to a pre-determined set of grid points for \eqn{\mu^*}. The highest element of this vector is the output in \emph{ll}.}
#' @references Deb, N., Saha, S., Guntuboyina, A. and Sen, B., 2018. Two-component Mixture Model in the Presence of Covariates. arXiv preprint arXiv:1810.07897.
#' @references Koenker, R. and Mizera, I., 2014. Convex optimization, shape constraints, compound decisions, and empirical Bayes rules. Journal of the American Statistical Association, 109(506), pp.674-685.
#' @export
#'

marg2 = function(y, x, nlslambda = 1e-6/length(y), level = 0.05){
  # for an user given sequence pi0 (overall non-null proportion)
  # solve non-linear least squares on marginal_transform to get initial value of beta
  # x=cbind(1,x);
  pi0grid = seq(from = 0.01, to = 0.49, by = 0.01);
  mt = abs(y) - mean(abs(rnorm(1e4)));
  verbose = FALSE;
  solvef1 = TRUE;
  n = length(y); f0y = dnorm(y);
  mugrid = mean(mt)/pi0grid;

  computefn = function (mu) {
    # this function uses non-linear least squares to solve for beta for a fixed value of mu
    binit = lm(mosaic::logit(pmax(pmin(mt/mu, 0.9), 0.1)) ~ 0 + x)$coefficients
    nlso = m2boptim(mt/mu, x, binit, nlslambda)
    er = mean((mt - mu * nlso$p)^2)
    return(list(nlso = nlso, er = er, mu = mu, pi0 = mean(mt)/mu))
  }

  if (verbose) {
    res = pblapply(mugrid, computefn)
  } else {
    res = lapply(mugrid, computefn)
  }

  ers = sapply(res, function (robj) robj$er)
  bi = which.min(ers)
  regres = res[[bi]]


  if(solvef1){
    pi0grid = pi0grid[bi]
    #robj = m1(y, x, pi0grid = pi0grid[bi], mt = mt)
    robj = newmarg1(y, x, level = level)
  } else {
    robj = list(pi0grid = pi0grid, mugrid = mugrid,
                pi0 = pi0grid[bi], mu = mugrid[bi],
                regres = regres, res = res)
  }

  return(robj)
}

# solve marginal regression (logistic model) assuming mean under alternate is 1
m2boptim = function(y, x, binit, lambda = 1e-2/length(y)){
  # by default this provides a l2-regularization
  # equivalent to gaussian prior N(0,100) on all parameters
  # x = cbind(1, x);
  # defining objective function
  lsobj = function(bb, y, x, lambda){
    pivec = as.vector(1/(1 + exp(-x %*% bb)))
    mean((y - pivec)^2) + lambda*sum(bb^2)
  }

  # defining first gradient
  lsfirst = function(bb, y, x, lambda){
    pivec = as.vector(1/(1 + exp(-x %*% bb)))
    wts = -2 * (y - pivec) * pivec * (1-pivec)
    apply(x, 2, function(vec) mean(vec * wts)) + 2*lambda*bb
  }

  optimres = optim(par = binit, fn = lsobj, gr = lsfirst,
                   y = y, x = x, lambda = lambda,
                   method = 'BFGS')

  b = optimres$par
  p = as.vector(1/(1 + exp(-x %*% b)))

  return(list(b = b, p = p))
}

rmse = function(vec, truth) sqrt(mean((vec - truth)^2))
deciles = function(vec) quantile(vec, seq(from = 0.1, to = 0.9, by = 0.1))

rmsepath = function(obj, dd){
  rmses = numeric(length(obj$ll))
  with(obj , {
    for (i in 1:length(rmses)){
      p = obj$lp_list[[i]]$p
      f1 = obj$kw_list[[i]]$f1y
      lfdr = 1 - (p * f1)/((1-p) * f0y + p * f1)
      rmses[i] = rmse(dd$localfdr, lfdr)
    }
  })
  return(rmses)
}

# from a full mle solution with lots of iterations, extract an earlier iteration for comparison
extract_index = function(obj, ii){
  with(obj, {
    rr = list(atoms = kw_list[[ii]]$atoms, probs = kw_list[[ii]]$probs,
                f1y = kw_list[[ii]]$f1y, f0y = f0y,
                b = lp_list[[ii]]$b, p = lp_list[[ii]]$p,
                ll = ll[ii], rejset = rejset)
  })
  rr$den = with(rr, (1-p)*f0y + p * f1y)
  rr$localfdr = with(rr, (1-p)*f0y/((1-p)*f0y + p * f1y))
  return(rr)
}

add_alpha = function(cols, alpha = 0.7)
  rgb(red = t(col2rgb(cols))/255, alpha = alpha)

rmlast = function(ss) substr(ss, 1, nchar(ss) - 1)

# find good indices
fgi = function(vec, alpha = 0.5) which(vec >= max(vec) - alpha * sd(vec))

# find index upto which vec increases
fc = function(vec, ivec = 1:length(vec)){
  first_decrease = which(sign(c(diff(vec), -1)) == -1)[1]
  ivec[first_decrease]
}

# relocate vec so that minimum is zero
rmmin = function(vec) {vec - min(vec)}

# mod fdrreg
modf = function(obj){
  obj$b = obj$model$coef
  obj$p = obj$priorprob
  obj$f0y = obj$M0
  obj$f1y = obj$M1
  obj$den = (1-obj$p) * obj$f0y + obj$p * obj$f1y
  obj$ll = mean(log(obj$den))
  return(obj)
}

# extract initialization
extract_init = function(obj){
  return(list(p = obj$p, b = obj$b, w = 1 - obj$localfdr))
}


###############################################################
# Implementing the npmleEM() method
###############################################################

#' Implements the full likelihood approach based on the EM algorithm for estimating signal proportion and density
#'
#' This function estimates the signal proportion and the signal density by using
#' the full likelihood of the sample, followed by an EM algorithm based approach.
#' It returns the vector of estimated local false discovery rates and the corresponding
#' rejection set at a prespecified level for the false discovery rate.
#' @param y The observed vector of z-scores.
#' @param x The \eqn{n\times p} data matrix, where \eqn{n} mist be equal to thelength of y. If you are interested in the intercept, you must add a column of \eqn{1's} to \eqn{x}.
#' @param level The level at which the false discovery rate is to be controlled. Should be a scalar in \eqn{[0,1]}. Default set to \eqn{0.05}.
#' @param initp The initialization method for the EM algorithm. It should be either \eqn{1,2,3} or \eqn{4}. \emph{1} indicates a marg1() initialization, \emph{2} indicates a marg2() initialization, \emph{3} indicates a FDRreg() initialization (see \strong{Details} and \strong{References}) and \emph{4} chooses that initialization among marg1(), marg2() and FDRreg() which yields the highest sample likelihood. Default is set to \emph{1}.
#' @details The key observation in the full likelihood approach is that the M-step of
#' the EM algorithm results in two decoupled optimization problems, one involving \eqn{\pi^*(\cdot)}
#' and the other involving \eqn{\phi_1(\cdot)}. These two individual problems are then
#' solved using the BFGS algorithm and the Rmosek optimization suite, as has been discussed previously in
#' the \strong{Details} sections of the methods marg1() and marg2().\cr
#' The FDRreg() method was introduced in Scott et al (2015). We recommend using the version of
#' the FDRreg() package available in \url{https://github.com/jgscott/FDRreg/tree/master/R_pkg}.
#' @return This function returns a list consisting of the following:
#' @return \item{atoms}{The vector of means for the Gaussian distributions used to approximate \eqn{G(\cdot)}.}
#' @return \item{probs}{The vector of probabilities for each Gaussian component used to approximate \eqn{G(\cdot)}.}
#' @return \item{f1y}{The vector of estimated signal densities evaluated at the data points.}
#' @return \item{f0y}{The vector of null densities evaluated at the data points.}
#' @return \item{b}{The estimates for the coefficient vector in the logistic function.}
#' @return \item{p}{The estimated prior probabilities, i.e., \eqn{\hat\pi(\cdot)} evaluated at the data points.}
#' @return \item{ll}{The log-likelihood evaluated at the estimated optima.}
#' @return \item{rejset}{The vector of \eqn{1}s and \eqn{0}s where \eqn{1} indicates that the corresponding hypothesis is to be rejected.}
#' @return \item{den}{The vector of estimated conditional densities evaluated at the data points.}
#' @return \item{localfdr}{The vector of estimated local false discovery rates evaluated at the data points.}
#' @references Deb, N., Saha, S., Guntuboyina, A. and Sen, B., 2018. Two-component Mixture Model in the Presence of Covariates. arXiv preprint arXiv:1810.07897.
#' @references Scott, J.G., Kelly, R.C., Smith, M.A., Zhou, P. and Kass, R.E., 2015. False discovery rate regression: an application to neural synchrony detection in primary visual cortex. Journal of the American Statistical Association, 110(510), pp.459-471.
#' @export
#'

npmleEM = function(y, x, level = 0.05, initp = 1){
  if(initp == 1)
  {
    m1n_ = marg1(y, x, level = level)
    init_list = list(m1n_ = m1n_)
    init_bi = which.max(sapply(init_list, function (ro) ro$ll))
    init_best = extract_init(init_list[[init_bi]])
    init_best_name = names(init_bi)
  }
  if(initp == 2)
  {
    m2n_ = marg2(y, x, level = level)
    init_list = list(m2n_ = m2n_)
    init_bi = which.max(sapply(init_list, function (ro) ro$ll))
    init_best = extract_init(init_list[[init_bi]])
    init_best_name = names(init_bi)
  }
  if(initp == 3)
  {
    ff=FDRreg(y,x[,-1])
    f_ = modf(ff)
    init_list = list(f_ = f_)
    init_bi = which.max(sapply(init_list, function (ro) ro$ll))
    init_best = extract_init(init_list[[init_bi]])
    init_best_name = names(init_bi)
  }
  if(initp == 4)
  {
    ff=FDRreg(y,x[,-1])
    f_ = modf(ff)
    m1n_ = marg1(y, x, level = level)
    m2n_ = marg2(y, x, level = level)
    init_list = list(f_ = f_, m1n_ = m1n_, m2n_ = m2n_)
    init_bi = which.max(sapply(init_list, function (ro) ro$ll))
    init_best = extract_init(init_list[[init_bi]])
    init_best_name = names(init_bi)
  }
  # EM starting from here, run at most 500 iterations
  em_ = lgem(y, x, weights = init_best$w, binit = init_best$b, timed = Inf, maxit = 100, level = level)
  em_ = extract_index(em_, length(em_$time_list))
  return(em_)
}

#######################################################
# Similar to marginal 1
#######################################################

newmarg1 = function(y, x, blambda = 1e-6/length(y), level = 0.05){
  # x=cbind(1, datax)
  mt = abs(y) - mean(abs(rnorm(1e4)));
  pi0grid = seq(from = 0.01, to = 0.49, by = 0.01);
  verbose=FALSE;
  grid_len = max(100, round(sqrt(length(y))));
  n = length(y); f0y = dnorm(y);
  kwm = kwprimal_weights(y, num_atoms = grid_len);
  fmy = kwm$f1y;

  m1step = function (pt) {
    weights = pmax(pmin(1 - (1-pt)*f0y/fmy, 0.99), 0.01)
    muinit = mean(mt)/pt
    binit = lm(mosaic::logit(pmax(pmin(mt/muinit, 0.9), 0.1)) ~ 0 + x)$coefficients;
    robj = lgstep(y, x, weights, binit, grid_len, blambda, level)
    robj$pi0 = pt
    return(robj)
  }

  if (verbose) {res = pblapply(pi0grid, m1step)}
  else {res = lapply(pi0grid, m1step)}

  ll_list = sapply(res, function (robj) robj$ll)
  bi = which.max(ll_list)
  robj = res[[bi]]
  robj$ll_list = ll_list

  return(robj)
}
