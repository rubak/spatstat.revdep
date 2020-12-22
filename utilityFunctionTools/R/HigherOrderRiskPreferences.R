# Authors
# Sebastian Schneider, sschneider@coll.mpg.de; sebastian@sebastianschneider.eu
# Giulia Baldini, giulia.baldini@uni-bonn.de

# Copyright (C) 2020 Sebastian O. Schneider & Giulia Baldini

# This program is free software; you can redistribute it and/or
# modify it under the terms of the GNU General Public License
# as published by the Free Software Foundation; either version 3
# of the License, or (at your option) any later version.

# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.

# You should have received a copy of the GNU General Public License
# along with this program; if not, write to the Free Software
# Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.


##########################################################
#### Utility Tools for Higher Order Risk Preferences #####
##########################################################

#' Truncated p-th power function. Helper function for creating the B-Spline basis (Code by Paul Eilers, Package JOPS, http://statweb.lsu.edu/faculty/marx/JOPS_0.1.0.tar.gz)
#' @param x Function value.
#' @param t Point of truncation.
#' @param p degree of the truncated polynomial function.
#' @return Returns a piece-wise defined basis functions for x > t.
#' @examples 
#' tpower(1, 2, 3)
#' @export
tpower <- function(x, t, p) {
  return ((x - t) ^ p * (x > t))
}


#' Constructs a B-spline basis of degree 'deg' (Code by Paul Eilers, Package JOPS, http://statweb.lsu.edu/faculty/marx/JOPS_0.1.0.tar.gz).
#'
#' @param x values for the x axis.
#' @param xl minimum value, default is the minimum value of the x-values.
#' @param xr maximum value, default is maximum value of the x-values.
#' @param ndx number of intervals to partition the distance between xl and xr.
#' @param deg degree of the B-spline basis.
#' @return a B-spline basis of degree deg and ndx + 1 internal knots.
#' @examples
#' x_finegrid <- seq(0.001, 1.0, (1.0 - 0.001) / 1000)
#' bbase(x_finegrid)
#' @export
bbase <- function(x,
                  xl = min(x),
                  xr = max(x),
                  ndx = 20,
                  deg = 6) {
  dx <- (xr - xl) / ndx
  knots <- seq(xl - deg * dx, xr + deg * dx, by = dx)
  P <- outer(x, knots, tpower, deg)
  n <- dim(P)[2]
  D <- diff(diag(n), diff = deg + 1) / (gamma(deg + 1) * dx ^ deg)
  B <- (-1) ^ (deg + 1) * P %*% t(D)
  return(B)
}

#' Estimates the model
#'
#' @param xi a vector containing the certainty equivalents (x-values of utility points) for a given participant in each use case.
#' @param yi can be a vector or a matrix representing the corresponding utility values (y-values of utility points).
#' @param lambda lambda is the penalization weight used to compute the initial estimate. The default value is 1.
#' @param n_penalty_dimensions number of dimensions (i.e., derivatives) to penalize. Possible values are 1 or 2. The default value is 1.
#' @param penalty_order highest dimension (i.e., derivative) to penalize. Must be lower than deg.
#' @param ndx number of intervals to partition the distance between the lowest and highest x-values of the utility points.
#' @param deg degree of the B-spline basis. Determines the degree of the function to be estimated. If deg = 2, the estimated utility function will consist of quadratic functions. 
#' @param cross_validation_mode determines which cross validation mode should be used. If 0, then the cross validation method is leave-one-third-out. If 1, then the cross validation method is a theoretical leave-one-out, i.e., based on a formula. The default value is 1.
#' @param return_estimate parameter that indicates whether or not to return the (initially) estimated coefficients. Default is false.
#' @param left_out_xi needed for cross validation: the x-values of the points that are left out for fitting the model, so that they can be predicted
#' @param left_out_yi needed for cross validation: the y-values of the points that are left out for fitting the model, so that they can be predicted
#' @return Returns the sum of residuals of the prediction of the left-out points using cross validation. If specified, additionally returns the estimated coefficients of the utility function (in the B-spline basis).
#' @examples
#' x <- c(0.0000000, 0.2819824, 0.3007812, 0.4375000, 0.5231934, 0.7784882, 0.8945312, 1.0000000)
#' y <- c(0.0000, 0.1250, 0.2500, 0.5000, 0.6250, 0.6875, 0.7500, 1.0000)
#' estimate_model(x, y, .5)
#' @importFrom stats lsfit
#' @export
estimate_model <- function(xi,
                           yi,
                           lambda = 1,
                           n_penalty_dimensions = 1,
                           penalty_order = 4,
                           ndx = 20,
                           deg = 6,
                           cross_validation_mode = 0,
                           return_estimate = 0,
                           left_out_xi = c(),
                           left_out_yi = c()) {
  if (n_penalty_dimensions > penalty_order) {
    n_penalty_dimensions = penalty_order
  }
  
  lambda <- lambda * 1000 ^ -(seq(0, n_penalty_dimensions - 1))
  
  # Enforce U(0) = 0 and U(1) = 1
  w1 <- 0
  w2 <- 0
  repeat {
    xi <- c(rep(0, w1) , xi, rep(1, w2))
    yi <- c(rep(0, w1) , yi, rep(1, w2))
    
    # Generate base and penalty
    B = bbase(xi, ndx = ndx, deg = deg)
    n = ncol(B)
    
    # Fit the model - poor algorithm, using one penalty
    P <- NULL
    sumP <- matrix(rep(0, n ^ 2), nrow = n)
    nix <- NULL
    # Fit the model using lsfit
    for (p in 1:n_penalty_dimensions) {
      D = diff(diag(n), diff = penalty_order - (p - 1))
      P = rbind(P, sqrt(lambda[p]) * D)
      sumP = sumP + (lambda[p] * t(D) %*% D)
      nix = c(nix, rep(0, n - (penalty_order - (p - 1))))
    }
    
    f = lsfit(rbind(B, P), c(yi, nix), intercept = F)
    
    a <- f$coef
    
    if (n_penalty_dimensions == 2) {
      # Re-adjust lambda
      lambda_ratio <-
        sum(abs(diff(a, diff = penalty_order))) / sum(abs(diff(a, diff = penalty_order - 1)))
      if (length(lambda) > 1) {
        lambda <- lambda[1]
      }
      lambda <- c(lambda, lambda * lambda_ratio / 5)
      
      P <- NULL
      sumP <- matrix(rep(0, n ^ 2), nrow = n)
      nix <- NULL
      for (p in 1:n_penalty_dimensions) {
        D = diff(diag(n), diff = penalty_order - (p - 1))
        P = rbind(P, sqrt(lambda[p]) * D)
        sumP = sumP + (lambda[p] * t(D) %*% D)
        nix = c(nix, rep(0, n - (penalty_order - (p - 1))))
      }
    }
    
    # Enforce monotonicity
    kappa <- 100000000
    for (i in 1:10) {
      D1 <- diff(diag(n), diff = 1)
      W <- as.vector(1 * (D1 %*% a < 0))
      P2 <- sqrt(kappa) * diag(W) %*% D1
      nix2 <- rep(0, n - 1)
      f_new <- lsfit(rbind(B, P, P2), c(yi, nix, nix2), intercept = F)
      a_new <- f_new$coef
      if (identical(a, a_new)) {
        break
      }
      a = a_new
    }
    
    # Predict at observed points
    z = B %*% a
    
    if (cross_validation_mode) {
      # == 1, so leave one out
      # Predict at missing points
      x_finegrid <- seq(min(xi), max(xi), (max(xi) - min(xi)) / 1000)
      y_hat <- bbase(x_finegrid, ndx = ndx, deg = deg) %*% a
    }
    
    # Increase weight at (0,0) and (-1,-1) if not predicted with the desired precision
    if ((round(z[1], 2) == 0) & round(z[length(z)], 2) == 1) {
      break
    } else {
      if (round(z[1], 2) != 0) {
        w1 <- w1 + 1
      }
      if (round(z[length(z)], 2) != 1) {
        w2 <- w2 + 1
      }
    }
  }
  
  if (return_estimate) {
    return (a)
  }
  
  # Model choice using Cross Validation
  if (cross_validation_mode) {
    # == 1, so leave one out
    lhs <- t(B) %*% B + sumP + kappa * t(D1) %*% diag(W) %*% D1
    H <- lsfit(lhs, diag(nrow = nrow(lhs)), intercept = F)
    h_new <- diag(B %*% H$coeff %*% t(B))
    r = (yi[abs(xi) > 0 &
              abs(xi) < 1] - z[abs(xi) > 0 &
                                 abs(xi) < 1]) / (1 - h_new[abs(xi) > 0 &
                                                              abs(xi) < 1])
    res <- sqrt(sum(r ^ 2))
  } else {
    # == 0, so 1/3
    # Predict at left-out points
    z <- bbase(c(left_out_xi, xi), ndx = ndx, deg = deg) %*% a
    res <- sum((left_out_yi - z[1:length(left_out_yi)]) ^ 2)
  }
  return(res)
}

#' Evaluates the cross validation function.
#'
#' @param xi a vector containing the certainty equivalents (x-values of utility points) for a given participant in each use case.
#' @param yi can be a vector or a matrix representing the corresponding utility values (y-values of utility points).
#' @param lambda lambda is the penalization weight used to compute the initial estimate. The default value is 1.
#' @param n_penalty_dimensions number of dimensions (i.e., derivatives) to penalize. Possible values are 1 or 2. The default value is 1.
#' @param penalty_order highest dimension (i.e., derivative) to penalize. Must be lower than deg.
#' @param ndx number of intervals to partition the distance between the lowest and highest x-values of the utility points.
#' @param deg degree of the B-spline basis. Determines the degree of the function to be estimated. If deg = 2, the estimated utility function will consist of quadratic functions. 
#' @param cross_validation_mode determines which cross validation mode should be used. If 0, then the cross validation method is leave-one-third-out. If 1, then the cross validation method is a theoretical leave-one-out, i.e., based on a formula. The default value is 1.
#' @return Returns, for the given utility points and (possibly default) settings, the predictive quality of the estimated utility function according to cross validation as a function of a specified penalty weight lambda.
#' @examples 
#' x <- c(0.0000000, 0.2819824, 0.3007812, 0.4375000, 0.5231934, 0.7784882, 0.8945312, 1.0000000)
#' y <- c(0.0000, 0.1250, 0.2500, 0.5000, 0.6250, 0.6875, 0.7500, 1.0000)
#' evaluate_cross_validation(x, y, .5)
#' @export
evaluate_cross_validation <- function(xi,
                                      yi,
                                      lambda = 1,
                                      n_penalty_dimensions = 1,
                                      penalty_order = 4,
                                      ndx = 20,
                                      deg = 6,
                                      cross_validation_mode = 0) {
  if (!cross_validation_mode) {
    # == 0, so 1/3 out
    var_xi <- xi[abs(xi) > 0 & abs(xi) < 1]
    var_yi <- yi[abs(xi) > 0 & abs(xi) < 1]
    number_xi <- length(var_xi)
    # Leave 1/3 out CV
    combinations <-
      combn(1:number_xi, number_xi - ceiling(1 * number_xi / 3))
    number_combn <- dim(combinations)[[2]]
    sum_ssr <- 0
    for (c in 1:number_combn) {
      sum_ssr <- sum_ssr + estimate_model(
        xi = c(xi[which(abs(xi) == 0 |
                          abs(xi) == 1)], var_xi[combinations[, c]]),
        yi = c(yi[which(abs(xi) == 0 |
                          abs(xi) == 1)], var_yi[combinations[, c]]),
        lambda = lambda,
        n_penalty_dimensions = n_penalty_dimensions,
        penalty_order = penalty_order,
        ndx = ndx,
        deg = deg,
        cross_validation_mode = 0,
        left_out_xi = var_xi[-combinations[, c]],
        left_out_yi = var_yi[-combinations[, c]]
      )
    }
    avg_ssr = sum_ssr / number_combn
    return (avg_ssr)
  } else {
    return(
      estimate_model(
        xi,
        yi,
        lambda = lambda,
        n_penalty_dimensions = n_penalty_dimensions,
        penalty_order = penalty_order,
        ndx = ndx,
        deg = deg,
        cross_validation_mode = 1
      )
    )
  }
}
#' Finds an optimal penalty weight lambda given the parameters
#' @param xi a vector containing the certainty equivalents (x-values of utility points) for a given participant in each use case.
#' @param yi can be a vector or a matrix representing the corresponding utility values (y-values of utility points).
#' @param lambda_max maximum lambda used for computing the optimal lambda. The default value is 10000.
#' @param n_penalty_dimensions number of dimensions (i.e., derivatives) to penalize. Possible values are 1 or 2. The default value is 1.
#' @param penalty_order highest dimension (i.e., derivative) to penalize. Must be lower than deg.
#' @param ndx number of intervals to partition the distance between the lowest and highest x-values of the utility points.
#' @param deg degree of the B-spline basis. Determines the degree of the function to be estimated. If deg = 2, the estimated utility function will consist of quadratic functions. 
#' @param cross_validation_mode determines which cross validation mode should be used. If 0, then the cross validation method is leave-one-third-out. If 1, then the cross validation method is a theoretical leave-one-out, i.e., based on a formula. The default value is 1.
#' @param grid_dim dimension of the search grid for the initial grid search before the actual optimization. Default value is 5.
#' @return the optimal lambda for the given set of utility points and (possibly default) settings according to the specified cross validation method. 
#' @examples 
#' x <- c(0.0000000, 0.2819824, 0.3007812, 0.4375000, 0.5231934, 0.7784882, 0.8945312, 1.0000000)
#' y <- c(0.0000, 0.1250, 0.2500, 0.5000, 0.6250, 0.6875, 0.7500, 1.0000)
#' find_optimal_lambda(x, y)
#' @importFrom stats optim
#' @importFrom utils combn
#' @importFrom spatstat.geom pairdist ppp
#' @export 
find_optimal_lambda <- function(xi,
                                yi,
                                lambda_max = 10000,
                                n_penalty_dimensions = 1,
                                penalty_order = 4,
                                ndx = 20,
                                deg = 6,
                                cross_validation_mode = 0,
                                grid_dim = 5) {
  # Take distance between points
  dist <- (pairdist(ppp(xi, yi)) <= .1) * 1
  colsums <- apply(dist, 1, sum)
  skip <- 0
  num_balls <- 0
  for (p in 1:length(colsums)) {
    if (skip == 0) {
      skip <- colsums[p]
      num_balls <- num_balls + 1
    }
    skip <- max(0, skip - 1)
  }
  lambda_min = max(0.01, (num_balls * (9 / length(xi)) - 1) ^ 2.5)
  
  interval <- c(lambda_min, lambda_max + 1)
  argmin <- 1
  lower.lim <- lambda_min
  upper.lim <- lambda_max
  
  # Do a grid search first
  for (i in 1:grid_dim) {
    vals <-
      exp(seq(log(interval[1]), log(interval[2]), length.out = grid_dim))
    min_grid <- 2e+308
    min_index <- 1
    for (j in 1:length(vals)) {
      est <- evaluate_cross_validation(
        xi = xi,
        yi = yi,
        lambda = vals[j],
        n_penalty_dimensions = n_penalty_dimensions,
        penalty_order = penalty_order,
        ndx = ndx,
        deg = deg,
        cross_validation_mode = cross_validation_mode
      )
      if (est < min_grid) {
        min_grid <- est
        min_index <- j
      }
    }
    argmin <- vals[min_index]
    lower.lim <- pmax(argmin - lambda_max / 10 ^ i, lambda_min)
    upper.lim <- pmin(argmin + lambda_max / 10 ^ i, lambda_max)
    interval <- c(lower.lim[1], upper.lim[1])
  }
  
  
  optim_argmin <-
    optim(
      argmin,
      evaluate_cross_validation,
      xi = xi,
      yi = yi,
      n_penalty_dimensions = n_penalty_dimensions,
      penalty_order = penalty_order,
      ndx = ndx,
      deg = deg,
      cross_validation_mode = cross_validation_mode,
      method = "L-BFGS-B",
      lower = lower.lim,
      upper = upper.lim,
      control = list(
        fnscale = 1,
        maxit = 30,
        trace = T
      )
    )
  
  return (optim_argmin$par)
}

#' Computes a continuous and smooth utility function from the given utility points
#'
#' @param x a matrix or dataframe containing the certainty equivalents (x-values of utility points) for a given participant in each use case.
#' @param y can be a vector or a matrix representing the corresponding utility values (y-values of utility points).
#' @param ids a list containing the IDs of the participants. If not given, a list with IDs from 1 to n_observations will be created.
#' @param mode an integer between 0, 1, 2 representing the three possible modes: multiple imputation, optimal classification or 'weak' classification. Default is optimal classification (1).
#' @param penalty_order highest dimension (i.e., derivative) to penalize. Must be lower than deg.
#' @param lambda_max maximum lambda used for computing the optimal lambda. It is used only in multiple imputation (mode = 0) and optimal (mode = 1). The default value is 10000.
#' @param current_lambda lambda considered in the current iteration. Only used in multiple imputation (mode = 0) to create the combinations and as actual lambda value in 'weak' classification mode (mode = 2). The default value is 1.
#' @param ndx number of intervals to partition the distance between the lowest and highest x-values of the utility points.
#' @param deg degree of the B-spline basis. Determines the degree of the function to be estimated. If deg = 2, the estimated utility function will consist of quadratic functions. 
#' @param verbose shows some information while the program is running.
#' @return A smooth and continuous utility function.
#' @examples
#' \donttest{
#' x <- matrix(c(24.60938,34.76074,78.75,81.86035,128.5156, 
#'               7.109375,80.4248,113.75,115.083,135.0781, 
#'               3.828125,7.211914,8.75,124.1064,131.7969, 
#'               1.640625,2.084961,8.75,36.94824,98.98438), nrow = 4, ncol = 5, byrow = TRUE)
#' y <- c(0.25, 0.375, 0.5, 0.625, 0.75)
#' compute_function(x, y, verbose = 1)
#' }
#' @export
compute_function <- function(x,
                             y,
                             ids = NULL,
                             mode = 1,
                             penalty_order = 4,
                             lambda_max = 10000,
                             current_lambda = 1,
                             ndx = 20,
                             deg = 6,
                             verbose = 0) {
  if (is.data.frame(x)) {
    # If the data is a dataframe
    x <- as.matrix(sapply(x, as.numeric))   # convert to matrix
  } else if (!is.matrix(x)) {
    stop("Please convert x to a dataframe or to a matrix before calling this function.")
  }
  
  if (is.data.frame(y)) {
    # If y is a dataframe
    if (nrow(y) == 1) {
      y <- as.vector(sapply(x, as.numeric))   # convert to vector
    } else {
      y <- as.matrix(sapply(x, as.numeric))   # convert to matrix
    }
  }
  
  if (!is.matrix(y) & !is.vector(y)) {
    stop("The accepted values for y are: matrix or vector.")
  }
  
  if (is.vector(y) && ncol(x) != length(y) || is.matrix(y) && ncol(x) != ncol(y)) {
    stop("The y values do not have the same number of columns as the x values.")
  }
  
  if (length(penalty_order) > 1) {
    stop("The order of penalty should be a single integer.")
  }
  
  if (!is.null(ids) & !is.vector(ids) & !is.list(ids)) {
    stop("Please convert the ids field to a list or a vector.")
  }
  
  if (!verbose) {
    if(.Platform$OS.type == "unix") {
      sink(file = "/dev/null") # use /dev/null in UNIX
    } else {
      sink(file = "NUL") # use /dev/null in UNIX
    }
  }
  
  if (is.null(ids)) {
    ids = seq(from = 1,
              to = nrow(x),
              by = 1)
  }

  if (deg <= penalty_order) {
    stop(
      paste(
        "A degree value of ",
        deg,
        " with penalty order",
        penalty_order,
        " is not valid. The degree must always be greater than the order of penalty.",
        sep = ""
      )
    )
  }
  
  if (penalty_order > 4) {
    warning(
      "This program has not been tested with order of penalties higher than 4 and it might produce wrong or unexpected results."
    )
  }
  
  if (mode == 0 || mode == 1) {
    if (mode == 0) {
      message("Using multiple imputation mode.")
    } else {
      message("Using optimal classification mode.")
      if (penalty_order != 3 & penalty_order != 4) {
        stop("The orders of penalty for optimal classification can either be 3 or 4.")
      }
    }
    n_penalty_dimensions = 2
  } else {
    message("Using weak classification mode.")
    n_penalty_dimensions = 1
  }
  
  xi <- x[1, ]
  xi <- xi[!is.na(xi)]
  min_length_xi <- length(xi)
  max_length_xi <- length(xi)
  for (i in 2:nrow(x)) {
    # Find the xi of minimum length
    xi <- x[i, ]
    xi <- xi[!is.na(xi)]
    min_length_xi <- min(min_length_xi, length(xi))
    max_length_xi <- max(max_length_xi, length(xi))
  }
  
  if (ndx < min_length_xi | ndx >= 2.5 * max_length_xi) {
    warning(
      paste(
        "The value of ndx should be larger than or equal to the minimum length of xi existing in the dataset",
        " and smaller than 2.5 times the maximum length of xi existing in the dataset.",
        "Other values have not been tested to yield reliable results.",
        "Consider allowing a higher (lower) degree of flexibility and increase (decrease) the value of ndx."
      )
    )
  }
  
  # deg is larger than 1, but still lower than 2.5 * min(length(xi))
  if (deg < 1 | deg >= 1.2 * min_length_xi)  {
    warning(
      paste(
        "The value of deg should be larger than or equal to 1 and smaller than 1.2 times the minimum length of xi existing in the dataset.",
        "The amount of data is likely too small to fit such a flexible model, consider lowering deg."
      )
    )
  }
  
  return (
    compute_function_aux(
      x,
      y,
      ids,
      mode,
      penalty_order,
      lambda_max,
      current_lambda,
      n_penalty_dimensions,
      ndx,
      deg
    )
  )
}

#' Computes a continuous and smooth function according to the given utility points
#' @keywords internal
#' 
#' @param x a matrix or dataframe containing the certainty equivalents (x-values of utility points) for a given participant in each use case.
#' @param y can be a vector or a matrix representing the corresponding utility values (y-values of utility points).
#' @param ids a list containing the IDs of the participants. If not given, a list with IDs from 1 to n_observations will be created.
#' @param mode an integer between 0, 1, 2 representing the three possible modes: multiple imputation, optimal classification or 'weak' classification. Default is optimal classification (1).
#' @param penalty_order highest dimension (i.e., derivative) to penalize. Must be lower than deg.
#' @param lambda_max maximum lambda used for computing the optimal lambda. It is used only in multiple imputation (mode = 0) and optimal (mode = 1). The default value is 10000.
#' @param current_lambda lambda considered in the current iteration. Only used in multiple imputation (mode = 0) to create the combinations and as actual lambda value in 'weak' classification mode (mode = 2). The default value is 1.
#' @param n_penalty_dimensions number of dimensions to penalise. Possible values are 1 or 2. The default value is 1.
#' @param ndx number of intervals to partition the distance between the lowest and highest x-values of the utility points.
#' @param deg degree of the B-spline basis. Determines the degree of the function to be estimated. If deg = 2, the estimated utility function will consist of quadratic functions. 
#' @return A smooth and continuous utility function.
compute_function_aux <- function(x,
                                  y,
                                  ids,
                                  mode = 1,
                                  penalty_order = 4,
                                  lambda_max = 10000,
                                  current_lambda = 1,
                                  n_penalty_dimensions = 1,
                                  ndx = 20,
                                  deg = 6) {
  nval_grid <- 1000
  coeffs <- data.frame(row.names = ids, matrix(NA, nrow = nrow(x), ncol = (ndx + deg)))
  x_finegrids <- data.frame(row.names = ids, matrix(0, nrow = nrow(x), ncol = (nval_grid + 1)))
  for (i in 1:nrow(x)) {
    cross_validation_mode <- 0 # means leave 1/3 out
    # Collect the needed observation
    xi <- x[i, ]
    yi <- y
    if (is.matrix(y)) {
      yi <- y[i, ]
    }
    id <- ids[i]
    
    
    # Consider only the columns which are not NA
    yi <- yi[!is.na(xi)]
    xi <- xi[!is.na(xi)]
    message(paste("Considering participant with id", id))
  
    # Rescale
    xi <- xi / max(abs(xi))
    yi <- yi / max(abs(yi[length(yi)]), abs(yi[1]))
    
    # Estimate utility curve if there is at least one point between (0,0) and (-1,-1) or (1,1), respectively
    if (length(xi[abs(xi) > 0 &
                  abs(xi) < 1]) < 1 || sum(is.na(xi))) {
      message(paste("The participant with id", id, "cannot be analysed."))
      next
    }
    # If total sum of points to combine is less than or equal to penalty order,
    # only leave-one-out CV possible according to formula
    if (length(xi) <= penalty_order) {
      if (mode == 0) {
        message(
          paste(
            "The participant with id ",
            id,
            " cannot be analysed with multiple imputation mode, and penalty order ",
            penalty_order,
            ".",
            sep = ""
          )
        )
        next
      }
      cross_validation_mode <- 1 # leave one out
    }
    
    if (penalty_order == 4 &
        n_penalty_dimensions == 1 & length(xi) == 3) {
      # In that case, cannot estimate the model. Solution: Duplicate the end-points
      xi <- c(0, xi, 1)
      yi <- c(0, yi, 1)
    }
    
    # If we are in optimal classification or in MI, we want to find the optimal lambda
    if (mode == 0 || mode == 1) {
      new_lambda = find_optimal_lambda(
        xi,
        yi,
        lambda_max,
        n_penalty_dimensions,
        penalty_order,
        ndx,
        deg,
        cross_validation_mode
      )
    } else {
      new_lambda = current_lambda
    }
    
    if (mode == 0) {
      # MI
      # Leave 1/3 out CV
      var_xi <- xi[abs(xi) > 0 & abs(xi) < 1]
      var_yi <- yi[abs(xi) > 0 & abs(xi) < 1]
      number_xi <- length(var_xi)
      combinations <-
        combn(1:number_xi, number_xi - ceiling(1 * number_xi / 3))
      number_combn <- dim(combinations)[[2]]
      if (number_combn < current_lambda) {
        current_lambda <- current_lambda %% number_combn
        current_lambda <- current_lambda + 1
      }
      coeff <-
        estimate_model(
          xi = c(xi[which(abs(xi) == 0 | abs(xi) == 1)], var_xi[combinations[, current_lambda]]),
          yi = c(yi[which(abs(xi) == 0 | abs(xi) == 1)], var_yi[combinations[, current_lambda]]),
          lambda = new_lambda,
          n_penalty_dimensions = n_penalty_dimensions,
          penalty_order = penalty_order,
          ndx = ndx,
          deg = deg,
          return_estimate = 1,
        )
    } else {
      coeff <-
        estimate_model(
          xi = xi,
          yi = yi,
          lambda = new_lambda,
          n_penalty_dimensions = n_penalty_dimensions,
          penalty_order = penalty_order,
          ndx = ndx,
          deg = deg,
          return_estimate = 1,
        )
    }
    
    x_finegrid <- seq(min(xi), max(xi), (max(xi) - min(xi)) / nval_grid)
    coeffs[id, ] <- coeff
    x_finegrids[id, ] <- x_finegrid
  }
  
  return(list(x_finegrids, coeffs))
}

#' Given a set of smooth and continuous functions, computes predefined and user-defined measures.
#' @keywords internal
#' @param x_grids a dataframe of vectors of x-values for a smooth and continuous function.
#' @param coeffs a dataframe of coefficients for a smooth and continuous function for each participant.
#' @param ids a list containing the IDs of the participants. If not given, a list with IDs from 1 to n_observations will be created.
#' @param ndx number of intervals to partition the distance between the lowest and highest x-values of the utility points.
#' @param deg degree of the B-spline basis. Determines the degree of the function to be estimated. If deg = 2, the estimated utility function will consist of quadratic functions. 
#' @param measures a vector of measures to be computed.
#' @param ... additional parameters for user-defined measures.
#' @return A set of measurements.
compute_measures_aux <- function(x_grids,
                                  coeffs,
                                  ids,
                                  ndx = 20,
                                  deg = 6,
                                  measures = c("risk-arrow-pratt", "crainich-eeckhoudt", "denuit-eeckhoudt"),
                                  ...)  {
  # TODO: Check if there is a problem with column names and functions
  output_measures <- data.frame(row.names = ids, matrix(0, nrow = length(ids), ncol = length(measures)))
  for (i in 1:nrow(coeffs)) {
    coeff <- as.numeric(coeffs[i,])
    x_finegrid <- as.numeric(x_grids[i,])
      
    if (all(is.na(coeff))){
      message(paste("The participant with id", ids[i], "cannot be analysed because all the coefficients are NAs."))
      next
    }
    
    # Computation of first derivative
    dy_rd <- derivative(x_finegrid, coeff, 1, ndx, deg)
    for (j in 1:length(measures)) {
      measure <- measures[[j]]
      if (mode(measure) == "function") {
        mes <- measure(x_finegrid, coeff, ndx, deg, ...)
        colnames(output_measures)[j] <- paste("custom-", j, sep="")
      } else {
        colnames(output_measures)[j] <- measure
        if (measure == "risk-arrow-pratt") {
          # Computation of second derivative
          ddy_rd <- derivative(x_finegrid, coeff, 2, ndx, deg)
          # Compute Risk Aversion measures by Pratt / Arrow
          mes <- -mean(ddy_rd, na.rm = T) / mean(dy_rd, na.rm = T)
        } else if (measure == "crainich-eeckhoudt") {
          # Computation of third derivative
          dddy_rd <- derivative(x_finegrid, coeff, 3, ndx, deg)
          # Compute Prudence intensity
          mes <- mean(dddy_rd, na.rm = T) / mean(dy_rd, na.rm = T)
        } else if (measure == "denuit-eeckhoudt") {
          # Computation of fourth derivative
          ddddy_rd <- derivative(x_finegrid, coeff, 4, ndx, deg)
          # Compute Temperance intensity
          mes <- -mean(ddddy_rd, na.rm = T) / mean(dy_rd, na.rm = T)
        } else {
          stop("The desired measure does not exist. Please use another measure.")
        }
      }
      output_measures[i, j] <- mes
    }
  }
  return (output_measures)
}

#' Given a set of smooth and continuous functions, computes predefined and user-defined measures.
#'
#' @param x_grids a dataframe of vectors of x values for a smooth and continuous function.
#' @param coeffs a dataframe of coefficients for a smooth and continous function for each participant.
#' @param ids a list containing the IDs of the participants. If not given, a list with IDs from 1 to n_observations will be created.
#' @param ndx number of intervals to partition the distance between the lowest and highest x-values of the utility points.
#' @param deg degree of the B-spline basis. Determines the degree of the function to be estimated. If deg = 2, the estimated utility function will consist of quadratic functions. 
#' @param measures a vector of measures to be computed.
#' @param ... additional parameters for user-defined measures.
#' @return A set of measurements.
#' @examples 
#' x <- rbind(seq(0.000002, 1.0, (1.0 - 0.000002) / 1000),
#'            seq(0.001, 1.0, (1.0 - 0.001) / 1000),
#'            seq(0.0004, 1.0, (1.0 - 0.0004) / 1000))
#' y <- rbind(seq(0.000002, 1.0, (1.0 - 0.000002) / 15),
#'            seq(0.001, 1.0, (1.0 - 0.001) / 15),
#'            seq(0.0004, 1.0, (1.0 - 0.0004) / 15))
#' compute_measures(x, y, ndx = 10, deg = 6)

#' # x_finegrid, coeff, ndx, deg are always there to be used
#' # The function should have additional unknown arguments (...) if the given parameters are not used
#' risk_arrow_pratt <- function(x_finegrid, coeff, ndx, deg){ 
#'   dy_rd <- derivative(x_finegrid, coeff, 1, ndx, deg)
#'   ddy_rd <- derivative(x_finegrid, coeff, 2, ndx, deg)
#'   return (-mean(ddy_rd, na.rm = TRUE) / mean(dy_rd, na.rm = TRUE))
#' }
#' measures = c("crainich-eeckhoudt", "denuit-eeckhoudt", risk_arrow_pratt)
#' compute_measures(x, y, ndx = 10, deg = 6, measures=measures)
#' @export
compute_measures <- function(x_grids,
                             coeffs,
                             ids = NULL,
                             ndx = 20,
                             deg = 6,
                             measures = c("risk-arrow-pratt", "crainich-eeckhoudt", "denuit-eeckhoudt"),
                             ...) {
  if (!is.data.frame(x_grids) & !is.matrix(x_grids)) {
    stop("Please convert x_grids to a dataframe or to a matrix before calling this function.")
  }
  
  if (!is.data.frame(coeffs) & !is.matrix(coeffs)) {
    stop("Please convert coeffs to a dataframe or to a matrix before calling this function.")
  }
  
  if (nrow(x_grids) != nrow(coeffs)){
    stop("The number of participants in the coefficients matrix and in the x matrix do not correspond.")
  }
  
  if (!is.null(ids) & !is.vector(ids) & !is.list(ids)) {
    stop("Please convert the ids field to a list or a vector.")
  }
  
  if (!is.null(ids) & length(ids) != ncol(x_grids)) {
    stop("The number of participants in the ids vector and in the x matrix do not correspond.")
  }
  
  if (ncol(coeffs) != (ndx + deg)){
    stop("The number of coefficients is not the same as ndx + deg.")
  }
  
  if (is.null(ids)) {
    ids = seq(from = 1,
              to = nrow(x_grids),
              by = 1)
  }
  
  return(compute_measures_aux(x_grids, coeffs, ids, ndx, deg, measures, ...))
}

#' Computes the derivative of a function
#'
#' @param x the x values for which the derivative should be computed.
#' @param coeffs the coefficient.
#' @param degree the degree of the derivative.
#' @param ndx number of intervals to partition the distance between the lowest and highest x-values of the utility points.
#' @param deg degree of the B-spline basis. Determines the degree of the function to be estimated. If deg = 2, the estimated utility function will consist of quadratic functions. 
#' @return the derivative of the specified degree.
#' @examples 
#' coeffs <- seq(0.000002, 1.0, (1.0 - 0.000002) / 25)
#' x <- seq(0.01, 1.0, (1.0 - 0.01) / 5)
#' derivative(x, coeffs)
#' @export
derivative <- function(x,
                       coeffs,
                       degree = 1,
                       ndx = 20,
                       deg = 6) {
  dB = bbase(x = x, ndx = ndx, deg = deg - degree)
  for (i in 1:degree) {
    coeffs <- (diff(coeffs) / ndx)
  }
  dy = dB %*% coeffs
  range = abs(min(dy) - max(dy))
  dy_rd <- round(dy / (range / 1000)) * (range / 1000)
  return (dy_rd)
}

#' Computes a continuous and smooth function according to the given utility points
#'
#' @param x a matrix or dataframe containing the certainty equivalents (x-values of utility points) for a given participant in each use case.
#' @param y can be a vector or a matrix representing the corresponding utility values (y-values of utility points).
#' @param ids a list containing the IDs of the participants. If not given, a list with IDs from 1 to n_observations will be created.
#' @param mode an integer between 0, 1, 2 representing the three possible modes: multiple imputation, optimal classification or 'weak' classification. Default is optimal classification (1).
#' @param penalty_orders vector or constant that contains the derivates that will be smoothened. The values in this vector should not be larger than 4.
#' @param ndx number of intervals to partition the distance between the lowest and highest x-values of the utility points.
#' @param deg degree of the B-spline basis. Determines the degree of the function to be estimated. If deg = 2, the estimated utility function will consist of quadratic functions. 
#' @param measures the utility based (intensity) measures to be computed.
#' @param ... additional parameters for user-defined measures.
#' @param root_filename filename containing the location of where the output files are going to be saved.
#' @param verbose shows some information while the program is running.
#' @return A smooth and continuous function.
#' @examples
#' \donttest{
#' x <- matrix(c(24.60938,34.76074,78.75,81.86035,128.5156, 
#'               7.109375,80.4248,113.75,115.083,135.0781, 
#'               3.828125,7.211914,8.75,124.1064,131.7969, 
#'               1.640625,2.084961,8.75,36.94824,98.98438), nrow = 4, ncol = 5, byrow = TRUE)
#' y <- c(0.25, 0.375, 0.5, 0.625, 0.75)
#' compute_higher_order_risk_preferences(x, y, mode = 1)
#' 
#' # could be used with root_filename argument: 
#' # outfile <- paste(dirname(getwd()), "/out/output", sep="")
#' compute_higher_order_risk_preferences(x, y, mode = 2, verbose = 1)
#' }
#' @importFrom utils write.csv
#' @export
compute_higher_order_risk_preferences <- function(x,
                                                  y,
                                                  ids = NULL,
                                                  mode = 0,
                                                  penalty_orders = c(4),
                                                  ndx = 20,
                                                  deg = 6,
                                                  measures = c("risk-arrow-pratt", "crainich-eeckhoudt", "denuit-eeckhoudt"),
                                                  ...,
                                                  root_filename = NULL,
                                                  verbose = 0) {
  if (is.data.frame(x)) {
    # If the data is a dataframe
    x <- as.matrix(sapply(x, as.numeric))   # convert to matrix
  } else if (!is.matrix(x)) {
    stop("Please convert x to a dataframe or to a matrix before calling this function.")
  }
  
  if (is.data.frame(y)) {
    # If y is a dataframe
    if (nrow(y) == 1) {
      y <- as.vector(sapply(x, as.numeric))   # convert to vector
    } else {
      y <- as.matrix(sapply(x, as.numeric))   # convert to matrix
    }
  }
  
  if (!is.matrix(y) & !is.vector(y)) {
    stop("The accepted values for y are: matrix or vector.")
  }
  
  if (is.vector(y) && ncol(x) != length(y) || is.matrix(y) && ncol(x) != ncol(y)) {
    stop("The y values do not have the same number of columns as the x values.")
  }
  
  if (!is.vector(penalty_orders) ||
      !is.numeric(penalty_orders)) {
    stop("The accepted values for the order of penalty are either a constant or a vector.")
  }
  
  if (!is.null(ids) & !is.vector(ids) & !is.list(ids)) {
    stop("Please convert the ids field to a list or a vector.")
  }
  
  if (!verbose) {
    if(.Platform$OS.type == "unix") {
      sink(file = "/dev/null") # use /dev/null in UNIX
    } else {
      sink(file = "NUL") # use /dev/null in UNIX
    }
  }
  
  if (is.null(ids)) {
    ids = seq(from = 1,
              to = nrow(x),
              by = 1)
  }
  
  for (order in penalty_orders) {
    if (deg <= order) {
      stop(
        paste(
          "A degree value of ",
          deg,
          " with penalty order",
          order,
          " is not valid. The degree must always be greater than the order of penalty.",
          sep = ""
        )
      )
    }
    if (order > 4) {
      warning(
        "This program has not been tested with order of penalties higher than 4 and it might produce wrong or unexpected results."
      )
    }
  }
  
  
  xi <- x[1, ]
  xi <- xi[!is.na(xi)]
  min_length_xi <- length(xi)
  max_length_xi <- length(xi)
  for (i in 2:nrow(x)) {
    # Find the xi of minimum length
    xi <- x[i, ]
    xi <- xi[!is.na(xi)]
    min_length_xi <- min(min_length_xi, length(xi))
    max_length_xi <- max(max_length_xi, length(xi))
  }
  
  if (ndx < min_length_xi | ndx >= 4 * (max_length_xi+2)) {
    warning(
      paste(
        "The value of ndx should be larger than or equal to the minimum length of xi existing in the dataset",
        " and not too large compared to the the maximum length of xi existing in the dataset (say, less than 4 times this number).",
        "Other values have not been tested to yield reliable results.",
        "Consider allowing a higher (lower) degree of flexibility and increase (decrease) the value of ndx."
      )
    )
  }
  
  # deg should be larger than 1, but still lower than 1.2 * min(length(xi))
  if (deg < 1 | deg > 1.2 * min_length_xi)  {
    warning(
      paste(
        "The value of deg should be larger than or equal to 1 and not exceed 1.2 times the minimum length of xi existing in the dataset.",
        "The amount of data is likely too small to fit such a flexible model, consider lowering deg."
      )
    )
  }
  
  mode_txt <- "weak"
  if (mode == 0 || mode == 1) {
    # Set range for lambda (minimum is determined by the data)
    lambda_max = 10000
    if (mode == 0) {
      mode_txt <- "MI"
      message("Using multiple imputation mode.")
      lambda_fix_loop_lambdas <- 1:15
    } else {
      mode_txt <- "opt"
      message("Using optimal classification mode.")
      for (order in penalty_orders) {
        if (order != 3 & order != 4) {
          stop("The orders of penalty for optimal classification can either be 3 or 4.")
        }
      }
      lambda_fix_loop_lambdas <- 1
    }
    n_penalty_dimensions = 2
  } else {
    message("Using weak classification mode.")
    lambda_fix_loop_lambdas <- c(.1, 1, 10, 20, 50, 100, 500, 750, 1000, 2000, 5000, 10000)
    n_penalty_dimensions = 1
  }
  
  # Loop over the penalization order
  for (penalty_order in penalty_orders) {
    # Start smoothing & classification
    message(paste("Smoothing over the", penalty_order, "order of penalty."))
    
    # Iterate over the lambdas, only once for optimization
    for (lambda_fix_loop in lambda_fix_loop_lambdas) {
      # Start smoothing & classification
      message("Computing the function for all individuals.")
      fct <-
        compute_function_aux(
          x,
          y,
          ids,
          mode,
          penalty_order,
          lambda_max,
          lambda_fix_loop,
          n_penalty_dimensions,
          ndx,
          deg
        )
      
      x_grids <- fct[[1]]
      coeffs <- fct[[2]]
      message("Computing the measures for all individuals.")
      out_measures <- compute_measures_aux(x_grids, coeffs, row.names(coeffs), ndx, deg, measures)
      if (!is.null(root_filename)){
        addition_filename <-
          paste(
            "_",
            penalty_order,
            "_",
            lambda_fix_loop,
            "_",
            mode_txt,
            "_",
            format(Sys.time(), "%y.%m.%d_%H:%M:%OS"),
            ".csv",
            sep = ""
          )
        
        filesave <- paste(root_filename, "_", "x", addition_filename, sep="")

        dir.create(dirname(filesave), showWarnings = FALSE)
        write.csv(x_grids, file = filesave)
        write.csv(coeffs, file = paste(root_filename, "_", "coeffs", addition_filename, sep=""))
        write.csv(out_measures, file = paste(root_filename, "_", "measures", addition_filename, sep=""))
      }
    }
  }
  
  if (!verbose) {
    sink()
  }
}
