# calculate q
qhat <- function( i, m, mean_log_ps, Adeltamean, Athetamean) {
  ksi <- log(sum(sapply(1:m, function(u) {exp(mean_log_ps[u] - Athetamean[u]
                                              + Adeltamean[u] - 1)})))
  q_i <- exp(-ksi + mean_log_ps[i] + Adeltamean[i] - Athetamean[i] - 1)
  return(q_i)
}

# partial derivative of function A
dAdth <- function(thetas, par, comp) {
  t1 <- thetas[1, comp]
  t2 <- thetas[2, comp]
  t3 <- thetas[3, comp]
  t4 <- thetas[4, comp]
  t5 <- thetas[5, comp]
  if(par == 1) {
    RVAL <- 0.5*(t5^2 - t2)/(t1*t2 - t3^2) - 0.5*t2*(t1*t2 - t3^2)^(-2)*
      (t2*t4^2 + t1*t5^2 - 2*t3*t4*t5)
  }
  if(par == 2) {
    RVAL <- 0.5*(t4^2 - t1)/(t1*t2 - t3^2) - 0.5*t1*(t1*t2 - t3^2)^(-2)*
      (t2*t4^2 + t1*t5^2 - 2*t3*t4*t5)
  }
  if(par == 3) {
    RVAL <- (t3 - t4*t5)/(t1*t2 - t3^2) + t3*(t1*t2 - t3^2)^(-2)*
      (t2*t4^2 + t1*t5^2 - 2*t3*t4*t5)
  }
  if(par == 4) {
    RVAL <- (t1*t2 - t3^2)^(-1)*(t2*t4 - t3*t5)
  }
  if(par == 5) {
    RVAL <- (t1*t2 - t3^2)^(-1)*(t1*t5 - t3*t4)
  }
  return(RVAL)
}

# Build objective function PEL
# function A
A <- function(thetas, i, is_theta) {
  if(is_theta == 1) {
    t1 <- thetas[[1]][,i]
    t2 <- thetas[[2]][,i]
    t3 <- thetas[[3]][,i]
    t4 <- thetas[[4]][,i]
    t5 <- thetas[[5]][,i]
  } else {
    t1 <- thetas[1, i]
    t2 <- thetas[2, i]
    t3 <- thetas[3, i]
    t4 <- thetas[4, i]
    t5 <- thetas[5, i]
  }
  RVAL <- -0.5*log(t1*t2 - t3^2) +
    0.5*(t1*t2 - t3^2)^(-1)*(t2*t4^2 + t1*t5^2 - 2*t3*t4*t5)
  return(RVAL)
}

get_theta <- function(sigmas, mus, intvl, m) {
  theta1 <- matrix(0, (intvl), m)
  theta2 <- matrix(0, (intvl), m)
  theta3 <- matrix(0, (intvl), m)
  theta4 <- matrix(0, (intvl), m)
  theta5 <- matrix(0, (intvl), m)
  for(i in 1:(intvl)) {
    temp <- sapply(sigmas[i, ],solve)
    theta1[i, ] <- temp[1, ]
    theta2[i, ] <- temp[4, ]
    theta3[i, ] <- temp[2, ]
    theta4[i, ] <- temp[1, ]*mus[, 1, i]+temp[2, ]*mus[, 2, i]
    theta5[i, ] <- temp[2, ]*mus[, 1, i]+temp[4, ]*mus[, 2, i]
  }
  RVAL <- list(theta1,theta2,theta3,theta4,theta5)
}

# define the function of calculating KL estimate of ps
getqs <- function(post1, burnin, intvl, m,i) {
  mus = post1$genmus[ , ,(burnin+intvl*(i-1)+1):(burnin+intvl*(i)), drop= FALSE]
  sigmas = post1$gensigmas[(burnin+intvl*(i-1)+1):(burnin+intvl*(i)), ,
                           drop = FALSE]
  ps = post1$genps[(burnin+intvl*(i-1)+1):(burnin+intvl*(i)), ,
                   drop = FALSE]

  thetas <- get_theta(sigmas, mus, intvl = intvl, m=m)

  # calculate deltas
  deltas <- t(matrix(sapply(thetas, colMeans),nrow = m, ncol = 5,
                     dimnames=list(paste("comp",1:m), paste("theta",1:5))))


  # calculate the mean of A(theta)
  Athetamean=rep(0,m)
  Atheta=matrix(0,(intvl),m)
  for(i in 1:m) {
    Athetamean[i] = mean(A(thetas, i, 1))
    Atheta[,i] = A(thetas, i, 1)
  }
  # calculate A(delta)
  Adeltamean=rep(0,m)
  for(i in 1:m) {
    Adeltamean[i]=A(deltas, i, 0)
  }


  dAdel <- matrix(0, 5, m)
  for(i in 1:m) {
    for(j in 1:5) {
      dAdel[j, i] <- dAdth(deltas, j, i)
    }
  }

  mean_log_ps <- colMeans(log(ps))

  qs <- sapply(1:m,function(x) qhat(x,m,mean_log_ps,Adeltamean,Athetamean))
  return(qs)
}

#' Retrieve the surface of Kullback-Leibler (KL) estimators
#'
#' @description
#' This function calculates the Kullback-Leibler
#' estimators of the parameters of the components of the
#' mixture intensity, based on a DAMCMC or BDMCMC fit. This is a decision
#' theoretic estimator of the parameters, meaning that, we compute the
#' Posterior Expected Loss (PEL) using the KL loss function and then
#' find the parameter values that minimize the PEL.
#'
#' For examples see
#'
#' \url{http://faculty.missouri.edu/~micheasa/sppmix/sppmix_all_examples.html
#' #GetKLEst}
#'
#' @param fit Object of class \code{damcmc_res} or \code{bdmcmc_res}.
#' @param burnin Number of initial realizations to discard. By default, it is 1/10 of the total number of iterations.
#' @param fixLS Logical requesting to check and fix label switching (if present). Default is FALSE.
#' @param approx Logical flag to request use of the identifiability constraint
#' to permute all realizations. Same parameter as in
#' function \code{\link{FixLS_da}}.
#' @param segment Number of segments to split the posterior realizations into. Each portion of posterior
#' realizations is used to calculate a single Kullback-Leibler realization. The KL estimator is the average of
#' all the KL realizations. Default is 50.
#' @return An object of class \code{intensity_surface}.
#' @author Jiaxun Chen, Sakis Micheas
#' @seealso \code{\link{rmixsurf}},
#' \code{\link{rsppmix}},
#' \code{\link{est_mix_damcmc}},
#' \code{\link{GetPMEst}},
#' \code{\link{GetMAPEst}},
#' \code{\link{CompareSurfs}}
#' @examples
#' \donttest{
#' #generate a surface
#' truemix_surf <- rmixsurf(m = 3, lambda=100, xlim = c(-3,3), ylim = c(-3,3))
#' plot(truemix_surf,main="True IPPP intensity surface")
#' #generate a point pattern
#' genPPP=rsppmix(intsurf = truemix_surf, truncate = FALSE)
#' #fit the IPPP model using DAMCMC
#' fit = est_mix_damcmc(genPPP, m = 3,L=20000)
#' #get the surfaces of posterior means, MAP and KL estimates
#' Meansest=GetPMEst(fit)
#' MAPest=GetMAPEst(fit)
#' KLest=GetKLEst(fit)
#' #plot all fitted surfaces
#' plot(Meansest,main="IPPP intensity surface of posterior means")
#' plot(MAPest,main="IPPP intensity surface of MAP estimates")
#' plot(KLest,main="IPPP intensity surface of KL estimates")
#' #fix labels (if label switching is detected)
#' KLestLSFixed=GetKLEst(fit,fixLS=TRUE,approx=FALSE)
#' plot(KLestLSFixed,main="IPPP intensity surface of KL estimates (LS fixed)")
#' #compare the four estimates against the truth
#' CompareSurfs(truemix_surf, Meansest, LL = 100, truncate = FALSE)
#' CompareSurfs(truemix_surf, MAPest, LL = 100, truncate = FALSE)
#' CompareSurfs(truemix_surf, KLest, LL = 100, truncate = FALSE)
#' CompareSurfs(truemix_surf, KLestLSFixed, LL = 100, truncate = FALSE)}
#'
#' @export
GetKLEst <- function(fit,burnin=floor(fit$L/10),
                     fixLS=FALSE,approx=FALSE, segment = 50)
{
  intvl <- floor((fit$L - burnin)/segment)
  nintl <- (fit$L - burnin)/intvl
  burnined <- drop_realization(fit, burnin)
  if(class(fit)=="bdmcmc_res")
  {
    tab=GetBDTable(burnined,FALSE)
    burnined=GetBDCompfit(burnined,tab$MAPcomp)$BDgens
  }
  if(fixLS)
  {
    #permute the realizations if there is label switching
    if(check_labels(burnined,burnin=0))
      fit_burnined = FixLS_da(burnined,approx=approx,burnin=0,run_silent=TRUE)
    else
      fit_burnined=burnined
  }
  else
    fit_burnined=burnined
  gen_qs <- matrix(0, nintl, fit$m)
  for ( i in 1:nintl) {
    gen_qs[i, ] <- getqs(post1 = fit_burnined, burnin = 0, intvl = intvl, m = fit$m, i = i)
  }
  post_ps <- colMeans(gen_qs)
  mus <- apply(fit_burnined$genmus[1:fit$m, , ,drop = FALSE], 1:2, mean)

  sigmas <- apply(fit_burnined$gensigmas[, 1:fit$m, drop = FALSE], 2,
                  function(mats) Reduce(`+`, mats) / length(mats))

  mean_lambda <- mean(fit_burnined$genlamdas)

  post_mus <- post_sigmas <- vector("list", fit$m)
  for (i in seq_along(post_mus)) {
    post_mus[[i]] <- mus[i, ]
    post_sigmas[[i]] <- matrix(sigmas[, i], 2, 2)
  }
  normmix(post_ps, post_mus, post_sigmas, mean_lambda,
          fit$data$window, estimated = TRUE)

}
