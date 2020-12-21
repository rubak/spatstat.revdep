#' Fix Label Switching
#'
#' @description
#' Permutes the posterior realizations in order to
#' fix the labels by either applying
#' an identifiability constraint or by minimizing
#' the squared error loss to find the best
#' permutation.
#'
#' For examples see
#'
#' \url{http://faculty.missouri.edu/~micheasa/sppmix/sppmix_all_examples.html
#' #FixLS_da}
#'
#' @param fit Object of class \code{damcmc_res} or \code{bdmcmc_res}.
#' @param burnin Number of initial realizations to discard. By default, it is 1/10 of the total number of iterations.
#' @param xlab The label for the x-axis.
#' @param ylab The label for the y-axis.
#' @param approx Logical flag to request use of the identifiability constraint
#' to permute all realizations. If FALSE, minimizing
#' the loss function can be very slow for moderate to large
#' number of components (m>10), since the algorithm
#' goes through all \code{m!} permutations for each
#' posterior realization.
#' @param plot_result Logical flag for requesting plots of the point pattern
#' and intensity surface based on the permuted
#' realizations. The default is FALSE.
#' @param run_silent Logical flag to hide progress messages. Default is FALSE.
#' @references Jasra, A., Holmes, C.C. and
#' Stephens, D. A. (2005). Markov Chain Monte
#' Carlo Methods and the Label Switching
#' Problem in Bayesian Mixtures. Statistical
#' Science, 20, 50-67.
#' @author Jiaxun Chen, Sakis Micheas
#' @seealso \code{\link{normmix}},
#' \code{\link{rsppmix}},
#' \code{\link{est_mix_damcmc}},
#' \code{\link{plot_chains}},
#' \code{\link{check_labels}}
#' @examples
#' \donttest{
#' # generate data
#' mix1 <- normmix(ps=c(.4, .2,.4), mus=list(c(0.3, 0.3), c(.5,.5),c(0.7, 0.7)),
#'  sigmas = list(.02*diag(2),.05*diag(2), .02*diag(2)),lambda = 100, win = spatstat::square(1))
#' #plot the true mixture
#' plot(mix1,main = "True Poisson intensity surface (mixture of normal components)")
#' pp1 <- rsppmix(mix1)
#' # Run Data augmentation MCMC and get posterior realizations
#' postfit = est_mix_damcmc(pp1, m=3, truncate=TRUE)
#' #plot the chains
#' plot_chains(postfit)
#' plot_chains(postfit,separate = FALSE)
#' # get the intensity of posterior means
#' post_mean = GetPMEst(postfit)
#' # plot the estimated intensity surface
#' plot(post_mean)
#' #check labels
#' check_labels(postfit)
#' # Fix label switching, start with approx=TRUE
#' post_fixed = FixLS_da(postfit, plot_result = TRUE)
#' plot_chains(post_fixed)
#' plot_chains(post_fixed, separate = FALSE)
#' #this one works better in theory
#' post_fixed = FixLS_da(postfit, approx=FALSE, plot_result = TRUE)
#' plot_chains(post_fixed)
#' plot_chains(post_fixed, separate = FALSE)}
#'
#' @export
FixLS_da<- function(fit, burnin = floor(fit$L / 10),
                 xlab = "x",ylab = "y",
                 approx = TRUE,
                 plot_result = FALSE,
                 run_silent = FALSE)
{
  if(burnin>=fit$L)
  {
    if(!run_silent)
    cat("\nBad burnin value. Using the default.")
    burnin = floor(fit$L/10)
  }
  fit <- drop_realization(fit, burnin)
  win <- spatstat::domain(fit$data)
  if(class(fit)=="bdmcmc_res")
  {
    tab=GetBDTable(fit,F)
    fit=GetBDCompfit(fit,tab$MAPcomp)$BDgens
    burnin = fit$L / 10
  }
  m <- fit$m
  if(m>5 && approx ==FALSE && !run_silent)
    cat("\nWARNING: m>5, this will take a long long time to complete...\n")
  xlims1 <- c(win$xrange)
  ylims1 <- c(win$yrange)
  L <- dim(fit$genps)[1]
  if (approx == TRUE)
  {
    if(!run_silent)
      cat("\nStarting the approximate relabeling algorithm.\n")
    iter=1
    lagnum=lagstart=min(.05*(fit$L-burnin),50)+iter
    while(check_labels(fit,lagnum=lagnum,showmessages=FALSE))
    {
      if(!run_silent)
        cat("\nIteration",iter,", lag=",lagnum,"\n")
      iter=iter+1
      lagnum=lagstart+iter*10
      permgens = PostGenGetBestPermIdenConstraint_sppmix(fit)
      fit$allgens_List = permgens$allgens_List
      fit$genps = permgens$genps
      fit$genmus = permgens$genmus
      fit$gensigmas = permgens$gensigmas
      fit$genzs = permgens$genzs
      fit$ApproxCompMass=PermuteZs_sppmix(fit$ApproxCompMass ,permgens$best_perm)
  #      fit$genlamdas = permgens$genlamdas
    }
  } else
  {
    if(!run_silent)
      cat("\nStarting the decision theoretic relabeling algorithm.\n")
    iter=1
    lagnum=lagstart=min(.05*(fit$L-burnin),50)+iter
    while(check_labels(fit,lagnum=lagnum,showmessages=FALSE))
    {
      if(!run_silent)
        cat("\nIteration",iter,", lag=",lagnum,"\n")
      iter=iter+1
      lagnum=lagstart+iter*10
      permgens = PostGenGetBestPerm_sppmix(fit$allgens_List)
      fit$allgens_List = permgens$permuted_gens
      fit$genps = permgens$permuted_ps
      fit$genmus = permgens$permuted_mus
      fit$gensigmas = permgens$permuted_sigmas
      fit$genzs=PermuteZs_sppmix(fit$genzs ,permgens$best_perm)
      fit$ApproxCompMass=PermuteZs_sppmix(fit$ApproxCompMass ,permgens$best_perm)
      #      fit$genlamdas = fit$genlamdas
    }
  }
  post_ps <- colMeans(fit$genps)
  mus <- apply(fit$genmus, 1:2, mean)
  mean_mat <- function(mats) Reduce("+", mats) / length(mats)
  sigmas <- apply(fit$gensigmas, 2, mean_mat)
  mean_lambda <- mean(fit$genlamdas)
  post_mus <- post_sigmas <- vector("list", m)
  for (i in 1:m) {
    post_mus[[i]] <- mus[i, ]
    post_sigmas[[i]] <- matrix(sigmas[, i], 2, 2)
  }
  post_normix = normmix(post_ps, post_mus, post_sigmas, lambda = mean_lambda,
                        win = win)
  if (plot_result == TRUE)
  {
    print(plot.sppmix(fit$data, post_normix$mus))
    print(plotmix_2d(post_normix,fit$data))
    if (approx == TRUE)
      plot.intensity_surface(post_normix,
                           main = "Surface of posterior means (Identifiability Constraint)")
    else
      plot.intensity_surface(post_normix,
                             main = "Surface of posterior means (Permuted Labels)")
  }
  return(fit)
}

#' Check for label switching
#'
#' @description
#' Checks if there is label switching present in the
#' posterior realizations using the chains for
#' mu. The algorithm is heuristic and works as follows;
#' for each chain of mu for a given
#' component, we look for sharp changes
#' in the generated values that lead to dramatically
#' different variability from the
#' variability observed in the past
#' history of the chain. The lag history is 5\% of the
#' total number of realizations, excluding the burnin realizations.
#'
#' For examples see
#'
#' \url{http://faculty.missouri.edu/~micheasa/sppmix/sppmix_all_examples.html
#' #check_labels}
#'
#' @param fit Object of class \code{damcmc_res} or \code{bdmcmc_res}.
#' @param burnin Number of initial realizations to discard. By default, it is 1/10 of the total number of iterations.
#' @param lagnum Number of past realizations to use for the detection algorithm. Default is .05 of the total realizations after burnin realizations are removed.
#' @param showmessages Logical variable requesting to show informative messages (TRUE) or suppress them (FALSE). Default is TRUE.
#' @references Jasra, A., Holmes, C.C. and Stephens, D. A. (2005). Markov Chain Monte Carlo Methods and the Label Switching Problem in Bayesian Mixtures. Statistical Science, 20, 50-67.
#' @author Jiaxun Chen, Sakis Micheas
#' @return Logical value indicating if label switching was detected.
#' @seealso \code{\link{normmix}},
#' \code{\link{rsppmix}},
#' \code{\link{est_mix_damcmc}},
#' \code{\link{plot_chains}},
#' \code{\link{FixLS_da}}
#' @details To avoid the label switching problem one can plot
#' the average of the intensities of the posterior realizations, i.e.,
#' the average of the surfaces instead of the surface of the posterior averages. However,
#' by doing so we lose the ability to perform mixture deconvolution,
#' namely, work with the posterior means of the ps, mus and sigmas of the
#' mixture intensity. In general, the average of posterior surfaces
#' for each realization of the ps, mus and sigmas, and the surface based on
#' the posterior means of the ps, mus and sigmas, need not be the same.
#'
#' For a DAMCMC fit, avoiding label switching can be accomplished using
#' function \code{\link{plot_avgsurf}}, which plots the average posterior surface.
#' The surface of posterior means is plotted using
#' the function \code{\link{plot.damcmc_res}} on the returned value of \code{\link{GetPMEst}}.
#'
#' When working with a BDMCMC fit, it is recommended that you use \code{\link{GetBDCompfit}} to retrieve
#' the realizations corresponding to a specific number of components and then fix the labels.
#' Of course the best estimate in this case, is the Bayesian model average intensity
#' obtained via \code{\link{GetBMA}}, but it can be very slow to compute. The BMA
#' is an average on surfaces based on each of the posterior realizations, and as such
#' does not suffer from the label switching problem.
#' @examples
#' \donttest{
#' # generate data
#' mix1 = normmix(ps=c(.4, .2,.4), mus=list(c(0.3, 0.3), c(.5,.5),c(0.7, 0.7)),
#'  sigmas = list(.02*diag(2),.05*diag(2), .02*diag(2)),lambda = 100,
#'  win = spatstat::square(1))
#' plot(mix1,main="True Poisson intensity surface (mixture of normal components)")
#' pp1 = rsppmix(mix1)
#' # Run Data augmentation MCMC and get posterior realizations
#' postfit = est_mix_damcmc(pp1,m=5)
#' #plot the chains
#' plot_chains(postfit)
#' #check labels
#' check_labels(postfit)
#' #plot the average posterior surface
#' plot(GetPMEst(postfit))
#' #plot the surface of posterior means, can be slow for large LL
#' avgsurf=plot_avgsurf(postfit, LL = 50, burnin = 1000)
#' # Fix label switching, start with approx=TRUE
#' post_fixed = FixLS_da(postfit, plot_result = TRUE)
#' plot_chains(post_fixed)
#' plot_chains(post_fixed,separate = FALSE)
#' #this one works better in theory
#' post_fixed = FixLS_da(postfit,approx=FALSE, plot_result = TRUE)
#' plot_chains(post_fixed)
#' plot_chains(post_fixed,separate = FALSE)}
#'
#' @export
check_labels<- function(fit, burnin = floor(fit$L/10),lagnum=floor(.05*(fit$L-burnin)+1),showmessages=TRUE ) {
  if(class(fit)=="bdmcmc_res")
  {
    tab=GetBDTable(fit,F)
    fit=GetBDCompfit(fit,tab$MAPcomp)$BDgens
    burnin = floor(fit$L/10)
    lagnum=floor(.05*(fit$L-burnin)+1)
  }
  m <- fit$m
  if(burnin>=fit$L)
  {
    cat("\nBad burnin or lagnum value. Using the default.")
    burnin = floor(fit$L/10)
    lagnum=floor(.05*(fit$L-burnin)+1)
  }
#  cat("\nburnin=",burnin)
  genmus <- drop_realization(fit, burnin)$genmus
  if(showmessages)
    cat("\nChecking for label switching...\n")
  for (i in 1:m) {
    if(Check4LabelSwitching_sppmix(genmus[i, 1, ],lagnum)) {
      if(showmessages)
        cat("Label switching present.\nPermute the labels to get a better fit,\nor obtain the average of the surfaces\n ")
      return(TRUE)
    }
    if(Check4LabelSwitching_sppmix(genmus[i, 2, ],lagnum)) {
      if(showmessages)
        cat("Label switching present.\nPermute the labels to get a better fit,\nor obtain the average of the surfaces\n ")
      return(TRUE)
    }
  }
  if(showmessages)
    cat("No Label switching detected. If the plots suggest label switching\nyou could rerun the detection algorithm with a smaller lag value.")
  return(invisible(FALSE))
}

