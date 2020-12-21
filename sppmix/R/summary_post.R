#' Summarize DAMCMC results
#'
#' @description
#' Prints a brief summary of the results of a DAMCMC fit.
#'
#' For examples see
#'
#' \url{http://faculty.missouri.edu/~micheasa/sppmix/sppmix_all_examples.html
#' #summary.damcmc_res}
#'
#' @param object Object of class \code{damcmc_res}.
#' @param burnin Number of initial
#' realizations to discard. By default, it
#' is 1/10 of the total number of iterations.
#' @param alpha Level alpha for the credible
#' sets. Default is 0.05, for 95% credible
#' sets of the mixture parameters.
#' @param dgt Number of digits to use
#' (formatting the output).
#' @param ... Additional arguments for the S3 method.
#'
#' @author Jiaxun Chen, Sakis Micheas, Yuchen Wang
#' @seealso \code{\link{rnormmix}},
#' \code{\link{to_int_surf}},
#' \code{\link{rsppmix}},
#' \code{\link{est_mix_damcmc}},
#' \code{\link[spatstat.geom]{owin}}
#' @examples
#' \donttest{
#' # generate data
#' truemix<- rnormmix(m = 3, sig0 = .1, df = 5, xlim= c(0, 5), ylim = c(0, 5))
#' summary(truemix)
#' intsurf=to_int_surf(truemix, lambda = 100, win =spatstat.geom::owin( c(0, 5),c(0, 5)))
#' pp1 = rsppmix(intsurf = intsurf)# draw points
#' #Run DAMCMC and get posterior realizations
#' postfit=est_mix_damcmc(pp1,m=3)
#' #summary of the posterior results
#' summary(postfit)}
#'
#' @export
#' @method summary damcmc_res
summary.damcmc_res <- function(object, burnin = object$L / 10, alpha=0.05,dgt = 4,...) {
  fit <- drop_realization(object, burnin)
  if(fit$L==0)
  {
    cat("\nNo realizations left for",fit$m,
        "component(s) after burn-in (",
        burnin,") is applied.\n")
    stop()
  }
  m = fit$m
  for (i in 1:m) {
    #true value and credible sets
    poststats = GetStats_sppmix(fit$genps[,i], alpha=alpha)
    poststats <- sapply(poststats, format, digits = dgt)

    cat("---------------- Component ", i , "----------------\n")
    cat("Probability: posterior mean =",
              poststats$Mean,"\n")
    cat(poststats$CredibleSetConfidence,
              "% Credible Set:\n[",poststats$CredibleSet[1],
              ",",poststats$CredibleSet[2],"]\n" )

    poststats = GetStats_sppmix(fit$genmus[i,1,],alpha=alpha)
    poststats <- sapply(poststats, format, digits = dgt)

    cat("\nMean vector, x-coord: post mean =",
              poststats$Mean,"\n")
    cat(poststats$CredibleSetConfidence,
              "% Credible Set:\n[",poststats$CredibleSet[1],
              ",",poststats$CredibleSet[2],"]\n" )

    poststats = GetStats_sppmix(fit$genmus[i,2,],alpha=alpha)
    poststats <- sapply(poststats, format, digits = dgt)

    cat("\nMean vector, y-coord: post mean =",
              poststats$Mean,"\n")
    cat(poststats$CredibleSetConfidence,
              "% Credible Set:\n[",poststats$CredibleSet[1],
              ",",poststats$CredibleSet[2],"]\n" )

    sigs=fit$gensigmas[,i];
    L=fit$L;
    sigs11=vector("double",L)
    sigs12=vector("double",L)
    sigs22=vector("double",L)
    for(j in 1:L)
    {
      sigs11[j]=sigs[[j]][1,1]
      sigs12[j]=sigs[[j]][1,2]
      sigs22[j]=sigs[[j]][2,2]
    }
    poststats=GetStats_sppmix(sigs11,alpha=alpha)
    poststats <- sapply(poststats, format, digits = dgt)
    cat("\nCovariance, (1,1):\npost mean =",poststats$Mean,"\n")
    cat(poststats$CredibleSetConfidence, "% Credible Set:\n[",poststats$CredibleSet[1],
          ",",poststats$CredibleSet[2],"]\n" )

    poststats=GetStats_sppmix(sigs12,alpha=alpha)
    poststats <- sapply(poststats, format, digits = dgt)
    cat("\nCovariance, (1,2) and (2,1):\npost mean =",poststats$Mean,"\n")
    cat(poststats$CredibleSetConfidence, "% Credible Set:\n[",poststats$CredibleSet[1],
              ",",poststats$CredibleSet[2],"]\n" )

    poststats=GetStats_sppmix(sigs22,alpha=alpha)
    poststats <- sapply(poststats, format, digits = dgt)
    cat("\nCovariance, (2,2):\npost mean =",poststats$Mean,"\n")
    cat(poststats$CredibleSetConfidence, "% Credible Set:\n[",poststats$CredibleSet[1],
              ",",poststats$CredibleSet[2],"]\n" )

  }
}

#' Summarize BDMCMC results
#'
#' @description
#' Prints a brief summary of the results of a BDMCMC fit.
#'
#' For examples see
#'
#' \url{http://faculty.missouri.edu/~micheasa/sppmix/sppmix_all_examples.html
#' #summary.bdmcmc_res}
#'
#' @param object Object of class \code{bdmcmc_res}.
#' @param num_comp Number of components requested. Only
#' the posterior realizations that have this many components will be returned. The function
#' fails if the BDMCMC chain never visited this number of components.
#' We can also pass a vector of integer values and present the posterior means summary. If
#' this argument is missing, the MAP estimator is chosen by default.
#' @param burnin Number of initial
#' realizations to discard. By default, it
#' is 1/10 of the total number of iterations.
#' @param alpha Level alpha for the credible
#' sets. Default is 0.05, for 95% credible
#' sets of the mixture parameters.
#' @param dgt Number of digits to use
#' (formatting the output).
#' @param ... Additional arguments for the S3 method.
#' @examples
#' \donttest{
#' # generate data
#' truemix<- rnormmix(m = 3, sig0 = .1, df = 5, xlim= c(0, 5), ylim = c(0, 5))
#' summary(truemix)
#' intsurf=to_int_surf(truemix, lambda = 100, win =spatstat.geom::owin( c(0, 5),c(0, 5)))
#' pp1 <- rsppmix(intsurf = intsurf)# draw points
#' #Run BDMCMC and get posterior realizations
#' postfit=est_mix_bdmcmc(pp1,m=5)
#' #summary of the posterior results
#' summary(postfit)
#' summary(postfit, num_comp=2)
#' summary(postfit, num_comp=c(2,4))}
#'
#' @author Sakis Micheas
#' @seealso \code{\link{rnormmix}},
#' \code{\link{to_int_surf}},
#' \code{\link{rsppmix}},
#' \code{\link{est_mix_bdmcmc}},
#' \code{\link[spatstat.geom]{owin}}
#'
#' @export
#' @method summary bdmcmc_res
summary.bdmcmc_res <- function(object, num_comp,burnin = floor(object$L / 10), alpha=0.05,dgt = 4,...) {
  tab<-GetBDTable(object,FALSE)
  fit <- drop_realization(object, burnin)
  if(missing(num_comp))
  {
    num_comp=tab$MAPcomp
    cat("\nMissing num_comp. Will show the MAP realizations.\n")
  }
  if(length(num_comp)==1 && is.numeric(num_comp))
  {
    if(num_comp<1 || num_comp>fit$maxnumcomp)
    {
      cat("\nBad number of components requested\n")
      stop()
    }

    if(tab$FreqTab[num_comp]==0)
    {
      GetBDTable(fit)
      cat("\nNo realizations for",num_comp,"component(s).\n")
      stop()
    }
    keep=(fit$numcomp==num_comp)
    if(!any(keep==TRUE))
    {
      cat("\nNo realizations left for",num_comp,
          "component(s) after burn-in (",
          burnin,") is applied.\n")
      stop()
    }
    m = num_comp
    cat("\nNumber of components is",m,"\n")
    BDfitcomp=GetBDCompfit(fit,m,burnin = 0)
    summary.damcmc_res(BDfitcomp$BDgens,
        burnin=burnin,alpha=alpha,dgt = dgt)
  }
  else
  if(is.vector(num_comp))
  {
      for(k in 1:length(num_comp))
      {
        if(num_comp[k]<1 || num_comp[k]>fit$maxnumcomp)
        {
          cat("\nBad number of components requested\n")
          next
        }
        if(tab$FreqTab[num_comp[k]]==0)
        {
          cat("\nNo realizations for",num_comp[k],"mixture component(s).\n")
          next
        }
        keep=(fit$numcomp==num_comp[k])
        if(!any(keep==TRUE))
        {
          cat("\nNo realizations left for",num_comp[k],
              "component(s) after burn-in (",
              burnin,") is applied.\n")
          next
        }
        m=num_comp[k]
        cat("\n-----------------------------")
        cat("\nNumber of components is",m,"\n")
        BDfitcomp=GetBDCompfit(fit,m,burnin = 0)
        summary.damcmc_res(BDfitcomp$BDgens,
            burnin=burnin,alpha=alpha,dgt = dgt)
      }
  }
  else
    stop("Pass an integer or vector of integers to this function.")
}
