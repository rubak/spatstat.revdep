MakeMixtureList<- function(allgens_List,burnin=1000)
{
  #takes the DAMCMC output and makes a mixture list
  #based on its means, sppmix::MakeMixtureList(gens)
  post_means=GetAllMeans_sppmix(allgens_List,burnin);
  m=length(post_means$meanps);
  mix=vector("list", m);
  for(i in 1:m)
  {
    mix[[i]]=list(p = post_means$meanps[i],
                  mu = post_means$meanmus[i,],
                  sigma = post_means$meansigmas[,,i]);
  }
  return (mix)
}

#' Retrieve parts of a BDMCMC fit
#'
#' @description
#' The function can be used to obtain the realizations and
#' the corresponding surface of posterior means, for a specific
#' number of components. Use \code{\link{GetPMEst}} if you want just the surface.
#'
#' For examples see
#'
#' \url{http://faculty.missouri.edu/~micheasa/sppmix/sppmix_all_examples.html
#' #GetBDCompfit}
#'
#' @param BDfit Object of class \code{damcmc_res}.
#' @param num_comp Number of components requested. Only
#' the posterior realizations that have this many components will be returned. The function
#' fails if the BDMCMC chain never visited this number of components.
#' @param burnin Number of initial realizations to discard. By default, it is 1/10 of the total number of iterations.
#'
#' @return A list containing the following:
#'  \item{BDgens}{Realizations corresponging to this many mixture components. This is a \code{damcmc_res} object (same as the result of a \code{\link{est_mix_damcmc}} call). All realizations for the requested number of components are returned, that is, burnin is not applied to this object.}
#'  \item{BDsurf}{For the requested \code{num_comp}, this is the Poisson intensity surface based on the corresponding posterior means (label switching might be present).}
#'  \item{BDnormmix}{For the requested \code{num_comp}, this is a \code{\link{normmix}} object containing the corresponding ps, mus and sigmas (label switching might be present).}
#' @author Sakis Micheas
#' @seealso \code{\link{est_mix_bdmcmc}},
#' \code{\link{GetBDTable}},
#' \code{\link{plot.damcmc_res}},
#' \code{\link{plot.normmix}}
#' @examples
#' \donttest{
#' fit <- est_mix_bdmcmc(pp = spatstat::redwood, m = 7)
#' GetBDTable(fit)
#' #retrieve all BDMCMC realizations corresponding to a mixture with 5 components
#' BDfit5comp=GetBDCompfit(fit,5)
#' plot(BDfit5comp$BDsurf,main="Mixture intensity surface with 5 components")
#' #plot with the correct window
#' plot(BDfit5comp$BDnormmix,xlim =BDfit5comp$BDsurf$window$xrange,ylim =
#'  BDfit5comp$BDsurf$window$yrange )
#' plot(BDfit5comp$BDgens)}
#'
#' @export
GetBDCompfit<- function(BDfit,num_comp,burnin=floor(BDfit$L/10))
{
  if(burnin>=BDfit$L)
  {
    burnin = floor(BDfit$L/10)
  }
  tab<-GetBDTable(BDfit,FALSE)
  if(missing(num_comp))
  {
    num_comp=tab$MAPcomp
    cat("\nMissing num_comp. Will collect the MAP realizations.\n")
  }
  if(length(num_comp)==1 && is.numeric(num_comp))
  {
    if(num_comp<1 || num_comp>max(tab$FreqTab))
    {
      cat("\nBad number of components requested\n")
      stop()
    }
    if(tab$FreqTab[num_comp]==0)
    {
      GetBDTable(BDfit)
      cat("\nNo realizations for",num_comp,"component(s).\n")
      stop()
    }
  }
  else
  {
    GetBDTable(BDfit)
    cat("\num_comp must be single integer treated as the number of components.\n")
    stop()
  }

  keep=(BDfit$numcomp==num_comp)
  L=sum(keep)
  if(L==0)
  {
    cat("\nNo realizations left for",num_comp,
        "component(s) after burn-in (",
        burnin,") is applied.\n")
    stop()
  }
  if(L<burnin)
  {
    cat("\nNumber of realizations smaller than the burn-in applied. Burnin set to zero.\n")
#    stop()
    burnin=0
  }
  gens=BDfit
  #  L=length(gens$genlamdas[keep])
  newgens<-list(allgens_List=gens$allgens_List[keep],
                genps=as.matrix(gens$genps[keep,],nrow=L,
                                ncol=num_comp),
                genmus=gens$genmus[,,keep],
                gensigmas=gens$gensigmas[keep,],
                genzs=as.matrix(gens$genzs[keep,],nrow=L,
                                ncol=BDfit$data$n),
                genlamdas=as.matrix(gens$genlamdas[keep],L,1),
                ApproxCompMass=gens$ApproxCompMass[keep,],
                Badgen=as.matrix(gens$Badgen[keep],L,1),
                data=BDfit$data,
                L=L,
                m=num_comp
  )
  #trim all realizations to length num_comp
  newgens$genps=newgens$genps[,1:num_comp, drop = FALSE]
  newgens$genmus=newgens$genmus[1:num_comp,,, drop = FALSE]
  newgens$gensigmas=newgens$gensigmas[,1:num_comp, drop = FALSE]
  newgens$ApproxCompMass=newgens$ApproxCompMass[,1:num_comp, drop = FALSE]
  for(i in 1:newgens$L)
  {
    allgens<-vector("list", num_comp);
    currentgen=newgens$allgens_List[[i]]
    for(j in 1:num_comp)
    {
      allgens[[j]]=currentgen[[j]]
    }
    newgens$allgens_List[[i]]=allgens
  }
  class(newgens) <- "damcmc_res"

  mix_of_postmeans=MakeMixtureList_sppmix(
    newgens$allgens_List,burnin)
  BDnormmix=MakeNormMixFromMixtureList(mix_of_postmeans)
  BDsurf=to_int_surf(BDnormmix,
                     lambda=mean(newgens$genlamdas[(burnin+1):newgens$L]),
                     win=BDfit$data$window)
  return(list(BDgens=newgens,BDsurf=BDsurf,BDnormmix=BDnormmix))
}

#' Retrieve the MAP and distribution of the number of components
#'
#' @description
#' The function can be used to obtain the MAP estimate (mode of the posterior)
#' along with the frequency table for the
#' number of components, based on a BDMCMC fit from \code{\link{est_mix_bdmcmc}}.
#'
#' For examples see
#'
#' \url{http://faculty.missouri.edu/~micheasa/sppmix/sppmix_all_examples.html
#' #GetBDTable}
#'
#' @param BDfit A BDMCMC fit obtain from \code{\link{est_mix_bdmcmc}}.
#' @param showtable Logical variable requesting to display the frequency table. Default is TRUE.
#' @return A list containing the following:
#'  \item{MAPcomp}{The MAP number of mixture components.}
#'  \item{FreqTab}{Frequency table for the number of components.}
#'  \item{MeanComp}{The posterior mean for the number of components.}
#' @author Sakis Micheas
#' @seealso \code{\link{est_mix_bdmcmc}}
#' @examples
#' \donttest{
#' fit <- est_mix_bdmcmc(pp = spatstat::redwood, m = 7)
#' GetBDTable(fit)}
#'
#' @export
GetBDTable<- function(BDfit,showtable=TRUE)
{
  if(showtable)
  {
    cat("Frequency table for the number of components")
    print(table(BDfit$numcomp))
  }
  tab=tabulate(BDfit$numcomp,nbins=BDfit$maxnumcomp)
  MAPcompList=GetMax_sppmix(tab)
  MAPcomp=MAPcompList$pos+1;
  if(max(BDfit$numcomp)==1)
    MeanComp=sum(1:max(BDfit$numcomp)*tab/BDfit$L)
  else
    MeanComp=1
  return(list(MAPcomp=MAPcomp,FreqTab=tab,MeanComp=MeanComp))
}

#' Retrieve the MAP estimates for the component labels
#'
#' @description
#' The function returns the Maximum A Posteriori (MAP)
#' estimates of the component labels (membership indicator variables)
#' based on a \code{damcmc_res} object (output from \code{\link{est_mix_damcmc}}) or
#' a \code{bdmcmc_res} object (output from \code{\link{est_mix_bdmcmc}}) for the
#' chain corresponding to MAP number of components.
#'
#' For examples see
#'
#' \url{http://faculty.missouri.edu/~micheasa/sppmix/sppmix_all_examples.html
#' #GetMAPLabels}
#'
#' @param fit Object of class \code{damcmc_res} or \code{bdmcmc_res}.
#' @return A vector with size equal to the number of points, containing
#' the MAP estimators of the component labels (or
#' membership indicator variables). This the most likely
#' component we would classify a point in.
#' @seealso \code{\link{normmix}},
#' \code{\link{to_int_surf}},
#' \code{\link{rsppmix}},
#' \code{\link{est_mix_damcmc}}
#' @examples
#' \donttest{
#' truemix <- normmix(ps=c(.4, .2,.4), mus=list(c(0.3, 0.3), c(.5,.5),c(0.7, 0.7)),
#'  sigmas = list(.02*diag(2), .05*diag(2), .01*diag(2)))
#' intsurf=to_int_surf(truemix,lambda = 100, win = spatstat.geom::square(1))
#' pp1 <- rsppmix(intsurf)
#' plot(pp1)
#' plot(pp1, mus = intsurf$mus)#plot the mixture means as well
#' #plot the points with different colors depending on the true component label
#' plot(pp1, colors = TRUE)
#' #plot the points with different colors depending on the estimated component label
#' fit <- est_mix_damcmc(pp1, 3)
#' est_comp <- GetMAPLabels(fit)
#' plot(pp1, estcomp = est_comp, colors = TRUE)
#' fitBD <- est_mix_bdmcmc(pp1, 5)
#' est_compBD <- GetMAPLabels(fitBD)
#' plot(pp1, estcomp = est_compBD, colors = TRUE)}
#'
#' @author Jiaxun Chen
#' @export
GetMAPLabels<-function(fit)
{
  if(class(fit)=="bdmcmc_res")
  {
    tab=GetBDTable(fit,F)
    fit=GetBDCompfit(fit,tab$MAPcomp)$BDgens
  }
  est_comp <- apply(fit$genzs, 2, function(x) which.max(table(x)))
  return(est_comp)
}

#' Retrieves basic Bayesian estimates from a generated chain
#'
#' @description
#' The function returns the posterior mean and Credible Set
#' for a parameter based on a chain of posterior realizations.
#'
#' For examples see
#'
#' \url{http://faculty.missouri.edu/~micheasa/sppmix/sppmix_all_examples.html
#' #GetStats}
#'
#' @param chain A Markov Chain (a vector) containing the posterior realizations of the parameter.
#' @param alpha Level to use for the credible set.
#' @return A list containing the min, max, mean, Credible Set and CredibleSetConfidence level.
#' @seealso \code{\link{normmix}},
#' \code{\link{to_int_surf}},
#' \code{\link{rsppmix}},
#' \code{\link{est_mix_damcmc}}
#' @examples
#' \donttest{
#' truemix <- normmix(ps=c(.4, .2,.4), mus=list(c(0.3, 0.3), c(.5,.5),c(0.7, 0.7)),
#'  sigmas = list(.02*diag(2), .05*diag(2), .01*diag(2)))
#' intsurf=to_int_surf(truemix,lambda = 100, win = spatstat.geom::square(1))
#' pp1 <- rsppmix(intsurf)
#' fit <- est_mix_damcmc(pp1, 3)
#' p1=GetStats(fit$genps[,1])
#' p1$Mean
#' p1$CredibleSet
#' p2=GetStats(fit$genps[,2])
#' p2$Mean
#' p2$CredibleSet
#' p3=GetStats(fit$genps[,3])
#' p3$Mean
#' p3$CredibleSet}
#'
#' @author Sakis Micheas
#' @export
GetStats<-function(chain,alpha=0.05)
{
  GetStats_sppmix(chain,alpha)
}
