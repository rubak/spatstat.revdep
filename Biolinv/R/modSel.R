#' Wrapper for function spatSim() which allows use on multiple datasets.
#'
#' This function uses spatSim() function to compare many simulated datasets with one, typically observed, data set.
#'
#' @importFrom utils flush.console
#'
#' @param M0 data frame containing 'y' and 'x' coordinates (projected coordinate system) to compare with M2.
#' @param M2 list of dataframes with same structure of M0, typically generated with simulacro() function.
#' @param WIN window of observation for the point patterns (MOD0 and MOD2)(see ?spatstat.geom::owin). Object of class 'sp::owin'.
#' @param RAD numeric vector of search distances for the K-function.
#' @param AV numeric vector of the Alpha values of the simulated datasets in the same order as in the list of argument M2. Used to save the output data frame.
#'
#' @return data frame with two columns. 'dissimilarity': dissimilarity values as computed by spatSim() function. 'compAlpha': same as AV.
#'
#' @author Luca Butikofer
#'
#' @export
#'
#' @examples
#' data(nzw)
#' data(frogsEM)  #see EM().
#' data(frogsLacro)  #see simulacro().
#'
#' \dontrun{
#' frogsSum<- modSel(WIN= nzw, M0= frogsEM, M2= frogsLacro,
#'  AV= c(2,3,4.5,7.5,11,15,20,25), RAD= seq(0,30000,1000))
#'  }


modSel<-function(WIN,M0,M2,AV,RAD){

  AV<- unique(AV)

  if(class(M2[[1]])=='list') {
    M3<- unlist(M2, recursive = F)

    alpha<- NULL
    for(i in 1:length(AV)){
      alpha[[i]]<- rep(AV[i],length(M2[[1]]))
    }
    AV<- unlist(alpha)

    }
  sv<-rep(NA,length(M3))

  for(i in 1:length(M3)){  #comparison loop
    flush.console()
    print(paste(i," of ",length(M3)))
    sv[i]<-spatSim(M0,M3[[i]],WIN,R=RAD)
  }

  #sv<-sv/min(sv,na.rm=T)  #re-scale sv with its minimum becoming 1
  results<-data.frame(dissimilarity=sv,compAlpha=AV)
  return(results)
}
