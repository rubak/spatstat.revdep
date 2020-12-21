#' Jackknifes a dataset.
#'
#' This function performes jackknife resampling on a dataset.
#'
#' @importFrom utils write.table
#'
#' @param DF data frame, matrix or any R object wich responds to function rownames().
#' @param N number of desired jackknifed datasets.
#' @param PR proportion of entries to DF that will be kept in the jackknifed datasets. Default is 0.85.
#' @param DIR directory whre to save the jackknifed datasets. If FALSE (default) will not save to disk.
#'
#' @return list of jacknifed datasets. If DIR is specified also a folder in directory DIR containing one .RDS file per jackknifed dataset (with extension .jds - Jackknifed Data Set) will be created.
#'
#' @author Luca Butikofer
#'
#' @export
#'
#' @examples
#' data('frogs')
#' frogsJK<- jackKnife(DF= frogs, N= 10)


jackKnife<- function(DF, N, PR=.85, DIR=F){

  if(DIR!=FALSE) { dir.create(paste0(DIR,'/Jackknifed_datasets/')) }

  PR<- .85  #resampling proportion
  rownames(DF)<- NULL
  jackdata<- NULL  #bootstrapped data

  for (i in 1:N){
    jack<- DF[sample(rownames(DF),PR*nrow(DF),replace=F),]
    jack<- jack[order(jack$year),]
    jackdata[[i]]<- jack
    if(DIR!=FALSE) { write.table(jack,file = paste(DIR,'/Jackknifed_datasets/Jack',i,'.jds',sep=''),row.names = F) }
  }

  return(jackdata)

}
