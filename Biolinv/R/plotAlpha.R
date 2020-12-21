#' Plots the output of modSel().
#'
#' This function makes a plot of the data frame outputted by function modSel(). Useful to choose wich of the simulated datasets is the most similar to the observed one.
#'
#' @importFrom grDevices colorRampPalette
#' @importFrom graphics par points
#'
#' @param SSIM data frame output of modSel.
#' @param REP number of replicates.
#' @param BP when FALSE (default) plots a line representing the average of all replicates. When TRUE plots boxplots of all replicates for each alpha value.
#'
#' @return plot.
#'
#' @author Luca Butikofer
#'
#' @export
#'
#' @examples
#' data(frogsSum)
#' plotAlpha(SSIM= frogsSum, REP= 10)
#' plotAlpha(SSIM= frogsSum, REP= 10, BP=TRUE)


plotAlpha<- function(SSIM, REP, BP=FALSE){

  #graphics
  par.original <- par(no.readonly=TRUE)  #original par values
  par(family='Courier', font.axis=3, font.main=3,
      cex.main=1.4, lwd=.85, bty='n')

  pal <- colorRampPalette(c('grey30','grey80'))
  colori<-pal(REP)

  #preparation
  UAV<- unique(SSIM$compAlpha)  #unique alpha values
  NAV<- length(UAV)  #number of alpha values

    lines<- matrix(nrow=NAV, ncol=REP) #one Alpha value per row, one replicate per column

    for (i in 1:NAV){

      s<- SSIM[SSIM$compAlpha== UAV[i],]
      #s$similarity<- (s$similarity-min(s$similarity))/max(s$similarity-min(s$similarity))
      lines[i,]<- s$dissimilarity

    }

    av<- apply(lines,1,mean)  #average
    low<- which.min(av)  #minimum point

    #boxplot



    #first replica
    plot(UAV,lines[,1],type='n',col=colori[1],
         main='Dissimilarity Plot',
         xlab='Alpha', ylab='Dissimilarity',
         ylim=c(min(lines),max(lines)))

    if (BP=='TRUE'){
      for(i in 1:nrow(lines)){
        boxplot(lines[i,], add=T, at=UAV[i],boxwex=1.5,axes=F)
      }
    }else{
      lines(UAV, av)
      points(UAV[low], av[low])
    }

  par(par.original)  #restore par values

}

