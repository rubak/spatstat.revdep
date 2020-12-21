#' Plots the file type used as input and output of simulacro().
#'
#' This function makes a plot of the data frame used for argument INIDIST of function simulacro() and of the output of simulacro(). The plot shows the points as circles and crosses to represent their natural and anthropic origin respectively and has an iformative legend.
#'
#' @importFrom graphics axis legend par points title
#'
#' @param x simulacro-type file (as outputted by simulacro() or inputted as argument INIDIST).
#' @param tr threshold for Pnat above which points are considered of human nature.
#' @param outline geographical boundaries to x. Object of class sp::SpatialPolygons or data frame or matrix with the outline's corners.
#' @param main main title of the plot, by default is 'Map of ...' where ...= x$species (the name of the modelled species).
#'
#' @return plot.
#'
#' @author Luca Butikofer
#'
#' @export
#'
#' @examples
#' data(frogsEM)
#' plotlacro(x= frogsEM, outline= nzp)


plotlacro<-function(x, tr=.5, outline, main='std'){

  par.original <- par(no.readonly=TRUE)  #original par values
  par(mar= c(5, 4, 4, 13) + 0.1,xpd=T,family='Courier')

  if(class(outline)=='data.frame' | class(outline)=='matrix'){
    plot(outline, type='l')
  }else{
    plot(outline)
  }

  if (main=='std'){  #standard title
    title(paste('Map of',x[1,4]),font.axis=3,
          font.main=3,cex.main=1.4, lwd=.85)
  }else{
    title(main,font.axis=3,
          font.main=3,cex.main=1.4, lwd=.85)
  }


  points(x[x$Pnat>=tr,3],x[x$Pnat>=tr,2],pch=1)
  points(x[x$Pnat<tr,3],x[x$Pnat<tr,2],pch=3)
  axis(1);axis(2)



  legend(x= (extent(outline)[2] + (extent(outline)[2]-extent(outline)[1])/10),
         y= extent(outline)[4],
         c(paste('Natural, n=',nrow(x[x$Pnat>=tr,])),
           paste('Human, n=',nrow(x[x$Pnat<tr,]))),
         pch=c(1,3),
         bty='n',cex=1.1)

  par(par.original)  #restore par values
}
