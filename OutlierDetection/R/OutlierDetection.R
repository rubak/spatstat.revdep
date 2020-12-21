#' Outlier Detection(Intersection of all the methods)
#'
#' Takes a dataset and finds its outliers using combination of different method
#' @param x dataset for which outliers are to be found
#' @param k No. of nearest neighbours to be used for for outlier detection using bootstrapping, default value is 0.05*nrow(x)
#' @param cutoff Percentile threshold used for distance, default value is 0.95
#' @param Method Distance method, default is Euclidean
#' @param rnames Logical value indicating whether the dataset has rownames, default value is False
#' @param depth  Logical value indicating whether depth based method should be used or not, default is False
#' @param dense  Logical value indicating whether density based method should be used or not, default is False
#' @param distance  Logical value indicating whether distance based methods should be used or not, default is False
#' @param dispersion  Logical value indicating whether dispersion based methods should be used or not, default is False

#' @details OutlierDetection finds outlier observations for the data using different methods and based on all the methods considered, labels an observation as outlier(intersection of all the methods). For bivariate data, it also shows the scatterplot of the data with labelled outliers.
#' @return Outlier Observations: A matrix of outlier observations
#' @return Location of Outlier: Vector of Sr. no. of outliers
#' @author Vinay Tiwari, Akanksha Kashikar

#' @examples
#' OutlierDetection(iris[,-5])

OutlierDetection=function(x,k=0.05*nrow(x),cutoff=.95,Method="euclidean",rnames=FALSE,depth=FALSE,dense=FALSE,distance=FALSE,dispersion=FALSE)
{

  data=as.data.frame(x)
  out=OutlierDetection::maha(x,cutoff=.95)$'Location of Outlier'
  if(dispersion==FALSE)
  {out=out
  }else
  {
    out2=OutlierDetection::disp(x,cutoff=.95)$'Location of Outlier'

    out=intersect(out2,out)
  }


  if(distance==FALSE)
  {out=out
  }else
  {
    out3=OutlierDetection::nn(x,k=0.05*nrow(x),cutoff=.95,Method="euclidean",rnames=F)$'Location of Outlier'
    out4=OutlierDetection::nnk(x,k=0.05*nrow(x),cutoff=.95,Method="euclidean",rnames=F)$'Location of Outlier'
    out=intersect(intersect(out3,out),intersect(out4,out))
  }

  if(depth==FALSE)
  {out=out
  }else
  {
  out5=OutlierDetection::depthout(x)$'Location of Outlier'
  out=intersect(out5,out)
  }
  if(dense==FALSE)
  {out=out
  }else
  {
    out5=OutlierDetection::dens(x)$'Location of Outlier'
    out=intersect(out5,out)
  }
  d=1:nrow(data)
  Class=c()
  Class=rep("Usual",length(d))
  for (i in 1:length(out)) {

    Class[d==out[i]]="Outlier"
  }
  cols <- c("Outlier" = "red", "Usual" = "blue")

  if(ncol(x)==2)
  {

    if(rnames==T)
    {
    s=subset(data,Class=="Outlier")
    gplot=ggplot2::ggplot(data,aes(data[,1],data[,2]))+geom_point(aes(colour=Class,pch=Class))+geom_text(data=s,aes(x=s[,1],y=s[,2],label=rownames(s)),colour="Red", hjust = "inward",check_overlap = T)+ggtitle("Outlier plot")+xlab("Variable1")+ylab("Variable2")+scale_color_manual(values=cols)

    }else{
      dd=cbind(data,1:nrow(data))
    s=subset(dd,Class=="Outlier")
    gplot=ggplot2::ggplot(data,aes(data[,1],data[,2]))+geom_point(aes(colour=Class,pch=Class))+geom_text(data=s,aes(x=s[,1],y=s[,2],label=s[,3]),colour="Red", hjust = "inward",check_overlap = T)+ggtitle("Outlier plot")+xlab("Variable1")+ylab("Variable2")+scale_color_manual(values=cols)
    }
    Out=x[out,]
    l=list("Outlier Observations"=Out,"Location of Outlier"=out,"Scatter plot"=gplot)
  }else if(ncol(x)==3)
  {


    plot=plotly::plot_ly(x=data[,1],y=data[,2],z=data[,3],type="scatter3d",mode="markers",color=Class,colors=c("Red","Blue"))
    Out=x[out,]

    l=list("Outlier Observations"=Out,"Location of Outlier"=out,"3Dplot"=plot)
  }else if(ncol(x)==4)
  {

    plot=plotly::plot_ly(x=data[,1],y=data[,2],z=data[,3],size = data[,4],type="scatter3d",mode="markers",color=Class,colors=c("Red","Blue"))
    Out=x[out,]

    l=list("Outlier Observations"=Out,"Location of Outlier"=out,"3Dplot"=plot)
  }else{
    Out=x[out,]
  l=list("Outlier Observations"=Out,"Location of Outlier"=out)}
return(l)
}
