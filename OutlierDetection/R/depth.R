#' Outlier detection using depth based method
#'
#' Takes a dataset and finds its outliers using depth-based method
#' @param x dataset for which outliers are to be found
#' @param cutoff Percentile threshold used for depth, default value is 0.05
#' @param rnames Logical value indicating whether the dataset has rownames, default value is False
#' @param boottimes Number of bootsrap samples to find the cutoff, default is 100 samples

#' @details depthout computes depth of an observation using depthTools package and based on the bootstrapped cutoff, label an observation as outlier. Outlierliness of the labelled 'Outlier' is also reported and it is the bootstrap estimate of probability of the observation being an outlier. For bivariate data, it also shows the scatterplot of the data with labelled outliers.
#' @return Outlier Observations: A matrix of outlier observations
#' @return Location of Outlier: Vector of Sr. no. of outliers
#' @return Outlier probability: Vector of proportion of times an outlier exceeds local bootstrap cutoff
#' @references Johnson, T., Kwok, I., and Ng, R.T. 1998. Fast computation of 2-dimensional depth contours. In Proc. Int. Conf. on Knowledge Discovery and Data Mining (KDD), New York, NY. Kno
#' @author Vinay Tiwari, Akanksha Kashikar

#' @examples
#' #Create dataset
#' X=iris[,1:4]
#' #Outlier detection
#' depthout(X,cutoff=0.05)


depthout=function(x,rnames=FALSE,cutoff=0.05,boottimes=100)
{

  data=as.data.frame(x)
  a=depthTools::MBD(data)
  d=a$MBD
  quanorig=quantile(d,cutoff)
  #ub=quantile(d,cutoff)
  k=density(d)
  aa=which(k$x<=quanorig)
  a=max(aa)
  b=which.max(k$x>=quanorig)
  f=((k$y[b]-k$y[a])/(k$x[b]-k$x[a]))*(quanorig-k$x[a])+k$y[a]

  varorig=((1-cutoff)*cutoff)/f^2

  bootlbnorm=c();f=0;k=0
  for (j in 1:boottimes) {
    s=sample(1:length(d),length(d),replace = T)
    bootdata=d[s]
    bootlb=quantile(bootdata,cutoff)
    k=density(bootdata)
    aa=which(k$x<=quanorig)
    a=max(aa)
    b=which.max(k$x>=quanorig)
    f=((k$y[b]-k$y[a])/(k$x[b]-k$x[a]))*(quanorig-k$x[a])+k$y[a]

    v=((1-cutoff)*cutoff)/f^2
    bootlbstand=(bootlb-quantile(d,cutoff))/sqrt(v)
    bootlbnorm[j]=bootlbstand*sqrt(varorig)+quanorig

  }

  lb=quantile(bootlbnorm,cutoff)
  wh=which(d<lb)
  out=data[wh,]
  loc=wh


  p=c()                             #outlier probability
  for (i in wh) {
    p[i]=length(which(bootlbnorm>d[i]))/length(bootlbnorm)
  }

if(ncol(x)==2)
{
    Class=as.factor(ifelse(d<lb,"Outlier","Normal"))
    cols <- c("Outlier" = "red", "Normal" = "blue")

    if(rnames==TRUE)
    {
      s=subset(data,Class=="Outlier")
      gplot=ggplot2::ggplot(data,aes(data[,1],data[,2]))+geom_point(aes(colour=Class,pch=Class))+geom_text(data=s,aes(x=s[,1],y=s[,2],label=rownames(s)),colour="Red", hjust = "inward",check_overlap = T)+ggtitle("Outlier plot using depth")+xlab("Variable1")+ylab("Variable2")+scale_color_manual(values=cols)
      gplot
      }else
    {dd=cbind(data,1:nrow(data))
    s=subset(dd,Class=="Outlier")
    gplot=ggplot2::ggplot(data,aes(data[,1],data[,2]))+geom_point(aes(colour=Class,pch=Class))+geom_text(data=s,aes(x=s[,1],y=s[,2],label=s[,3]),colour="Red", hjust = "inward",check_overlap = T)+ggtitle("Outlier plot using depth")+xlab("Variable1")+ylab("Variable2")+scale_color_manual(values=cols)
    gplot
    }
  l=list("Outlier Observations"=out,"Location of Outlier"=loc,"Outlier Probability"=p[is.na(p)==F],"Scatter plot"=gplot)
}else if(ncol(x)==3)
{
  Class=as.factor(ifelse(d<lb,"Outlier","Usual"))

  plot=plotly::plot_ly(x=data[,1],y=data[,2],z=data[,3],type="scatter3d",mode="markers",color=Class,colors=c("Red","Blue"))

  l=list("Outlier Observations"=out,"Location of Outlier"=loc,"Outlier Probability"=p[is.na(p)==F],"3Dplot"=plot)
}else if(ncol(x)==4)
{
  Class=as.factor(ifelse(d<lb,"Outlier","Usual"))

  plot=plotly::plot_ly(x=data[,1],y=data[,2],z=data[,3],size = data[,4],type="scatter3d",mode="markers",color=Class,colors=c("Red","Blue"))

  l=list("Outlier Observations"=out,"Location of Outlier"=loc,"Outlier Probability"=p[is.na(p)==F],"3Dplot"=plot)
}else
  l=list("Outlier Observations"=out,"Location of Outlier"=loc,"Outlier Probability"=p[is.na(p)==F])
   return(l)
  }
