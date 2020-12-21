#' Univariate Outlier Detection(Intersection of all the methods)
#'
#' Takes a vector and finds its outliers using combination of different methods
#' @param x vector for which outliers are to be found
#' @param k No. of nearest neighbours to be used for distance methods, default value is 0.05*nrow(x)
#' @param cutoff Percentile threshold used for outlier detection using bootstrapping, default value is 0.95
#' @param Method Distance method, default is euclidean
#' @param rnames Logical value indicating whether the dataset has rownames, default value is False
#' @param depth  Logical value indicating whether depth based method should be used or not, default is False
#' @param dens  Logical value indicating whether density based method should be used or not, default is False
#' @param dist  Logical value indicating whether distance based methods should be used or not, default is False
#' @details UnivariateOutlierDetection finds outlier observations for an univariate data using different methods and based on all the methods, labels an observation as outlier(intersection of all the methods). It also shows the scatterplot of the data with labelled outliers with observation no. as x-axis.
#' @return Outlier Observations: A vector of outlier observations
#' @return Location of Outlier: Vector of Sr. no. of outliers
#' @author Vinay Tiwari, Akanksha Kashikar

#' @examples
#' #Create dataset
#' X=iris[,1:4]
#' #Outlier detection
#' depthout(X,cutoff=0.05)
#' UnivariateOutlierDetection(iris[,1],cutoff=.95,Method="euclidean",rnames=FALSE)

  UnivariateOutlierDetection=function(x,k=0.05*length(x),cutoff=.95,dist=FALSE,dens=FALSE,depth=FALSE,Method="euclidean",rnames=FALSE)
  {
    data=x
  out1=mahauni(x)
  out2=dispuni(x)
  out=intersect(out1,out2)

  if(dist==FALSE)
  {out=out
  }else
  {
    out3=nnuni(x,Method=Method,k,cutoff);out4=nnkuni(x,Method=Method,k,cutoff)
    out=intersect(intersect(out3,out),intersect(out4,out))

  }
  if(depth==FALSE)
  {out=out
  }else
  {
    out5=depthuni(x,cutoff)
    out=intersect(out5,out)
  }
  if(dens==FALSE)
  {out=out
  }else
  {
    out6=densuni(x,k=5,cutoff)
    out=intersect(out6,out)
  }


    d=1:length(data)
    Class=c()
    Class=rep("Normal",length(d))
    for (i in 1:length(out)) {

      Class[d==out[i]]="Outlier"
    }
    cols <- c("Outlier" = "red", "Normal" = "blue")
    data=as.data.frame(cbind(1:length(x),x))
    if(rnames==T)
    {
      s=subset(data,Class=="Outlier")
      gplot=ggplot2::ggplot(data,aes(data[,1],data[,2]))+ggplot2::geom_point( aes(colour=Class,pch=Class))+ggplot2::geom_text(data=s,aes(x=s[,1],y=s[,2],label=rownames(s)),colour="Red", hjust = "inward",check_overlap = T)+ggplot2::ggtitle("Outlier plot")+ggplot2::xlab("Sr.No.")+ggplot2::ylab("Variable")+ggplot2::scale_color_manual(values=cols)

    }else{
      dd=cbind(data,1:nrow(data))
      s=subset(dd,Class=="Outlier")
      gplot=ggplot2::ggplot(data,aes(data[,1],data[,2]))+ggplot2::geom_point(aes(colour=Class,pch=Class))+ggplot2::geom_text(data=s,aes(x=s[,1],y=s[,2],label=s[,3]),colour="Red", hjust = "inward",check_overlap = T)+ggplot2::ggtitle("Outlier plot")+ggplot2::xlab("Sr.No.")+ggplot2::ylab("Variable")+ggplot2::scale_color_manual(values=cols)
    }
    Out=x[out]
    l=list("Outlier Observations"=Out,"Location of Outlier"=out,"Scatter plot"=gplot)
  return(l)
}

  mahauni=function(x,cutoff=0.95)
  {
    s=scale(x)
    cut=qnorm(cutoff)
    wh=which(abs(s)>cut)
    out=x[wh]
    loc1=wh
    return(loc1)
  }

  nnkuni=function(x,k=0.05*nrow(x),cutoff=0.95,Method="euclidean",rnames=F)
  {
    data=x
    dis=as.matrix(dist(data,diag=T,upper = T,method=Method))
    d=c()
    for (i in 1:length(data)) {
      temp=dis[,i]
      d[i]=sort(temp)[k]
    }
    bootub=c()
    for (j in 1:1000) {
      s=sample(1:length(d),length(d),replace = T)
      bootdata=d[s]
      bootub[j]=quantile(bootdata,cutoff)

    }
    ub=mean(bootub)
    loc=which(d>ub)
    return(loc)
  }

  nnuni=function(x,k=0.05*nrow(x),cutoff=0.95,Method="euclidean",rnames=F)
  {
    data=x
    dis=as.matrix(dist(data,diag=T,upper = T,method=Method))
    d=c()
    for (i in 1:length(data)) {
      temp=dis[,i]
      neighbour=sort(temp)[1:k]
      d[i]=mean(neighbour)
    }
    bootub=c()
    for (j in 1:1000) {
      s=sample(1:length(d),length(d),replace = T)
      bootdata=d[s]
      bootub[j]=quantile(bootdata,cutoff)

    }
    ub=mean(bootub)
    loc=which(d>ub)
    return(loc)
  }

  dispuni=function(x,cutoff=0.95,rnames=F)
  {
    data=x;d=c()
    for (i in 1:length(data)) {
      remdata=data[-i]
      remvar=var(remdata)
      d[i]=abs((remvar)-(var(x)))
    }
    bootub=c()
    for (j in 1:1000) {
      s=sample(1:length(d),length(d),replace = T)
      bootdata=d[s]
      bootub[j]=quantile(bootdata,0.95)

    }
    ub=mean(bootub)
    loc=which(d>ub)
    return(loc)
  }

  densuni=function(x,k=10,cutoff=0.95,rnames=F)
  {
    data=x
    t=ldbod::ldbod(x,k)
    d=t$rkof
    bootub=c()
    for (j in 1:1000) {
      s=sample(1:length(d),length(d),replace = T)
      bootdata=d[s]
      bootub[j]=quantile(bootdata,cutoff)

    }
    ub=mean(bootub)
    loc=which(d>ub)
    return(loc)
  }

  depthuni=function(x,rnames=F,cutoff=0.05)
  {
    data=as.matrix(x)
    d=c()
    for (i in 1:length(x)) {
      d[i]=depth::depth(x[i],x)

    }
    bootlb=c()
    for (j in 1:1000) {
      s=sample(1:length(d),length(d),replace = T)
      bootdata=d[s]
      bootlb[j]=quantile(bootdata,cutoff)

    }
    lb=mean(bootlb)
    loc=which(d<lb)
    return(loc)
  }



