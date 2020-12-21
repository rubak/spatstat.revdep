# defast = function(x){
#   if(x>0.1){ast = " "}else
#   {
#     if(x>0.05){ast = "."}else
#     {
#       if(x>0.01){ast = "*"}else
#       {
#         if(x>0.001){ast = "**"}else
#         {/
#           {
#             ast = "***"
#           }
#         }
#       }
#     }
#   }
#   return(ast)
# }


########################################################################
#EM ALGORITHM
########################################################################
# 
# data(ais)
# attach(ais)
# 
# ##Setting
# y<-BMI
# x<-cbind(1,LBM,Sex)
# 
# cc = rep(0,length(y))
# LI = LS = rep(NA,length(y))
# ind = sample(x = c(0,1),size = length(y),replace = TRUE,prob = c(0.95,0.05))
# ind1 = (ind == 1)
# 
# cc[ind1] = 1
# LI[ind1] = y[ind1] - 10
# LS[ind1] = y[ind1] + 10
# 
# y[ind1] = NA
# y
# cc
# LI
# LS


cens.lqr <- function(y,x,cc,LL,UL,p=0.5,nu = NULL,precision = 1e-6,envelope=FALSE){
  
  LI = LL
  LS = UL
  
  pi = 3.14159
  
  loglikT = function(x,mu,sigma,nu,p)
  {
    return(sum(log(
      ifelse(test=x<mu,
             yes=(4*p*(1-p)*gamma((nu+1)/2)/(gamma(nu/2)*sigma*sqrt(nu*pi)))*((4*((x-mu)^2)/(nu*sigma^2))*(1-p)^2 +1)^(-(nu+1)/2),
             no=(4*p*(1-p)*gamma((nu+1)/2)/(gamma(nu/2)*sigma*sqrt(nu*pi)))*((4*((x-mu)^2)/(nu*sigma^2))*(p)^2 +1)^(-(nu+1)/2))
    )))
  }
  
  likefun = function(nu,...){
    loglik = loglikT(y[!ind1],as.numeric(x[!ind1,]%*%teta[1:d]),sigma = teta[d+1],nu=nu,p=p)
    if(sum(cc)>0){
    Lvec   = vecpSKT(a = LI[ind1],b = LS[ind1],mu = as.numeric(x[ind1,]%*%beta),sigma = sqrt(sigma2),nu = nu,p = p)
    }else{
      Lvec = 1
    }
    return(loglik + sum(log(Lvec)))
  }
  
  # Initial -----------------------------------------------------------------
  
  
  n <- length(y)
  d = dim(x)[2]
  ind1 = (cc == 1)
  #beta<-solve(t(x[!ind1,])%*%x[!ind1,])%*%t(x[!ind1,])%*%y[!ind1]
  beta <- suppressWarnings(quantreg::rq(formula = y ~ x - 1,tau = p,na.action = "na.omit")$coefficients)
  sigma2<-sum((y[!ind1]-x[!ind1,]%*%beta)^2)/(n-d-sum(cc))
  teta <- c(beta,sigma2)
  
  beta.old <- beta
  sigma2.old <- sigma2
  nu.old = nu
  if(is.null(nu.old)){nu = 4}
  
  
  # preparation -------------------------------------------------------------
  
  
  EU = rep(1,n)
  pvec = rep(1,n)
  criterio <- 1
  like.old = -1
  count <- 0
  
  yest   = y
  ypred  = y
  ypred2 = ypred^2
  if(sum(cc)>0){
    yest[ind1] = vecmeanSKT(LI[ind1],LS[ind1],as.numeric(x[ind1,]%*%beta),sqrt(sigma2),nu,p)
  }
  xi = rep(NA,n)
  
  # Lvec  = vecpSKT(a = LI[ind1],b = LS[ind1],mu = as.numeric(x[ind1,]%*%beta),sigma = sqrt(sigma2),nu = nu,p = p)
  # Lvec2 = matrix(apply(vecexpSKT(LI[ind1],LS[ind1],as.numeric(x[ind1,]%*%beta),sqrt(nnu2*sigma2),nu + 2,p),2,unlist),ncol = 3,byrow = TRUE)
  
  while(criterio > precision){
    
    count <- count + 1
    
    #print(count)
    
    # EM algorithm ------------------------------------------------------------
    
    ### E-step: calculando ui, tui###
    
    izq = yest - as.numeric(x%*%beta) <= 0
    
    xi[izq] = 1-p
    xi[!izq] = p
    zi = as.numeric((yest - x%*%beta)/sqrt(sigma2))
    
    EU = (nu + 1)/(nu + 4*xi^2*zi^2)
    
    if(sum(cc)>0){
    nnu2  = nu/(nu + 2)
    Lvec  = vecpSKT(a = LI[ind1],b = LS[ind1],mu = as.numeric(x[ind1,]%*%beta),sigma = sqrt(sigma2),nu = nu,p = p)
    Lvec2 = matrix(apply(vecexpSKT(LI[ind1],LS[ind1],as.numeric(x[ind1,]%*%beta),sqrt(nnu2*sigma2),nu + 2,p),2,unlist),ncol = 3,byrow = TRUE)
    
    EU[ind1] = Lvec2[,1]/Lvec
    
    ypred[ind1] = Lvec2[,2]
    ypred2[ind1] = Lvec2[,3]
    }
    
    #M-step
    
    pvec = xi^2
    
    Omega = diag(c(EU*pvec))
    den = t(x)%*%Omega%*%x
    num = t(x)%*%Omega%*%ypred
    beta = solve(den)%*%(num)
    sigma2 = 4/n*sum(diag((Omega%*%(diag(ypred2) - (2*ypred - x%*%beta)%*%t(x%*%beta)))))
    
    param <- teta
    teta <- c(beta,sqrt(sigma2))
    
    # Likelihood --------------------------------------------------------------
    
    
    #like = likefun(nu)
    
    
    # nu estimation -----------------------------------------------------------
    
    # if(is.null(nu.old)){
    #   likeseq = sapply(X = seq(1,30),FUN = likefun)
    #   like = max(likeseq)
    #   nu = which.max(likeseq)
    # }
    
    if(is.null(nu.old)){
      l10 = likefun(10)
      l11 = likefun(11)
      if(l10-l11>0){
        
        grid0 = 1:9
        
        likeseq = c(sapply(X = grid0,FUN = likefun),l10)
        like = max(likeseq)
        nu = which.max(likeseq)
        
      }else{
        
        l20 = likefun(20)
        l21 = likefun(21)
        if(l20-l21>0){
          
          grid0 = 11:19
          
          likeseq = c(sapply(X = grid0,FUN = likefun),l20)
          like = max(likeseq)
          nu = which.max(likeseq) + 10
          
        }else{
          
          l30 = likefun(30)
          l31 = likefun(31)
          if(l30-l31>0){
            
            grid0 = 21:29
            
            likeseq = c(sapply(X = grid0,FUN = likefun),l30)
            like = max(likeseq)
            nu = which.max(likeseq) + 20
            
          }else{
            
            grid0 = c(seq(30,60,by = 5),seq(70,100,by = 10),200,300)
            
            likeseq = sapply(X = grid0,FUN = likefun)
            like = max(likeseq)
            nu = grid0[which.max(likeseq)]
          }
        }
      }
    }else{
      
      like = likefun(nu)
      
    }
    
    # end EM ------------------------------------------------------------------
    
    
    # convergence? ------------------------------------------------------------
    
    criterio <- abs(like/like.old - 1)
    
    like.old = like
    
    #print(criterio)
  }
  
  sigma = sqrt(sigma2)
  
  # Decimal nu value --------------------------------------------------------
  
  if(is.null(nu.old)){
    optnu  = optimize(f = likefun,interval = c(nu - 0.7,nu + 0.7),maximum = TRUE)
    like = optnu$objective #-380.0496
    nu     = round(optnu$maximum,2)
  }
  
  # loglik based criteria ---------------------------------------------------
  
  npar = length(c(teta)) + 1
  AIC  = -2*like +2*npar
  BIC  = -2*like +log(n)*npar
  HQ   = -2*like +2*log(log(n))*npar
  
  
  # estimated values --------------------------------------------------------
  
  if(sum(cc)>0){
  yest = vecmeanSKT(LI[ind1],LS[ind1],as.numeric(x[ind1,]%*%beta),sqrt(sigma2),nu,p)
  }
  
  #fitted and residuals values
  fitted.values    = c(x%*%teta[1:d])
  imputed.values   = y
  imputed.values[ind1]   = yest
  residuals        = imputed.values - fitted.values
  #hist(residuals)
  
  ###### Computing SE's
  
  #Individual score (by i)
  derbetai  = (4/sigma2)*diag(c(EU*pvec*(ypred - x%*%beta)))%*%x
  dersigmai = - 1/sigma + c(4/sigma^3)*diag(Omega%*%(diag(ypred2) - (2*ypred - x%*%beta)%*%t(x%*%beta)))
  
  gradi     = cbind(derbetai,dersigmai)
  
  #apply(gradi,2,sum)
  
  IEI = 0
  for(i in 1:n)
  {
    IEI = gradi[i,]%*%t(gradi[i,]) + IEI
  }
  
  SE         = sqrt(diag(solve(IEI)))
  
  table      = data.frame(beta,SE[1:d],beta/SE[1:d],2*pnorm(abs(beta/SE[1:d]),lower.tail = F))
  asteriscos = apply(X = table[4],MARGIN = 1,FUN = defast)
  table      = data.frame(round(table,5),asteriscos)
  rownames(table) = paste("beta",1:d)
  colnames(table) = c("Estimate","Std. Error","z value","Pr(>|z|)","")
  # 
  
  # ######## ENVELOPES: T
  # 
  if(envelope==TRUE){
    muc<- (imputed.values -x%*%beta)
    #muc<- (y -x%*%beta)
    #muc<- (y - y1)
    Ind<- (muc<0)+0
    d2s<- (muc/sigma)*(p-Ind)
    d2s=sort(d2s)
    
    # hist(d2s,breaks = 50,freq = F)
    # seqq2 = seq(0,4,length.out = 100)
    # dens3 = dent(x = seqq2,mu = 0,sigma2 = 1/4,nu = nu)*2
    # lines(seqq2,dens3,col="blue",lwd=2)
    
    qseq = seq(1/n,0.99,length.out = n)
    xq2 <- qt(qseq/2 + 0.5,df = nu)*0.5
    
    #plot(d2s,xq2);abline(coef = c(0,1))
    
    MM = 500
    
    Xsim<-matrix(0,MM,n)
    for(i in 1:MM){
      Xsim[i,]<-abs(rt(n,df = nu))*0.5
    }
    
    Xsim2<-apply(Xsim,1,sort)
    d21<-matrix(0,n,1)
    d22<-matrix(0,n,1)
    for(i in 1:n){
      d21[i]  <- quantile(Xsim2[i,],0.05)
      d22[i]  <- quantile(Xsim2[i,],0.95)
    }
    
    d2med <-apply(Xsim2,1,median)
    
    fy <- range(d2s,d21,d22)
    
    exp2 = bquote("Theoretical HT(0.5," * nu == .(nu) * ") quantiles")
    
    plot(xq2,d2s,xlab = exp2,
         ylim = fy,ylab="Sample values and envelope",pch=20,cex=0.8,type="n",
         main = paste("Student's t censored QR model (p = ",p,")",sep=""))
    
    polygon(c(rev(xq2), xq2), c(rev(d22), d21), col = 'grey80', border = NA)
    
    pchs = rep(1,n)
    pchs[ind1] = 4
    
    colors = rep(1,n)
    colors[ind1] = 4
    
    points(xq2,d2s,pch=pchs,cex=1,col = colors)
    
    lines(xq2,d21,type="l",ylim=fy,xlab="",ylab="",lwd=0.5)
    lines(xq2,d2med,type="l",ylim=fy,xlab="",ylab="")
    lines(xq2,d22,type="l",ylim=fy,xlab="",ylab="",lwd=0.5)
  }
  
  return(list(iter =count,
              criteria = criterio,
              #theta = teta,
              beta = c(beta),
              sigma = sqrt(sigma2),
              nu = nu,
              SE = SE,
              table = table,
              loglik = like,
              AIC=AIC,BIC=BIC,HQ=HQ,
              fitted.values = fitted.values,
              imputed.values = imputed.values,
              residuals = residuals))
}
# 
# out = cens.lqr(y,x,cc,LI,LS,p=0.5,nu = NULL,precision = 1e-6,envelope = TRUE)
# 
# out = cens.lqr(y,x,cc,LI,LS,p=0.5,nu = 8,precision = 1e-6,envelope = TRUE)
