# require(MomTrunc)
# 
# mu = 0
# sigma = 1
# nu = 4
# p = 0.1
# a = 1
# b = 2

expSKT = function(a,b,mu,sigma,nu,p){
  
  sd1 = sqrt(sigma^2/(1-p)^2/4)
  sd2 = sqrt(sigma^2/p^2/4)
  
  a1 = a
  b2 = b
  
  b1 = min(mu,b)
  a2 = max(mu,a)
  
  # Probabilities -----------------------------------------------------------
  
  L1 = ifelse(a1<b1,pmvnormt(lower = a1,upper = b1,mean = mu,sigma = sd1^2,nu = nu),0)
  L2 = ifelse(a2<b2,pmvnormt(lower = a2,upper = b2,mean = mu,sigma = sd2^2,nu = nu),0)
  
  L = 2*p*L1 + 2*(1-p)*L2
  # Expectations ------------------------------------------------------------
  
  if(a1<b1 & L1 > 0){
    E1 = meanvarTMD(lower = a1,upper = b1,mu = mu,Sigma = sd1^2,nu = nu,dist = "t")
  }else{
    E1 = list(mean = 0,EYY = 0,varcov = 0)
  }
  
  if(a2<b2 & L2 > 0){
    E2 = meanvarTMD(lower = a2,upper = b2,mu = mu,Sigma = sd2^2,nu = nu,dist = "t")
  }else{
    E2 = list(mean = 0,EYY = 0,varcov = 0)
  }
  E = 2/L*(p*L1*E1$mean + (1-p)*L2*E2$mean)
  
  Equad = 2/L*(p*L1*E1$EYY + (1-p)*L2*E2$EYY)
  
  return(list(L = L, E1 = E, E2 = Equad))
  
}

meanSKT = function(a,b,mu,sigma,nu,p){
  
  sd1 = sqrt(sigma^2/(1-p)^2/4)
  sd2 = sqrt(sigma^2/p^2/4)
  
  a1 = a
  b2 = b
  
  b1 = min(mu,b)
  a2 = max(mu,a)
  
  # Probabilities -----------------------------------------------------------
  
  L1 = ifelse(a1<b1,pmvnormt(lower = a1,upper = b1,mean = mu,sigma = sd1^2,nu = nu),0)
  L2 = ifelse(a2<b2,pmvnormt(lower = a2,upper = b2,mean = mu,sigma = sd2^2,nu = nu),0)
  
  L = 2*p*L1 + 2*(1-p)*L2
  # Expectations ------------------------------------------------------------
  
  if(a1<b1 & L1 > 0){
    E1 = onlymeanTMD(lower = a1,upper = b1,mu = mu,Sigma = sd1^2,nu = nu,dist = "t")
  }else{
    E1 = 0
  }
  
  if(a2<b2 & L2 > 0){
    E2 = onlymeanTMD(lower = a2,upper = b2,mu = mu,Sigma = sd2^2,nu = nu,dist = "t")
  }else{
    E2 = 0
  }
  E = 2/L*(p*L1*E1 + (1-p)*L2*E2)
  
  return(E)
  
}


pSKT = function(a,b,mu,sigma,nu,p){
  
  sd1 = sqrt(sigma^2/(1-p)^2/4)
  sd2 = sqrt(sigma^2/p^2/4)
  
  a1 = a
  b2 = b
  
  b1 = min(mu,b)
  a2 = max(mu,a)
  
  # Probabilities -----------------------------------------------------------
  
  L1 = ifelse(a1<b1,pmvnormt(lower = a1,upper = b1,mean = mu,sigma = sd1^2,nu = nu),0)
  L2 = ifelse(a2<b2,pmvnormt(lower = a2,upper = b2,mean = mu,sigma = sd2^2,nu = nu),0)
  
  L = 2*p*L1 + 2*(1-p)*L2
  return(L)
}

# pSKT(-1,2,0,1,4,0.1)
# expSKT(-1,2,0,1,4,0.1)

meanSKT(-1,2,0,1,4,0.1)

vecpSKT = Vectorize(pSKT)
vecexpSKT = Vectorize(expSKT)
vecmeanSKT = Vectorize(meanSKT)
