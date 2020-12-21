best.lqr.range <- function(y, x, p.min = 0.10, p.max = 0.9, by = 0.1, criterion = "AIC")
{
 start.time     <- Sys.time()
 if(p.min > p.max) {stop("The number p.min should be less than p.max.\n")}
 if(p.min == p.max) {stop("The number p.min should be different than p.max.\n")}
 if(criterion != "" && criterion != "AIC" && criterion != "BIC" && criterion != "HQ" && criterion != "loglik"){stop("The criterion must be AIC, BIC, HQ or loglik.")}

 p       <- seq(p.min, p.max, by = by)
 vdist   <- c("normal", "t", "laplace", "slash", "cont")
 vDIST   <- c("Normal", "Student-t", "Laplace", "Slash", "Cont. Normal")
 
 obj.out1 <- lqr(y,x,p = seq(from = p.min,to = p.max,by = by),dist = vdist[1])
 obj.out2 <- lqr(y,x,p = seq(from = p.min,to = p.max,by = by),dist = vdist[2])
 obj.out3 <- lqr(y,x,p = seq(from = p.min,to = p.max,by = by),dist = vdist[3])
 obj.out4 <- lqr(y,x,p = seq(from = p.min,to = p.max,by = by),dist = vdist[4])
 obj.out5 <- lqr(y,x,p = seq(from = p.min,to = p.max,by = by),dist = vdist[5])
 
 RES1 <- matrix(data = NA, nrow = 4, ncol = length(p)); for (k in 1:length(p)){RES1[1, k] = (obj.out1[[k]])$AIC; RES1[2, k] = (obj.out1[[k]])$BIC; RES1[3, k] = (obj.out1[[k]])$HQ; RES1[4, k] = (obj.out1[[k]])$loglik}
 RES2 <- matrix(data = NA, nrow = 4, ncol = length(p)); for (k in 1:length(p)){RES2[1, k] = (obj.out2[[k]])$AIC; RES2[2, k] = (obj.out2[[k]])$BIC; RES2[3, k] = (obj.out2[[k]])$HQ; RES2[4, k] = (obj.out2[[k]])$loglik}
 RES3 <- matrix(data = NA, nrow = 4, ncol = length(p)); for (k in 1:length(p)){RES3[1, k] = (obj.out3[[k]])$AIC; RES3[2, k] = (obj.out3[[k]])$BIC; RES3[3, k] = (obj.out3[[k]])$HQ; RES3[4, k] = (obj.out3[[k]])$loglik}
 RES4 <- matrix(data = NA, nrow = 4, ncol = length(p)); for (k in 1:length(p)){RES4[1, k] = (obj.out4[[k]])$AIC; RES4[2, k] = (obj.out4[[k]])$BIC; RES4[3, k] = (obj.out4[[k]])$HQ; RES4[4, k] = (obj.out4[[k]])$loglik}
 RES5 <- matrix(data = NA, nrow = 4, ncol = length(p)); for (k in 1:length(p)){RES5[1, k] = (obj.out5[[k]])$AIC; RES5[2, k] = (obj.out5[[k]])$BIC; RES5[3, k] = (obj.out5[[k]])$HQ; RES5[4, k] = (obj.out5[[k]])$loglik}
 
 p.qr <- paste("p=", p.min, sep = "")
 for (i in 1:length(p)) {
   k    <- paste("p=", p[i], sep = "")
   p.qr <- c(p.qr, k)
 }
 p.qr   <- p.qr[-1]
 colnames(RES1)=colnames(RES2)=colnames(RES3) = colnames(RES4) = colnames(RES5)= p.qr
 rownames(RES1)=rownames(RES2)=rownames(RES3) = rownames(RES4) = rownames(RES5)= c("AIC", "BIC", "HQ", "loglik")
 RES   <- rbind(RES1,RES2,RES3,RES4,RES5)
 
 if (criterion == "AIC"){index = which(RES[c(1,5,9,13,17), ] == min(RES[c(1,5,9,13,17), ]))}
 if (criterion == "BIC"){index = which(RES[c(2,6,10,14,18), ] == min(RES[c(2,6,10,14,18), ]))}
 if (criterion == "HQ") {index = which(RES[c(3,7,11,15,19), ] == min(RES[c(3,7,11,15,19), ]))}
 if (criterion == "loglik") {index = which(RES[c(4,8,12,16,20), ] == max(RES[c(4,8,12,16,20), ]))}
 
 p.best <- p[floor(index/5) + 1]
 dist   <- vdist[index - index%/%5*5]
 
 if((index - index%/%5*5) == 1){  obj.out <- obj.out1; RES=RES1}
 if((index - index%/%5*5) == 2){  obj.out <- obj.out2; RES=RES2}
 if((index - index%/%5*5) == 3){  obj.out <- obj.out3; RES=RES3}
 if((index - index%/%5*5) == 4){  obj.out <- obj.out4; RES=RES4}
 if((index - index%/%5*5) == 5){  obj.out <- obj.out5; RES=RES5}
 
 out    <- obj.out[[floor(index/5) + 1]]
 end.time      <- Sys.time()
 time.taken    <- end.time - start.time
 
 cat("\n")
 cat("--------------------------------------------------------------\n")
 cat("        Quantile Linear Regression using SKD family\n")
 cat("--------------------------------------------------------------\n")
 cat("\n")
 cat("Criterion =", criterion, "\n")
 cat("Best fit =", vDIST[index - index%/%5*5], "\n")
 cat("Quantile =", p.best, "\n")
 cat("\n")
 cat("--------------------------------\n")
 cat("Model Likelihood-Based criterion\n")
 cat("--------------------------------\n")
 cat("\n")
 print(RES)
 cat("\n")
 cat("---------\n")
 cat("Estimates\n")
 cat("---------\n")
 cat("\n")
 print(out$table)
 cat("---\n")
 cat("Signif. codes:  0 \"***\" 0.001 \"**\" 0.01 \"*\" 0.05 \".\" 0.1 \" \" 1\n")
 cat("\n")
 cat("sigma =", round(out$sigma, 5), 
     "\n")            
 if (dist == "normal" || dist == "laplace") {cat("\n")}
 if (dist == "t" || dist == "slash") {
   cat("nu    =", out$nu, "\n")
   cat("\n")
 }
 if (dist == "cont") {
   cat("nu    =", out$nu, "\n")
   cat("gamma =", out$gamma, "\n")
   cat("\n")
 }
 cat("Processing time =",time.taken,units(time.taken))
 cat('\n')
 if (dist == "normal" || dist == "laplace") {
   obj.out = list(iter = out$iter, criteria = out$criteria, 
                  beta = out$beta, sigma = out$sigma, SE = out$SE, table = out$table, loglik = out$loglik, 
                  AIC = out$AIC, BIC = out$BIC, HQ = out$HQ, fitted.values = out$fitted.values, 
                  residuals = out$residuals, time=time.taken, dist=dist, p= p.best)
 }
 if (dist == "t" || dist == "slash")
 {
   obj.out = list(iter = out$iter, criteria = out$criteria, beta = out$beta, sigma = out$sigma, nu = out$nu, SE = out$SE, table = out$table,loglik = out$loglik, AIC = out$AIC, BIC = out$BIC, 
              HQ = out$HQ, fitted.values = out$fitted.values, residuals = out$residuals , time=time.taken, dist=dist, p= p.best)
 }
 if (dist == "cont") {
   obj.out = list(iter = out$iter, criteria = out$criteria, 
                  beta = out$beta, sigma = out$sigma, nu = out$nu, gamma = out$gamma, SE = out$SE, 
                  table = out$table, loglik = out$loglik, AIC = out$AIC, 
                  BIC = out$BIC, HQ = out$HQ, fitted.values = out$fitted.values, 
                  residuals = out$residuals, time=time.taken, dist=dist, p= p.best)
 }
 class(obj.out) = "qr"
 return(obj.out)
 
}

#library(graphics)
#library(stats) 
#library(ghyp)
#library(spatstat)
#source("/home/lbenitesanchez/Dropbox/LQR package/Package codes/lqr/R/EM2.R")
#source("/home/lbenitesanchez/Dropbox/LQR package/Package codes/lqr/R/FUNCTION.R")
#source("/home/lbenitesanchez/Dropbox/LQR package/Package codes/lqr/R/auxiliar.R")
#source("/home/lbenitesanchez/Dropbox/LQR package/Package codes/lqr/R/pdqr.R")
#best_p =  best.lqrfull(y, x, p.min = 0.4, p.max = 0.5, by = 0.1, criterion = "AIC")
#model = best_p


