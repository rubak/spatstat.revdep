### R code from vignette source 'estimate_RACS_properties.Rnw'

###################################################
### code chunk number 1: loadbinarymap
###################################################
library(lacunaritycovariance)
load(system.file("extdata/egbinarymap.RData",
     package="lacunaritycovariance"))
# the following converts egbinarymap to logically-valued pixels
egbinarymap <- as.im(egbinarymap, eps = egbinarymap$xstep*8)
egbinarymap <- eval.im(egbinarymap > 0.5)
plot(egbinarymap,
     col = c("grey", "black"),
     main = "The Binary Map")


###################################################
### code chunk number 2: coverageprobability
###################################################
phat <- coverageprob(egbinarymap)
phat


###################################################
### code chunk number 3: estimatecovariance
###################################################
cvchat <- racscovariance(egbinarymap, estimators = "pickaH",
                         drop = TRUE)
plot(cvchat, main = "Estimated RACS Covariance", axes = TRUE)


###################################################
### code chunk number 4: paircorr
###################################################
pclnest <- eval.im(cvchat /  (phat^2))
plot(pclnest, main = "Pair Correlation Estimate")


###################################################
### code chunk number 5: gblestimation
###################################################
gblest <- gbl(xi = egbinarymap, seq(1, 200/4, by = 1), 
	      estimators = c("GBLcc.pickaH", "GBLemp"))
plot(gblest[[1]], main = "GBL Estimate")


