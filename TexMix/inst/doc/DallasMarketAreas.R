## ----setup, include = FALSE---------------------------------------------------
knitr::opts_chunk$set(collapse = TRUE, comment = "R>", cache=FALSE, fig.width=7.2)
rm(list=ls())                                        # start fresh
library(TexMix)              # Get data and support function plotBoxesByFactor.r
library(ClustGeo)                  # Hierarchical Cluster Analysis with Spatial Constraint
library(maptools)
library(sp)
library(DMwR)
library(spdep)

## ----Adjacency, echo=TRUE, fig.height=7.2-------------------------------------
data(tractShp)
tractShp <- tractShp[tractShp$NIGHTPOP!=0, ]         # Remove 2 airport tracts with NA's

##
## Identify spatial neighboring census tracts and process them
##
nb <- spdep::poly2nb(tractShp, queen=F)              # extract first order neighbors links
B <- spdep::nb2mat(nb, style="B")                    # convert neighbor list to binary matrix

## Visualize first order neighbors
plot(tractShp, col="palegreen3", border=grey(0.9), axes=T) 
plot(nb, coords=coordinates(tractShp), pch=19, cex=0.1, col="blue", add=T)
title("Spatial Neighbors Links among Tracts") 

## Transform the contiguity matrix into spatial dissimilarity matrix
geoDist <- 1-B; 
diag(geoDist) <- 0; 
geoDist <- as.dist(geoDist)

##
## Alternative spatial dissimilarity matrices
##
# # Calculate shortest path distances between the tracts
# BNa <- ifelse(B==0,NA,B)                             # recode 0's to NA's
# allPath <- e1071::allShortestPaths(BNa)              # calucate the shortest path among all nodes
# topoDist <- as.dist(allPath$length)                  # number of steps from origin to destination node
# e1071::extractPath(allPath, 1, 527)                  # example: path between tract 1 and tract 527
#
# # spherical distances matrix among tracts in km
# sphDist <- as.dist(sp::spDists(tractShp, longlat=T)) 

## ----DataPrep, echo=TRUE, fig.height=7.2, fig.width=7.2-----------------------
##
## Select variables for feature dissimilarity organized into three groups:
##   [a] demographic features
##   [b] socio-economic features
##   [c] housing infrastructure features
##

## Calculate additional features
tractShp$PRE1960 <- tractShp$PCTB1950+tractShp$PCTB1940+tractShp$PCTBPRE
tractShp$POST2000 <- tractShp$PCTB2000+tractShp$PCTB2010

## Feature list
varKeep <- c("PCTWHITE","PCTBLACK","PCTHISPAN","MEDAGE",
             "PCTUNIVDEG","HHMEDINC","PCTUNEMP","PCTNOHINS",
             "POPDEN","PCTDAYPOP","PCTHUVAC","MEDVALHOME","PRE1960","POST2000")
xVars <- tractShp@data
xVars <- xVars[varKeep]
summary(xVars)

## Replace missing values in MEDVALHOME by neighbors average
xVars <- DMwR::knnImputation(xVars, k=5, scale=T, meth="weighAve")
summary(xVars$MEDVALHOME)

## ID each row of xVars
row.names(xVars) <- 1:nrow(xVars)

## Calculate feature distance matrix
featDist <- dist(scale(xVars))

## ----AlphaMixture, echo=TRUE, fig.height=7.2----------------------------------
##
## Evaluate mixture of feature and spatial dissimilarity.
##
K <- 12                            # Number of distinct clusters
range.alpha <- seq(0, 1, by=0.1)   # Evaluation range
cr <- choicealpha(featDist, geoDist, range.alpha, K, graph=TRUE)

## ----ClusterAnalysis, echo=TRUE, fig.height=7.2, fig.width=7.2----------------
## 
## Perform spatially constrained cluster analysis
##
tree <- hclustgeo(featDist, geoDist, alpha=0.2)
plot(tree, hang=-1)
rect.hclust(tree, k=K)

## ----TractCount, echo=TRUE, fig.height=7.2------------------------------------
## Number of census tracts per market area
neighClus <- as.factor(cutree(tree, K))        # Determine cluster membership
table(neighClus)                               # number of tracts in each cluster

## ----AllMap, echo=TRUE, fig.height=7.2----------------------------------------
##
## Map Results
##
mapColorQual(neighClus, tractShp, 
             map.title ="Spatially constrained Cluster Analysis",
             legend.title="Cluster\nId.", legend.cex=0.9)
plot(lakesShp, col="skyblue", border="skyblue",add=T)
plot(hwyShp, col="cornsilk2", lwd=4, add=T)

## ----AreaCharacteristics, echo=TRUE, fig.height=7.2, fig.width=7.2------------
##
## Evaluate cluster characteristics
##
plotBoxesByFactor(xVars[,1:7], neighClus, ncol=2, zTrans=T, varwidth=F)
plotBoxesByFactor(xVars[,8:14], neighClus, ncol=2, zTrans=T, varwidth=F)

## ----SelectMap, echo=TRUE, fig.height=7.2-------------------------------------
## Select subset nLevels[c(4,10)] and set remaining clusters to NA
nLevels <- levels(neighClus)
neighClusSet <- as.factor(ifelse(neighClus %in% nLevels[c(4,10)], neighClus, NA))
mapColorQual(neighClusSet, tractShp, 
             map.title ="Clusters of Major Employment Centers",
             legend.title="Cluster\nId.", legend.cex=0.9)
plot(lakesShp, col="skyblue", border="skyblue",add=T)
plot(hwyShp, col="cornsilk2", lwd=4, add=T)

