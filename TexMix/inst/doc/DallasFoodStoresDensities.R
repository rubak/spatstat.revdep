## ----setup, include=FALSE-----------------------------------------------------
knitr::opts_chunk$set(collapse = TRUE, comment = "R>", cache=FALSE, fig.width=7.2)
rm(list=ls())
library(TexMix)
library(spatstat.core)        # Key library for spatial point pattern analysis
library(smacpod)         # Relative risk kernel densities based on statstat
library(colorspace)      # Used for thematic mapping of kernel densities
library(car)             # Regression diagnostics
library(effects)         # Effect plots to explore non-linearities
library(sp)
library(maptools)
library(rgdal)

## ----weightedKernel, echo=TRUE, fig.height=5----------------------------------
###################################################
### Impact of weights on Gaussian Kernel Density
###################################################
x <- c(-5,4.9,5,5.1)
nx <-length(x)
bw <-1
fx <- density(x , bw=bw, kernel="gaussian")
plot(fx,
     main="Impact on Weights on Gaussian Kernel Densities\nweight=1 versus weight=3 (equivalent to 3 points)")
abline(h=0)
points(x, rep(0, nx), col="red", cex=1.5, pch=20)
abline(v=c(-5-bw,-5+bw, 5-bw,5+bw), lty=2)

## ----LogRiskRatio, echo=TRUE, fig.height=5------------------------------------
######################################################
### Grovery versus Convenience Stores Kernel Densities
######################################################
bw <- 1
groc <- c(-1,-1.5,-2)
conv <- c(1,1.5,2)
both <- c(groc,conv)

f.both <- density(both , bw=bw, kernel="gaussian")
f.both$y <- f.both$y+0.006
f.groc <- density(groc , bw=bw, kernel="gaussian",
                  from=min(f.both$x), to=max(f.both$x))
f.groc$y <- f.groc$y/2
f.conv <- density(conv , bw=bw, kernel="gaussian",
                  from=min(f.both$x), to=max(f.both$x))
f.conv$y <- f.conv$y/2
f.risk <- log(f.conv$y/f.groc$y)/100

plot(f.both, ylim=c(-0.2,0.2), type="n",
     xlab="Store Locations & Cumulative Market Areas (bw=1)", ylab="Density & log(RR)",
     main="Market Areas by Store Types &\nlog(Relative Risks)")
lines(f.conv, lwd=4, col="red")
lines(f.groc, lwd=4, col="green")
lines(f.both, lwd=4, col="blue")
lines(f.both$x, f.risk, col="purple", lwd=4, lty=2)
abline(h=0, lwd=4)

points(conv, rep(0, length(conv)), col="red", cex=4, pch=20)
points(groc, rep(0, length(groc)), col="green", cex=4, pch=20)

legend("bottomright",
       legend = c("Grocery Stores",
                  "Convenience Stores",
                  "Joint Both Stores",
                  "log(Relative Risk)"),
       col=c("green","red","blue","purple"), lwd=rep(4,4),
       bty="o", title="Store Types & log(RR)")

## ----densityKernelFct, echo=TRUE----------------------------------------------
densityKernelFct <- function(kernel, tess){
  ##
  ## Split a kernel density image into individual tessellations (also type image)
  ## and calculated for each tessellation descriptive statistics of the kernel
  ## density pixel values
  ##
  tract <- split(kernel, tess) # see Baddeley p 122
  #sapply(tract, integral.im)  # kernel sum in tract
  n <- length(tract)           # number of tracts
  ## Initialized data-frame of results
  df <- data.frame(SeqId=1:n, nOfCells=NA, Sum=NA, Density=NA, varDens=NA)
  ## Cycle over tracts. Data values are in $v and may include NAs.
  ## $v is a matrix.
  for (i in 1:n){
    Sum <- sum(tract[[i]]$v, na.rm = T)
    Cells <- length(na.omit(as.vector(tract[[i]]$v)))
    Density <- mean(as.vector(tract[[i]]$v), na.rm = T)
    varDens <- var(as.vector(tract[[i]]$v), na.rm = T)
    df[i,2:5] <- c(Cells, Sum, Density, varDens)
  }
  return(df)
} ## end::densityKernelFct

## ----Reproject, echo=TRUE, fig.height=6.5-------------------------------------
## Check tracts unique sequence Id
plot(tractShp, axes=T, main="Check SeqId")
text(coordinates(tractShp), as.character(tractShp$SeqId), cex=0.5)

## Change the map coordinate system to UTM zone 14 (best fit for Texas)
proj4string(bndShp)                               # Current projection system
projUTM <- CRS("+proj=utm +zone=14  +units=m")    # New coordinate sytem definition
bndUTM <- spTransform(bndShp, projUTM)            # Re-project boundary
ptsUTM <- spTransform(foodStoresShp, projUTM)     # Re-project store locations
tractUTM <- spTransform(tractShp, projUTM)        # Re-project tracts

## Converting UTM shape-files to ppp object and assign marks
ptsDf <- as.data.frame(ptsUTM)
pts <- as.ppp(ptsUTM)
pts$marks <- NULL
pts$marks <- ptsDf$STORETYPE
#pts$marks <- factor(ptsDf$STORETYPE, labels=c("healhy","junk"))

## ----SetResoulution, echo=TRUE------------------------------------------------
## Define window with units and pixel resolution
win <- as.mask(as.owin(bndUTM), eps=200)   # Width and height of pixels = 200 m
unitname(win) <- list("meter","meters")
pts <- pts[win]        # assign owin to pts
summary(pts)

## Convert tracts to tesselation object
tractDf <- as.data.frame(tractUTM)          # Save the attribute informaton
tracts <- slot(tractUTM, "polygons")
tracts <- lapply(tracts, function(x) { SpatialPolygons(list(x)) })
tracts <- lapply(tracts, as.owin)
tractsTess <- tess(tiles=tracts)

## Converting tesselation to grid image with 200x200 meters cell size
tractsTessIm <- tess(tiles=tracts, window=win)

## ----Stores, echo=TRUE, fig.height=6.5----------------------------------------
##############################################################################
## Sample code evaluating the food store market areas by food sales weighted
## kernel density estimates
##############################################################################

## subset cases and controls
groc <- subset(pts, marks=="Healthy Food")
conv <- subset(pts, marks=="Junk Food")
grocFoodSales <- ptsDf[ptsDf$STORETYPE=="Healthy Food","FOODSALES"]
convFoodSales <- ptsDf[ptsDf$STORETYPE=="Junk Food","FOODSALES"]

## ----GroceryStores, echo=TRUE, fig.height=6.5---------------------------------
## Evaluate grocery food stores by census tract (b=4000)
goodFoodIm <- density(groc, weights=grocFoodSales, sigma=4000)
summary(goodFoodIm)
plot(goodFoodIm, main="Weighted Nutritious Food Store Kernel Density\nbw = 4000 m")
plot(tractsTess, add=T)
plot(groc, cex=0.5, pch=16, col="green", add=T)
box(); axis(1); axis(2)
grocFoodStats <- densityKernelFct(goodFoodIm, tractsTessIm)
summary(grocFoodStats)

tractShp$good4000D <- grocFoodStats$Density
tractShp$good4000V <- grocFoodStats$varDens
mapColorRamp(tractShp$good4000D, tractShp, breaks=8, legend.cex = 0.6,
             map.title="Integrated Nutritious Food Store WKD\nbw=4000",
             legend.title= "Grocery Food\nStores", add.to.map=F)
plot(lakesShp, col="skyblue", border="skyblue",add=T)
plot(hwyShp, col="cornsilk3", lwd=3, add=T)
plot(foodDesertShp, border="magenta",lwd=2, add=T)

## ----ConvenienceStores, echo=TRUE, fig.height=6.5-----------------------------
## Evaluate convenience food stores by census tract (bw=2000)
badFoodIm <- density(conv, weights=convFoodSales, sigma=2000)
plot(badFoodIm, main="Weighted Processed Food Store Kernel Density\nbw = 2000")
plot(tractsTess, add=T)
plot(conv, cex=0.5, pch=16, col="bisque", add=T)
box(); axis(1); axis(2)
badFoodStats <- densityKernelFct(badFoodIm, tractsTessIm)
summary(badFoodStats)

tractShp$bad2000D <- badFoodStats$Density
tractShp$bad2000V <- badFoodStats$varDens
mapColorRamp(tractShp$bad2000D, tractShp, breaks=8, legend.cex = 0.6,
             map.title="Integrated Processed Food Store WKD\nbw=2000",
             legend.title= "Convenience Food\nStores", add.to.map=F)
plot(lakesShp, col="skyblue", border="skyblue",add=T)
plot(hwyShp, col="cornsilk3", lwd=3, add=T)
plot(foodDesertShp, border="magenta",lwd=2, add=T)

## ----AllStores, echo=TRUE, fig.height=6.5-------------------------------------
## Evaluate both store types together (bw=3000)
bothStoresIm <- density(unmark(pts), weights=ptsDf$FOODSALES, sigma=3000)
plot(bothStoresIm, main="Weighted All Stores Kernel Density\nbw = 3000 m")
plot(tractsTess, col="grey", add=T)
plot(pts, cex=0.5, pch=16, col="green", add=T)
box(); axis(1); axis(2)
allFoodStats <- densityKernelFct(bothStoresIm, tractsTessIm)
summary(allFoodStats)

tractShp$nCells <- allFoodStats$nOfCells
tractShp$all3000D <- allFoodStats$Density
tractShp$all3000V <- allFoodStats$varDens
mapColorRamp(tractShp$all3000D, tractShp, breaks=8, legend.cex = 0.6,
             map.title="Integrated All Food Stores WKD\nbw=3000",
             legend.title= "All Food\nStores", add.to.map=F)
plot(lakesShp, col="skyblue", border="skyblue",add=T)
plot(hwyShp, col="cornsilk3", lwd=3, add=T)
plot(foodDesertShp, border="magenta",lwd=2, add=T)

## ----RiskMap, echo=TRUE, fig.height=6.5---------------------------------------
##
## Relative Risk evaluation with smacpod::logrr
##
riskMap <- logrr(pts, case=2, sigma=2000, sigmacon=4000,
                 weights=ptsDf$FOODSALES)
bound <- max(abs(range(riskMap$v, na.rm=TRUE)))
plot(riskMap, main="Relative Weighted Risk Surface\nlog(Convenience/Grocery)",
     breaks=seq(-bound, bound, length.out=16+1), col=diverge_hsv(16))
plot(tractsTess, col="grey", add=T)
box(); axis(1); axis(2)

## Evaluate risk by census tract
riskStats <- densityKernelFct(riskMap, tractsTessIm)
summary(riskStats)

tractShp$LRRmedD <- riskStats$Density
tractShp$LRRmedV <- riskStats$varDens
mapBiPolar(tractShp$LRRmedD, tractShp, neg.breaks=4, pos.breaks=6, legend.cex = 0.6,
           map.title="Integrated Relative Weighted Risk Surface with High BW",       
           legend.title= "Log Relative Risk", add.to.map=F)
plot(lakesShp, col="turquoise", border="turquoise",add=T)
plot(hwyShp, col="cornsilk3", lwd=3, add=T)
plot(foodDesertShp, border="magenta",lwd=2, add=T)

## ----RiskMapSignificance, echo=TRUE, fig.height=6.5---------------------------
## Identify areas exessive and lower risk
## decrease error prob: alpa=1-level
riskMapSig <- logrr(pts, case=2, level=0.95, alternative="two.sided",
                    sigma=2000, sigmacon=4000, nsim=99, weights=ptsDf$FOODSALES)
bound <- 5
plot(riskMapSig, 
     main="Significant Relative Risk Surface at 5%\nWhite denotes insignificance",
     breaks=seq(-bound, bound, length.out=16+1), col=diverge_hsv(16))
plot(tractsTess, add=T)
box(); axis(1); axis(2)

## Global Test
logrr.test(riskMapSig)

## ----MapAirports, echo=TRUE, fig.height=6.5-----------------------------------
## Identify two tracts without night population and map them
tractShp@data[tractShp$NIGHTPOP <= 0, "SeqId" ]
tractShp$airport <- 0
tractShp$airport[c(154,248)] <- 1
tractShp$airport <- factor(tractShp$airport, labels=c("no","airport"))

mapColorQual(tractShp$airport, tractShp, legend.cex=0.7, legend.title = "Airport\nTracts",
             map.title = "Airport Census Tracts\nto be Excluded from the Analysis")
plot(lakesShp, col="blue", border="blue",add=T)
plot(hwyShp, col="red", lwd=3, add=T)

## Remove both records from tractShp
tractShp <- tractShp[-c(154,248),]

## ----Regression, echo=TRUE, fig.height=7.2------------------------------------
##
## Save augmented (with Density Estimates) Census Tract shape file
##
# rgdal::writeOGR(tractShp, dsn=getwd(), layer="DallasFoodTracts",
#                 overwrite_layer=T, driver="ESRI Shapefile")
#
# tractShp <- rgdal::readOGR(dsn=".", layer="DallasFoodTracts",
#                            integer64="warn.loss")

## Preliminary regression model
car::scatterplotMatrix(~all3000D+log(POPDEN)+log(PCTDAYPOP)+log(PCTBLACK+1)+PCTUNIVDEG+
                         PCTBADENG +log(PCTHUVAC+1)+log(PCTNOVEH+1)+PCTNOHINS, data=tractShp,
                       main="Potential Variables Impacting the Store-Density",
                       pch=20, smooth=list(span = 0.35,lty.smooth=1, col.smooth="red", col.var="red"),
                       regLine=list(col="green"))
all3000.lm1 <- lm(all3000D~log(POPDEN)+log(PCTDAYPOP)+log(PCTBLACK+1)+PCTUNIVDEG+
                           PCTBADENG +log(PCTHUVAC+1)+log(PCTNOVEH+1)+PCTNOHINS, 
                           data=tractShp)
summary(all3000.lm1)
car::vif(all3000.lm1)
car::residualPlots(all3000.lm1)
all3000.lm2 <- lm(all3000D~log(POPDEN)+log(PCTDAYPOP)+I(log(PCTDAYPOP)^2)+log(PCTBLACK+1)+
                           PCTUNIVDEG+PCTBADENG +log(PCTHUVAC+1)+log(PCTNOVEH+1)+PCTNOHINS, 
                           data=tractShp)
summary(all3000.lm2)
plot(allEffects(all3000.lm2))

car::influenceIndexPlot(all3000.lm2, id=list(n=3))

