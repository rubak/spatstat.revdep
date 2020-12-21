## ----setup, include=FALSE-----------------------------------------------------
knitr::opts_chunk$set(echo = TRUE, fig.pos = 'center', fig.align = 'center')
#, fig.path = "images/"

## ----load_package, message=FALSE, warning=FALSE-------------------------------
library(ipdw)

## ----load_data, message = FALSE, results='hide'-------------------------------
library(rgdal)
pols <- readOGR(system.file("extdata/kattegat_coast.gpkg", package = "ipdw"))
pnts <- readOGR(system.file("extdata/kattegat_pnts.gpkg", package = "ipdw"))

## ----create_costraster--------------------------------------------------------
costras <- costrasterGen(pnts, pols, extent = "pnts", 
												 projstr = projection(pols))
# insert contiguous barrier
costras[160:170,1:80] <- 10000

## ----avg_nearest_neigh, message = FALSE---------------------------------------
# find average nearest neighbor
library(spatstat.core)

W              <- owin(range(coordinates(pnts)[,1]), range(coordinates(pnts)[,2]))
kat.pp         <- ppp(coordinates(pnts)[,1], coordinates(pnts)[,2], window = W)
mean.neighdist <- mean(nndist(kat.pp))

# grid building
gridsize       <- mean.neighdist * 2
grainscale.fac <- gridsize / res(costras)[1]
gridras        <- aggregate(costras, fact = grainscale.fac)
gridpol        <- rasterToPolygons(gridras)
gridpol$value  <- row.names(gridpol)

# spatial join
fulldataset.over    <- over(pnts, gridpol)
fulldataset.over    <- cbind(data.frame(fulldataset.over), 
	setNames(data.frame(pnts), 
	c("id", "salinity", "x.utm", "y.utm", "optional")))

# grid selection
set.seed(2)
gridlev <- unique(fulldataset.over$value)
for(i in seq_along(gridlev)){
  activesub <- subset(fulldataset.over, fulldataset.over$value == gridlev[i])
  selectnum <- gdata::resample(seq_len(nrow(activesub)), 1)
  if(i == 1){
    training <- activesub[selectnum,]
  }
  else{
    training <- rbind(training, activesub[selectnum,])
  }
}

## ----split_training_validation------------------------------------------------
validate             <- fulldataset.over[!(row.names(fulldataset.over) %in% 
															 	row.names(training)),]
xy                   <- cbind(training$x.utm, training$y.utm)
training             <- SpatialPointsDataFrame(xy, training)
xy                   <- cbind(validate$x.utm, validate$y.utm)
validate             <- SpatialPointsDataFrame(xy, validate)
projection(training) <- projection(pnts)
projection(validate) <- projection(pnts)

## ----plot_cost_raster, fig.cap = "<strong>Figure 1: Cost raster representing the high cost of travel through land areas. Training and validation points are shown in black and red respectively.</strong>", fig.height = 6, fig.width = 5----

plot(costras)
points(training)
points(validate, col = "red")

## ----interpolate, cache = FALSE, message = FALSE, results = 'hide'------------
paramlist <- c("salinity")
final.ipdw <- ipdw(training, costras, range = mean.neighdist * 10, paramlist,
									 overlapped = TRUE)

## ----plot_interpolation, fig.cap = "<strong>Figure 2: Interpolated salinity surface by IPDW.</strong>", fig.height = 6, fig.width = 5----
plot(final.ipdw, main = "Kattegat salinity (ppt)")

## ----create_idw, eval=FALSE---------------------------------------------------
#  idw.grid <- rasterToPoints(costras, fun = function(x){ x < 10000}, spatial = TRUE)
#  gridded(idw.grid) <- TRUE
#  kat.idw <- gstat::idw(salinity~1, training, idw.grid, maxdist = mean.neighdist*10,
#  											debug.level = 0)
#  final.idw <- raster(kat.idw)

## ----plot_ipdw_vs_idw, fig.cap = "<strong>Figure 3: Comparison between IPDW and IDW outputs. Note the overestimation of salinity on the upstream (south) side of the contiguous barrier.</strong>", fig.width = 6, fig.height = 4, eval=FALSE----
#  par(mfrow = c(1, 3), mar = c(5.1, 4.1, 4.1, 5.1))
#  plot(final.ipdw, main = "IPDW")
#  plot(final.idw, main = "IDW")
#  plot(final.idw-final.ipdw,  main = "IDW versus IPDW")

## ----plot_ipdw_vs_idw_img, echo=FALSE-----------------------------------------
knitr::include_graphics("images/plot_ipdw_vs_idw-1.png")

## ----generate_validation, eval=FALSE------------------------------------------
#  measured.spdf              <- data.frame(validate$salinity)
#  coordinates(measured.spdf) <- coordinates(validate)
#  
#  valid.ipdw <- errorGen(final.ipdw, measured.spdf, measured.spdf@data)
#  valid.idw  <- errorGen(final.idw, measured.spdf, measured.spdf@data)

## ----plot_validation, fig.cap = "<strong>Figure 4: Comparison between IPDW and IDW interpolation error.  A one-to-one line and best-fit line are shown in black and red respectively.</strong>", fig.width = 8, fig.height = 5, eval=FALSE----
#  par(mfrow = c(1, 2))
#  valid.ipdw <- errorGen(final.ipdw, measured.spdf, measured.spdf@data,
#  											 plot = TRUE, title = "IPDW")
#  valid.idw <- errorGen(final.idw, measured.spdf, measured.spdf@data,
#  											plot = TRUE, title = "IDW")

## ----plot_validation_img, echo=FALSE------------------------------------------
knitr::include_graphics("images/plot_validation-1.png")

