## ----setup, include = FALSE---------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  fig.width=12, 
  fig.height=8,
  out.width="100%")
library(foster)
library(ggplot2)
library(raster)
library(knitr)

## ---- echo=FALSE, out.width="80%", fig.cap="FOSTER workflow"------------------
knitr::include_graphics("https://raw.githubusercontent.com/mqueinnec/storage/main/FOSTER/FOSTER_workflow.png")

## -----------------------------------------------------------------------------
library(foster)
library(raster)

## ---- eval = FALSE------------------------------------------------------------
#  #Read single layer raster
#  raster_layer <- raster("path/to/raster_layer.tif")
#  #Read multi-layer raster
#  raster_stack <- stack("path/to/raster_stack.tif")

## ---- eval = FALSE------------------------------------------------------------
#  dem_samples <- rgdal::readOGR(dsn = "path/to/directory",
#                                layer = "layer_name")

## ---- eval=FALSE--------------------------------------------------------------
#  #Example to write output of calcIndices to disk using different options
#  ind <- calcIndices(x, indices = c("NDVI","TCG","TCW"),red = 3,nir=4,filename = "full/path/to/filename.tif")
#  ind <- calcIndices(x, indices = c("NDVI","TCG","TCW"),red = 3,nir=4,filename = "full/path/to/filename", format="GTiff", overwrite=TRUE, byLayer=TRUE)

## ---- eval=FALSE--------------------------------------------------------------
#  rasterOptions(tmpdir) <- "path/to/tempdir"

## -----------------------------------------------------------------------------
elev_p95 <- raster(system.file("extdata/vignette/ALS_metrics/ALS_metrics_p95.tif",package="foster"))
cover <- raster(system.file("extdata/vignette/ALS_metrics/ALS_metrics_cov_mean.tif",package="foster"))
Y_vars <- stack(elev_p95,cover)

#Set layers names
names(Y_vars) <- c("p95","cover")
Y_vars
plot(Y_vars)

## -----------------------------------------------------------------------------
spectral_2006 <- stack(system.file("extdata/vignette/Landsat_BAP/Landsat_BAP_2006.tif",package="foster"))
names(spectral_2006) <- c("blue","green","red","nir","swir1","swir2")

spectral_2006
plot(spectral_2006, col = grey.colors(255))


## -----------------------------------------------------------------------------
dem <- raster(system.file("extdata/vignette/topo/DEM.tif",package="foster"))
dem_slope <- raster(system.file("extdata/vignette/topo/DEM_slope.tif",package="foster"))

plot(dem)
plot(dem_slope)

## -----------------------------------------------------------------------------
mask_forest <- raster(system.file("extdata/vignette/landcover/VLCE_forest_2008.tif",package="foster"))
plot(mask_forest)

## -----------------------------------------------------------------------------
Y_vars_resampled <- matchResolution(x = Y_vars,
                                    ref = spectral_2006,
                                    method='bilinear')
Y_vars_resampled

## -----------------------------------------------------------------------------
filt <- matrix(1,nrow=3,ncol=3)
Y_vars_smooth <- focalMultiBand(Y_vars_resampled,
                                w=filt,
                                fun=mean,
                                pad=TRUE,
                                padValue=NA, 
                                na.rm=TRUE, 
                                keepNA = TRUE)
plot(Y_vars_smooth)

## ---- out.width="100%"--------------------------------------------------------
Y_vars_mask <- matchExtent(Y_vars_smooth,
                           mask_forest,
                           mask=TRUE,
                           maskValue = NA)
plot(Y_vars_mask)

# We do the same with the DEM and slope

dem_mask <- matchExtent(dem, 
                        mask_forest, 
                        mask = TRUE)

dem_slope_mask <- matchExtent(dem_slope, 
                        mask_forest, 
                        mask = TRUE)

## -----------------------------------------------------------------------------
indices_list <- c("NDVI","TCB","TCW","TCG")
# Example for one year
VI_2006 <- calcIndices(spectral_2006,
                   indices = indices_list, 
                   sat="Landsat5TM", 
                   red=3, 
                   nir=4)
plot(VI_2006[["TCG"]])

## -----------------------------------------------------------------------------
# List all filenames in time-series order: 2006 to 2008
spectral_ts_files <- list(system.file("extdata/vignette/Landsat_BAP/Landsat_BAP_2006.tif",package = "foster"), 
                    system.file("extdata/vignette/Landsat_BAP/Landsat_BAP_2007.tif",package = "foster"), 
                    system.file("extdata/vignette/Landsat_BAP/Landsat_BAP_2008.tif",package = "foster"))
  
# Open all Landsat images and place them in a list
spectral_ts <- lapply(spectral_ts_files, stack)
spectral_ts

## -----------------------------------------------------------------------------
VI_ts <- calcIndices(spectral_ts,
                     indices = indices_list, 
                     sat="Landsat5TM", 
                     red=3, 
                     nir=4)

# The output is a list with VI for each year
VI_ts

## -----------------------------------------------------------------------------
plot(VI_ts[[1]])

## -----------------------------------------------------------------------------
VI_ts_smooth <- focalMultiBand(VI_ts,
                             w=filt,
                             fun=mean,
                             na.rm=TRUE,
                             pad = TRUE,
                             keepNA = TRUE)

## -----------------------------------------------------------------------------
plot(VI_ts_smooth[[1]])

## -----------------------------------------------------------------------------
funSummary <- function(x){
  
  c(
    median = median(x,na.rm=TRUE),
    IQR = IQR(x,na.rm=TRUE),
    slope = theilSen(x)
  )
}

## -----------------------------------------------------------------------------
VI_ts_metrics <- temporalMetrics(VI_ts,
                               metrics="funSummary")

## -----------------------------------------------------------------------------
VI_ts_metrics

plot(VI_ts_metrics[["NDVI_median"]])

## -----------------------------------------------------------------------------
Y_vars_edges <- edges(Y_vars_mask,
                      w=3)
plot(Y_vars_edges)

## -----------------------------------------------------------------------------
nSamples = 230
nClasses = 5
mindist = 75

set.seed(1234) #For example reproducibility
sample_strata <- getSample(Y_vars_edges, 
                     layers = c("p95","cover"), 
                     n = nSamples, 
                     strata = nClasses,
                     mindist = mindist, 
                     norm = TRUE,
                     xy = TRUE)

## -----------------------------------------------------------------------------
# Sampled points
sampleLoc <- sample_strata$sample

# Map of strata
strata_map <- sample_strata$clusterMap

# k-NN model
kmeans_model <- sample_strata$model

## -----------------------------------------------------------------------------
plot(strata_map)
plot(sampleLoc,add=TRUE)

## -----------------------------------------------------------------------------
# Predictor variables
X_vars <- stack(VI_ts_metrics, 
                dem_mask, 
                dem_slope_mask)

# Response variables

Y_vars <- Y_vars_mask


## -----------------------------------------------------------------------------

# Extract values at sample
X_vars_sample <- getSampleValues(X_vars, sampleLoc)
Y_vars_sample <- getSampleValues(Y_vars, sampleLoc)

X_vars_sample
Y_vars_sample

## -----------------------------------------------------------------------------
#Create data partition
set.seed(1234) #for example reproducibility
train_idx <- partition(sampleLoc$cluster,
                       type="kfold", 
                       kfold = 5,
                       returnTrain = TRUE)

train_idx

## -----------------------------------------------------------------------------
set.seed(1234) #for example reproducibility
kNN <- trainNN(x = X_vars_sample,
              y=Y_vars_sample,
              inTrain = train_idx, 
              k = 1, 
              method = "randomForest",
              ntree = 200)

## -----------------------------------------------------------------------------
kNN_model <- kNN$model

kNN_preds <- kNN$preds
head(kNN_preds, 10)

## -----------------------------------------------------------------------------
accuracy(obs = kNN_preds$obs, 
         preds = kNN_preds$preds, 
         vars = kNN_preds$variable, 
         folds = kNN_preds$Fold)


## -----------------------------------------------------------------------------
plots_scatter <- scatter(obs = kNN_preds$obs, 
         preds = kNN_preds$preds, 
         vars = kNN_preds$variable)

plots_scatter$cov

plots_scatter$p95

## -----------------------------------------------------------------------------
imp <- varImp(kNN_model,scaled=FALSE,plot=TRUE,plotType="boxplot")
imp$plot
imp <- varImp(kNN$model,scaled=TRUE,plot=TRUE,plotType="grid")
imp$plot

## -----------------------------------------------------------------------------
Y_imputed <- predictTrgs(model=kNN_model,
                         x = X_vars, 
                         nnID = TRUE, 
                         nnDist = TRUE)
Y_imputed

## -----------------------------------------------------------------------------
plot(Y_imputed$p95)
plot(Y_imputed$cover)
plot(Y_imputed$nnID1,col=rainbow(length(unique(Y_imputed$nnID1))))
plot(Y_imputed$nnDist1)

