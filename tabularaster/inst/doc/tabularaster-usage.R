## ---- eval=FALSE--------------------------------------------------------------
#  cellnumbers(r, q)

## -----------------------------------------------------------------------------
library(raster)
library(tabularaster)
(r <- raster(volcano))
(cell <- cellnumbers(r, cbind(0.5, 0.5)))

## -----------------------------------------------------------------------------
xyFromCell(r, cell$cell_)

raster::extract(r, cell$cell_)

## -----------------------------------------------------------------------------
library(dplyr)
as_tibble(r)
b <- brick(r, r*2)
as_tibble(b)
as_tibble(b, cell = FALSE) %>% arrange(desc(dimindex)) ## leave out the cell index

## -----------------------------------------------------------------------------
btime <- raster::setZ(b, Sys.time() + c(1, 10))
as_tibble(btime) %>% group_by(dimindex) %>% summarize(n = n())

as_tibble(btime, split_date = TRUE)


## -----------------------------------------------------------------------------
library(tabularaster)
## https://gis.stackexchange.com/questions/102870/step-by-step-how-do-i-extract-raster-values-from-polygon-overlay-with-q-gis-or

library(raster)

# Create integer class raster
r <- raster(ncol=36, nrow=18)
r[] <- round(runif(ncell(r),1,10),digits=0)

# Create two polygons
cds1 <- rbind(c(-180,-20), c(-160,5), c(-60, 0), c(-160,-60), c(-180,-20))
cds2 <- rbind(c(80,0), c(100,60), c(120,0), c(120,-55), c(80,0))
polys <- SpatialPolygonsDataFrame(SpatialPolygons(list(Polygons(list(Polygon(cds1)), 1), 
                                                       Polygons(list(Polygon(cds2)), 2))),data.frame(ID=c(1,2)))

## do extraction in abstract terms
(cn <- cellnumbers(r, polys))

library(dplyr)
## now perform extraction for real
## and pipe into grouping by polygon (object_) and value, and
## calculate class percentage from class counts per polygon
cn %>% 
  mutate(v = raster::extract(r, cell_)) %>% 
  group_by(object_, v) %>% 
  summarize(count = n()) %>% 
  mutate(v.pct = count / sum(count)) 

## here is the traditional code used in the stackoverflow example
# Extract raster values to polygons                             
#( v <- extract(r, polys) )
# Get class counts for each polygon
#v.counts <- lapply(v,table)
# Calculate class percentages for each polygon
#( v.pct <- lapply(v.counts, FUN=function(x){ x / sum(x) } ) )



## -----------------------------------------------------------------------------
library(tabularaster)
data("ghrsst")  ## a RasterLayer
data("sst_regions") ## a polygon layer, contiguous with ghrsst

gcells <- cellnumbers(ghrsst, sst_regions) %>% mutate(object_ = as.integer(object_))

result <- gcells %>% mutate(sst = raster::extract(ghrsst, cell_)) %>% 
  group_by(object_) %>% 
  summarize_at(vars(sst), funs(mean(., na.rm = TRUE), sd(., na.rm = TRUE), length))



## -----------------------------------------------------------------------------
library(tabularaster)
library(raster)
library(dplyr)
data("rastercano")
data("polycano")
cells <- cellnumbers(rastercano, polycano[4:5, ])


cellnumbers(rastercano, as(polycano[4:5, ], "SpatialLinesDataFrame"))
cellnumbers(rastercano, as(as(polycano[4:5, ], "SpatialLinesDataFrame"), "SpatialPointsDataFrame"))

