## ---- include = FALSE---------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)
# load things here for a clean visual experience
library(sf)
library(sp)
library(raster)
library(geometr)
library(tibble)
library(dplyr)
library(readr)
library(magrittr)
getExtent(x = st_read(system.file("shape/nc.shp", package="sf")))

## ---- eval=FALSE--------------------------------------------------------------
#  library(sf)
#  library(sp)
#  library(raster)
#  library(geometr)
#  library(tibble)
#  library(dplyr)
#  library(readr)
#  library(magrittr)

## -----------------------------------------------------------------------------
nc_sf <- st_read(system.file("shape/nc.shp", package="sf"))
st_bbox(nc_sf)

nc_sp <- as_Spatial(nc_sf)
bbox(nc_sp)

ras <- raster(system.file("external/test.grd", package="raster"))
extent(ras)

## -----------------------------------------------------------------------------
myInput <- nc_sf
getExtent(x = myInput)

myInput <- nc_sp
getExtent(x = myInput)

myInput <- ras
getExtent(x = myInput)

## -----------------------------------------------------------------------------
gtGeoms$grid$categorical@point
getPoints(x = gtGeoms$grid$categorical)

## ---- fig.height=3.5----------------------------------------------------------
object.size(gtRasters$categorical)
object.size(gtGeoms$grid$categorical)

visualise(gtRasters$categorical, gtGeoms$grid$categorical)

## -----------------------------------------------------------------------------
gtGeoms$grid$categorical@feature
getFeatures(x = gtGeoms$grid$categorical)

## -----------------------------------------------------------------------------
getGroups(x = gtGeoms$grid$categorical)

## -----------------------------------------------------------------------------
# transform from class sf
(nc_geom <- gc_geom(input = nc_sf))

# make "by hand"
coords <- tibble(x = c(40, 70, 70, 50),
                 y = c(40, 40, 60, 70))
(somePoints <- gs_point(anchor = coords))

## -----------------------------------------------------------------------------
(aPoly <- gs_polygon(anchor = coords))

## -----------------------------------------------------------------------------
pond <- tibble(x = c(30, 30, 80, 80, 30),
               y = c(30, 80, 80, 30, 30),
               fid = 1)
island <- tibble(x = c(60, 65, 65, 60, 60),
                y = c(45, 45, 50, 50, 45),
                fid = 1)
temp <- bind_rows(pond, island, getPoints(aPoly))

(perforatedPoly <- gs_polygon(anchor = temp))
visualise(perforatedPoly, fillcol = fid)

## -----------------------------------------------------------------------------
new_geom <- gc_geom(input = nc_sf, group = TRUE)

getFeatures(x = new_geom)

getGroups(x = new_geom)

## -----------------------------------------------------------------------------
clg_dat <- read_delim(file = "http://www.pommerening.org/wiki/images/d/dc/Clg6.txt", 
                        delim = "\t", col_types = "iiidddd")
clg_spec <- read_delim(file = "http://www.pommerening.org/wiki/images/d/df/Clg6.species",
                       delim = "\t", col_types = "icccc") %>% 
  mutate_if(is.character, trimws)

# make table of locations, tree properties and species metadata
locations <- tibble(x = clg_dat$x,
                    y = clg_dat$y,
                    fid = clg_dat$Number)
prop <- tibble(fid = clg_dat$Number,
               gid = clg_dat$Species,
               dbh = clg_dat$dbh,
               height = clg_dat$ht)
species <- tibble(gid = unique(clg_dat$Species)) %>% 
  left_join(clg_spec, by = c("gid" = "NC")) %>% 
  select(gid, name = CN, scientific = BN) %>% 
  mutate(deciduous = c(T, T, F, F, F, T, T, F, T, T, F)) %>% 
  arrange(gid)

# set these as attribute tables and select only trees with a height value
trees <- locations %>%
  gs_point(anchor = .) %>%
  setFeatures(table = prop) %>%
  setGroups(table = species)

# make the title
validpoints <- getPoints(x = trees)
title <- paste0("Clocaenog - tree heights (", dim(validpoints)[1], ")")

# and visualise
visualise(!!title := trees, linecol = height)
visualise(!!title := getGroups(trees, deciduous), linecol = height)

## -----------------------------------------------------------------------------
visualise(aPoly)

reflPoly <- gt_reflect(geom = aPoly, angle = 45)

visualise(reflPoly, linecol = "deeppink", new = FALSE, trace = TRUE)

## -----------------------------------------------------------------------------
# make points and discard duplicates ...
(somePoints <- gs_point(anchor = aPoly) %>%
  getPoints(!(duplicated(x) & duplicated(y))))

# ... and recreate a polygon from that
(aPoly <- gs_polygon(anchor = somePoints))

## ---- fig.height=3.5----------------------------------------------------------
visualise('I made a polygon!' = aPoly, gtRasters$continuous)

## ---- fig.height=3.5----------------------------------------------------------
relPoly <- gs_rectangle(anchor = aPoly) %>% 
  setWindow(to = tibble(x = c(0, 100), y = c(0, 100))) %>%
  gt_scale(to = "relative")

visualise(gtRasters$categorical, gtRasters$continuous)
visualise(relPoly, linecol = "deeppink", new = FALSE)

## -----------------------------------------------------------------------------
zoom <- setWindow(x = relPoly, to = getExtent(x = gtRasters$categorical)) %>% 
  gt_scale(to = "absolute")
visualise('zoomed object' = gtRasters$categorical, window = getExtent(x = zoom))

## -----------------------------------------------------------------------------
anImage <- brick(system.file("external/rlogo.grd", package="raster"))
visualise('R logo' = anImage, image = TRUE)

## -----------------------------------------------------------------------------
gtTheme

## -----------------------------------------------------------------------------
treeTheme <- setTheme(scale = list(param = "pointsize", to = "dbh"), 
                      vector = list(pointsize = c(0.05, 0.77)))
visualise(`Clocaenog - tree diameters` = trees, theme = treeTheme)

## -----------------------------------------------------------------------------
treeTheme <- setTheme(vector = list(pointsize = c(0.05, 0.77)))
visualise(`Clocaenog - tree diameters and heights` = trees, theme = treeTheme,
          pointsize = dbh, linecol = height)

## -----------------------------------------------------------------------------
treeTheme <- setTheme(vector = list(linecol = c("#00204DFF", "#73D216FF", "#FFEA46FF", "#CC0000F0")))
visualise(`Clocaenog - tree heights` = trees, theme = treeTheme, linecol = height)

