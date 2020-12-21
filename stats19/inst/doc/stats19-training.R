## ---- message=FALSE, warning=FALSE, eval=FALSE--------------------------------
#  library(pct)      # access travel data from DfT-funded PCT project
#  library(sf)       # spatial vector data classes
#  library(stats19)  # get stats19 data
#  library(stplanr)  # transport planning tools
#  library(tidyverse)# packages for 'data science'
#  library(tmap)     # interactive maps

## ---- include = FALSE---------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  out.width = "50%"
)

## ----pkgs, warning=FALSE, echo=FALSE------------------------------------------
pkgs = c(
  "sf",          # spatial data package
  "stats19",     # downloads and formats open stats19 crash data
  "dplyr",       # a package data manipulation, part of the tidyverse
  "tmap"         # for making maps
)

## ----cite, echo=FALSE---------------------------------------------------------
knitr::write_bib(x = pkgs, "packages.bib")

## ---- eval=FALSE, echo=FALSE--------------------------------------------------
#  remotes::install_cran(pkgs)
#  # remotes::install_github("ITSLeeds/pct")

## ----rstudioui, echo=FALSE, out.width="70%"-----------------------------------
knitr::include_graphics("https://raw.githubusercontent.com/ITSLeeds/TDS/master/courses/2day/images/rstudio-ui.png")

## ----edit, eval=FALSE---------------------------------------------------------
#  file.edit("stats19-lesson-1.R")

## ---- eval=FALSE--------------------------------------------------------------
#  x = 1:5
#  y = c(0, 1, 3, 9, 18)
#  plot(x, y)

## -----------------------------------------------------------------------------
vehicle_type = c("car", "bus", "tank")
casualty_type = c("pedestrian", "cyclist", "cat")
casualty_age = seq(from = 20, to = 60, by = 20)
set.seed(1)
dark = sample(x = c(TRUE, FALSE), size = 3, replace = TRUE)
small_matrix = matrix(1:24, nrow = 12)
crashes = data.frame(vehicle_type, casualty_type, casualty_age, dark)

## ----summary------------------------------------------------------------------
summary(casualty_age)

## ----summary-answers, echo=FALSE, eval=FALSE----------------------------------
#  summary(vehicle_type)
#  class(vehicle_type)
#  typeof(vehicle_type)
#  dim(vehicle_type)
#  length(vehicle_type)

## ----tibble1, echo=FALSE, eval=FALSE------------------------------------------
#  tibble::tibble(
#    vehicle_type,
#    casualty_type,
#    casualty_age,
#    dark
#  )

## ----autocomp, echo=FALSE-----------------------------------------------------
knitr::include_graphics("https://raw.githubusercontent.com/ITSLeeds/TDS/master/courses/2day/images/autocomplete.jpg")

## ----help, echo=FALSE---------------------------------------------------------
knitr::include_graphics("https://raw.githubusercontent.com/ITSLeeds/TDS/master/courses/2day/images/fucntionhelp.jpg")

## -----------------------------------------------------------------------------
# Create vector objects (a whole line comment)
x = 1:5 # a seqence of consecutive integers (inline comment)
y = c(0, 1, 3, 9, 18.1) 

## ----debug, echo=FALSE, out.width="60%"---------------------------------------
knitr::include_graphics("https://raw.githubusercontent.com/ropensci/stats19/master/inst/rstudio-autocomplete.png")

## -----------------------------------------------------------------------------
saveRDS(crashes, "crashes.Rds")

## -----------------------------------------------------------------------------
crashes2 = readRDS("crashes.Rds")
identical(crashes, crashes2)

## ----readr-write, eval=FALSE--------------------------------------------------
#  readr::write_csv(crashes, "crashes.csv")
#  crashes3 = readr::read_csv("crashes.csv")
#  identical(crashes3, crashes)

## ---- eval=FALSE--------------------------------------------------------------
#  casualty_age[2:3] # second and third casualty_age
#  crashes[c(1, 2), ] # first and second row of crashes
#  crashes$vehicle_type # returns just one column
#  crashes[, c("casualty_type", "casualty_age")] # first and third columns

## ---- eval=FALSE, echo=FALSE--------------------------------------------------
#  crashes[, c(1, 3)] # first and third column of crashes by positional numbers
#  crashes[c(2), c(3)]
#  crashes[c(2), c(2, 3)]
#  class(crashes[, c(1, 3)])
#  class(crashes[c(2), c(3)])

## ---- eval=FALSE--------------------------------------------------------------
#  x[c(TRUE, FALSE, TRUE, FALSE, TRUE)] # 1st, 3rd, and 5th element in x
#  x[x == 5] # only when x == 5 (notice the use of double equals)
#  x[x < 3] # less than 3
#  x[x < 3] = 0 # assign specific elements
#  casualty_age[casualty_age %% 6 == 0] # just the ages that are a multiple of 6
#  crashes[crashes$dark == FALSE, ]

## ---- eval=FALSE, echo=FALSE--------------------------------------------------
#  casualty_age[casualty_age < 50] # the  casualty_age less than 50
#  crashes[crashes$vehicle_type == "tank", ] # rows where the name is tank
#  crashes$casualty_age[crashes$vehicle_type == "tank"] = 61

## ---- eval=FALSE--------------------------------------------------------------
#  z = c(4, 5, NA, 7)

## ---- eval=FALSE--------------------------------------------------------------
#  sum(z) # result is NA

## ---- eval=FALSE--------------------------------------------------------------
#  sum(z, na.rm = TRUE) # result is equal to 4 + 5 + 7

## ---- eval=FALSE--------------------------------------------------------------
#  is.na(z)
#  z_nona = z[!is.na(z)] # note the use of the not operator !
#  sum(z)

## ---- echo=FALSE, eval=FALSE--------------------------------------------------
#  crashes$vehicle_type = as.character(crashes$vehicle_type)
#  as.matrix(crashes)

## -----------------------------------------------------------------------------
z = c(1, 2, -1, 1, 3)
l = c(NA, "a", "b", "c") # labels in ascending order
z_factor = factor(z, labels = l)
z_charcter = as.character(z_factor)
z_charcter

## ----smile, out.width="30%", fig.align="center"-------------------------------
eyes = c(2.3, 4, 3.7, 4)
eyes = matrix(eyes, ncol = 2, byrow = T)
mouth = c(2, 2, 2.5, 1.3, 3, 1, 3.5, 1.3, 4, 2)
mouth = matrix(mouth, ncol = 2, byrow = T)
plot(eyes, type = "p", main = "RRR!", cex = 2, xlim = c(1, 5), ylim = c(0, 5))
lines(mouth, type = "l", col = "red")

## ---- eval=FALSE--------------------------------------------------------------
#  install.packages("sf")
#  # remotes::install_github("r-spatial/sf")

## -----------------------------------------------------------------------------
library(sf)

## ----tibble2, eval=FALSE------------------------------------------------------
#  crashes_tibble = tibble::tibble(
#    vehicle_type,
#    casualty_type,
#    casualty_age,
#    dark
#  )

## ---- message=FALSE, out.width="40%", eval=FALSE------------------------------
#  library(ggplot2)
#  ggplot(crashes) + geom_point(aes(x = casualty_type, y = casualty_age))

## ----gg-extend, echo=FALSE, message=FALSE, eval=FALSE-------------------------
#  library(ggplot2)
#  # install.packages("ggthemes")
#  g1 = ggplot(crashes) + geom_point(aes(x = casualty_type, y = casualty_age))
#  g2 = ggplot(crashes) + geom_point(aes(x = casualty_type, y = casualty_age)) +
#    ggthemes::theme_economist()
#  g3 = cowplot::plot_grid(g1, g2)
#  ggsave(filename = "inst/ggtheme-plot.png", width = 8, height = 2, dpi = 80)

## ----gg2, echo=FALSE, out.width="80%", fig.align="center"---------------------
library(ggplot2)
knitr::include_graphics("https://raw.githubusercontent.com/ropensci/stats19/b4c40ad4c134853007493a9eac116b00acd4ec5a/inst/ggtheme-plot.png")

## -----------------------------------------------------------------------------
library(dplyr)
class(crashes)       
crashes %>% class()

## ---- eval=FALSE--------------------------------------------------------------
#  crashes %>%
#    filter(casualty_age > 50) # filter rows
#  crashes %>%
#    select(casualty_type) # select just one column
#  crashes %>%
#    group_by(dark) %>%
#    summarise(mean_age = mean(casualty_age))

## ----dplyr, eval=FALSE, echo=FALSE--------------------------------------------
#  # answers
#  crashes %>%
#    arrange(desc(casualty_age))
#  crashes %>% filter(casualty_age > 21)
#  crashes %>%
#    mutate(birth_year = 2019 - casualty_age) %>%
#    filter(birth_year > 1969)

## ---- message=FALSE-----------------------------------------------------------
library(lubridate)

## -----------------------------------------------------------------------------
today()

## ---- eval=FALSE--------------------------------------------------------------
#  x = today()
#  day(x)
#  month(x)
#  year(x)
#  weekdays(x)

## ---- eval=FALSE--------------------------------------------------------------
#  as.Date("2019-10-17") # works
#  as.Date("2019 10 17") # fails
#  ymd("2019 10 17") # works
#  dmy("17/10/2019") # works

## -----------------------------------------------------------------------------
x = c("2009-01-01", "2009-02-02", "2009-03-03")
x_date = ymd(x)
x_date

## ---- echo=FALSE, eval=FALSE--------------------------------------------------
#  # 1. Extract the day, the year-day, the month and the weekday (as a non-abbreviated character vector) of each element of `x_date`.
#  day(x_date)
#  yday(x_date)
#  month(x_date)
#  weekdays(x_date, abbreviate = FALSE)
#  # 1. Modify the previous example to parse the following character string: `"09/09/1993"` and extract its weekday.
#  weekdays(dmy("09/09/93"))
#  wday(dmy("09/09/93"))

## -----------------------------------------------------------------------------
crashes$casualty_day = x_date

## -----------------------------------------------------------------------------
filter(crashes, day(casualty_day) < 7) # the events that ocurred in the first week of the month
filter(crashes, weekdays(casualty_day) == "Monday") # the events occurred on monday

## -----------------------------------------------------------------------------
x = c("18:23:35", "00:00:01", "12:34:56")
x_hour = hms(x)
x_hour

## -----------------------------------------------------------------------------
hour(x_hour)
minute(x_hour)
second(x_hour)

## -----------------------------------------------------------------------------
x = c("18:23", "00:00", "12:34")
(x_hour = hm(x))

## -----------------------------------------------------------------------------
crashes$casualty_hms = hms(c("18:23:35", "00:00:01", "12:34:56"))
crashes$casualty_hour = hour(crashes$casualty_hms)

## ---- eval=FALSE--------------------------------------------------------------
#  library(stats19)
#  crashes_2017 = stats19::get_stats19(year = 2017, type = "ac")
#  crashes_2017

## ---- eval=FALSE, echo=FALSE--------------------------------------------------
#  # solutions
#  crashes %>% filter(casualty_hour >= 12)
#  crashes %>% filter(casualty_hour > 15 & casualty_hour < 19)
#  
#  crashes_2017 %>%
#    mutate(my_weekdays = weekdays(date)) %>%
#    filter(my_weekdays == "Monday") %>%
#    nrow()
#  crashes_2017 %>%
#    mutate(my_weekdays = weekdays(date)) %>%
#    filter(my_weekdays == "Friday") %>%
#    nrow()
#  
#  crashes_2017 %>%
#    mutate(my_weekdays = weekdays(date)) %>%
#    group_by(my_weekdays) %>%
#    summarize(n = n()) %>%
#    ggplot() +
#    geom_col(aes(x = my_weekdays, y = n))
#  
#  crashes_2017 %>%
#    mutate(my_hours = hour(hm(time))) %>%
#    group_by(my_hours) %>%
#    summarize(n = n()) %>%
#    ggplot() +
#    geom_col(aes(x = my_hours, y = n))
#  
#  crashes_2017 %>%
#    mutate(my_weekdays = weekdays(date), my_hours = hour(hm(time))) %>%
#    group_by(my_weekdays, my_hours) %>%
#    summarise(n = n()) %>%
#    ggplot() +
#    geom_line(aes(x = my_hours, y = n, col = my_weekdays), size = 1.05)
#  # the legend needs some reordering

## ----crashes-sf, fig.height=2, fig.width=3------------------------------------
library(sf) # load the sf package for working with spatial data
crashes_sf = crashes # create copy of crashes dataset
crashes_sf$longitude = c(-1.3, -1.2, -1.1)
crashes_sf$latitude = c(50.7, 50.7, 50.68)
crashes_sf = st_as_sf(crashes_sf, coords = c("longitude", "latitude"), crs = 4326)
# plot(crashes_sf[1:4]) # basic plot
# mapview::mapview(crashes_sf) # for interactive map

## ----crashes-sf-ex, echo=FALSE, out.width="30%", fig.show='hold'--------------
plot(crashes_sf$geometry)
plot(crashes_sf["casualty_age"])
plot(crashes_sf[2:3, "dark"])
# st_distance(crashes_sf)
# Bembridge

# # updload geographic crash data
# write_sf(crashes_sf, "crashes_sf.geojson")
# piggyback::pb_upload("crashes_sf.geojson")

## ---- eval=FALSE--------------------------------------------------------------
#  write_sf(zones, "zones.geojson") # save geojson file
#  write_sf(zones, "zmapinfo", driver = "MapInfo file")
#  read_sf("zmapinfo") # read in mapinfo file

## -----------------------------------------------------------------------------
zones = pct::get_pct_zones("isle-of-wight")[1:9]

## ---- echo=FALSE--------------------------------------------------------------
# class(zones)
# names(zones)
zones[1:2, c(1, 5, 6, 7, 8)]

## ---- message=FALSE-----------------------------------------------------------
zones_containing_crashes = zones[crashes_sf, ]

## ----sp-ex, echo=FALSE, out.width="33%", fig.show='hold', message=FALSE, warning=FALSE----
plot(zones$geometry)
plot(zones_containing_crashes$geometry, col = "red", add = TRUE)
plot(zones$geometry)
plot(zones[crashes_sf[2, ], ], col = "blue", add = TRUE)
plot(zones$geometry)
plot(zones[zones_containing_crashes, ], col = "yellow", add = TRUE)
plot(crashes_sf$geometry, pch = 20, add = TRUE)

## ---- message=FALSE-----------------------------------------------------------
zones_joined = st_join(zones[1], crashes_sf)

## ----joinf, echo=FALSE, out.width="40%", fig.show='hold', message=FALSE-------
plot(zones_joined["casualty_age"])
zjd = st_join(zones[1], crashes_sf["dark"], left = FALSE)
plot(zjd)

## ----crs1---------------------------------------------------------------------
crashes_osgb = st_transform(crashes_sf, 27700)

## ---- eval=TRUE---------------------------------------------------------------
# load example dataset if it doesn't already exist
zones = pct::get_pct_zones("isle-of-wight")
sel = zones$all > 3000  # create a subsetting object
zones_large = zones[sel, ] # subset areas with a popualtion over 100,000
zones_2 = zones[zones$geo_name == "Isle of Wight 002",] # subset based on 'equality' query
zones_first_and_third_column = zones[c(1, 3)]
zones_just_all = zones["all"]

## ---- echo=FALSE, eval=FALSE--------------------------------------------------
#  # 1. Practice subsetting techniques you have learned on the `sf data.frame` object `zones`:
#  #      1. Create an object called `zones_small` which contains only regions with less than 3000 people in the `all` column
#  # in base R
#  zones_small = zones[zones$all < 3000, ]
#  # with dplyr
#  zones_small = zones %>%
#    filter(all < 3000)
#  #      1. Create a selection object called `sel_high_car` which is `TRUE` for regions with above median numbers of people who travel by car and `FALSE` otherwise
#  median_car = median(zones$car_driver)
#  sel_high_car = zones$car_driver > median_car
#  #      1. How many regions have the number '1' in the column 'geo_name'? What percentage of the regions in the Isle of Wight is this?
#  sel_region_name_contains_1 = grepl("1", x = zones$geo_name)
#  sum(sel_region_name_contains_1) / nrow(zones)
#  #      1. Create an object called `zones_foot` which contains only the foot attribute from `zones`
#  # using base R
#  zones_foot = zones["foot"]
#  # dplyr
#  zones_foot = zones %>%
#    select(foot)
#  #      1. Bonus: plot the result to show where walking is a popular mode of travel to work
#  plot(zones_foot)
#  #      1. Bonus: bulding on your answers to previous questions, use `filter()` from the `dplyr` package to subset small regions where high car use is high
#  zones_small_car_high = zones %>%
#    filter(all < 3000, car_driver > median_car)
#  # 1. Bonus: What is the population density of each region (hint: you may need to use the functions `st_area()`, `as.numeric()` and use the 'all' column)?
#  zones$area_km2 = as.numeric(st_area(zones)) /1000000
#  zones$population_density = zones$all / zones$area_km2
#  plot(zones["population_density"])
#  # in dplyr
#  zones_density = zones %>%
#    mutate(area_km2 = as.numeric(st_area(geometry)) / 1000000) %>%
#    mutate(population_density = all / area_km2)
#  plot(zones_density %>% select(population_density))
#  # 1. Bonus: Which zone has the highest percentage who cycle?
#  zones %>%
#    mutate(pcycle = bicycle / all) %>%
#    top_n(n = 1, wt = pcycle)
#  # 1. Bonus: Find the proportion of people who drive to work (`car_driver`) in areas in which more than 500 people walk to work
#  zones %>%
#    group_by(foot > 500) %>%
#    summarise(mean_car = sum(car_driver) / sum(all) )

## -----------------------------------------------------------------------------
library(tmap)
tmap_mode("plot")

## ----plot3, fig.show='hold', out.width="33%", echo=FALSE----------------------
plot(zones[c("all", "bicycle")])
tm_shape(zones) + 
  tm_polygons(c("all", "bicycle"))
tmap_mode("view")
m = tm_shape(zones_joined) + 
  tm_polygons(c("casualty_type")) +
  tm_scale_bar()
m
# knitr::include_graphics("tmap-zones-interactive.png")

## ---- eval=FALSE--------------------------------------------------------------
#  vignette(package = "stats19") # view all vignettes available on stats19
#  vignette("stats19") # view the introductory vignette

## ---- echo=FALSE, results='hide', message=FALSE, eval=FALSE-------------------
#  library(stats19)
#  library(dplyr)
#  library(sf)
#  a = get_stats19(2018, "ac")
#  asf = format_sf(a)
#  a_zones = asf %>%
#    filter(local_authority_district == "Isle of Wight")
#  nrow(a_zones)
#  zones = pct::get_pct_zones(region = "isle-of-wight")
#  zones_osbg = st_transform(zones, 27700)
#  a_zones_sf = a_zones[zones_osbg, ]
#  nrow(a_zones_sf)
#  # mapview::mapview(zones) +
#  #   mapview::mapview(a_zones)
#  class(a$date)
#  class(a$time)
#  a_zones$month = lubridate::month(a_zones$date)
#  a_zones_may = a_zones %>%
#    filter(month == 5)
#  a_agg = aggregate(a_zones_may["speed_limit"], zones_osbg, mean)
#  plot(a_agg)
#  class(a$date)

## -----------------------------------------------------------------------------
u = "https://github.com/ropensci/stats19/releases/download/1.1.0/roads_key.Rds"
roads_wgs = readRDS(url(u))
roads = roads_wgs %>% st_transform(crs = 27700)

## -----------------------------------------------------------------------------
u = "https://github.com/ropensci/stats19/releases/download/1.1.0/car_accidents_2017_iow.Rds"
crashes_iow = readRDS(url(u))

## ---- echo=FALSE, out.width="49%", fig.show='hold', message=FALSE-------------
plot(roads$geometry)
plot(crashes_iow["accident_severity"], add = TRUE)
roads_buffer = st_buffer(roads, 200, endCapStyle = "FLAT")
crashes_outside_roads = crashes_iow[roads_buffer, , op = sf::st_disjoint]
roads_agg = aggregate(crashes_iow[1], by = roads_buffer, FUN = length)
# plot(roads_agg, border = NA, main = "")
names(roads_agg)[1] = "N. Crashes"
tmap_mode("plot")
tm_shape(roads_agg) + tm_fill("N. Crashes") +
  tm_shape(crashes_outside_roads) + tm_dots(col = "blue")

## ----final-plot, echo=FALSE, out.width="100%"---------------------------------
# knitr::include_graphics("final-figure.png")

