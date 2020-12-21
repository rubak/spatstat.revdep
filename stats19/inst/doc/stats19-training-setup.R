## ---- include = FALSE---------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  out.width = "50%"
)

## ---- eval=FALSE--------------------------------------------------------------
#  source("https://git.io/JeaZH")

## ---- eval=FALSE--------------------------------------------------------------
#  install.packages("remotes")
#  pkgs = c(
#    "pct",         # package for getting travel data in the UK
#    "sf",          # spatial data package
#    "stats19",     # downloads and formats open stats19 crash data
#    "stplanr",     # for working with origin-destination and route data
#    "tidyverse",   # a package for user friendly data science
#    "tmap"         # for making maps
#  )
#  remotes::install_cran(pkgs)
#  # remotes::install_github("ITSLeeds/pct")

## ----message=FALSE, eval=FALSE------------------------------------------------
#  library(stats19)
#  library(tidyverse)
#  library(tmap) # installed alongside mapview
#  crashes = get_stats19(year = 2017, type = "ac")
#  crashes_iow = crashes %>%
#    filter(local_authority_district == "Isle of Wight") %>%
#    format_sf()
#  
#  # basic plot
#  plot(crashes_iow)

