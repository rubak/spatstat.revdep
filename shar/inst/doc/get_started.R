## ----setup, include = FALSE----------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----load_libraries_hidden, eval = TRUE, echo = FALSE, message = FALSE, results = "hide"----
library(dplyr)
library(shar)
library(spatstat)
library(raster)

## ----classify_habitats---------------------------------------------------
landscape_discrete <- classify_habitats(raster = landscape, classes = 5)

## ----randomize_raster, eval = FALSE--------------------------------------
#  torus_trans <- translate_raster(raster = landscape_discrete)
#  
#  random_walk <- randomize_raster(raster = landscape_discrete, n_random = 39)

## ----randomize_pp, eval = FALSE------------------------------------------
#  gamma_test <- fit_point_process(pattern = species_a, n_random = 39, process = "cluster")
#  
#  reconstruction <- reconstruct_pattern_homo(pattern = species_a, n_random = 39) # takes some time

## ----results-------------------------------------------------------------
results_habitat_association(pattern = species_a, raster = random_walk)

results_habitat_association(pattern = reconstruction, raster = landscape_discrete)

## ----plotting, eval = FALSE----------------------------------------------
#  plot_randomized_raster(random_walk)
#  
#  plot_randomized_pattern(reconstruction)

## ----energy, message = FALSE---------------------------------------------
calculate_energy(pattern = gamma_test, return_mean = TRUE)

calculate_energy(pattern = reconstruction, return_mean = TRUE)

