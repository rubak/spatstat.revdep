## ---- include = FALSE---------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----setup, message=FALSE-----------------------------------------------------
library(stats19)
library(dplyr)

## -----------------------------------------------------------------------------
v = get_stats19(year = 2018, type = "vehicles")
v

## -----------------------------------------------------------------------------
v = v %>% mutate(vehicle_type2 = case_when(
  grepl(pattern = "motorcycle", vehicle_type, ignore.case = TRUE) ~ "Motorbike",
  grepl(pattern = "Car", vehicle_type, ignore.case = TRUE) ~ "Car",
  grepl(pattern = "Bus", vehicle_type, ignore.case = TRUE) ~ "Bus",
  grepl(pattern = "cycle", vehicle_type, ignore.case = TRUE) ~ "Cycle",
  # grepl(pattern = "Van", vehicle_type, ignore.case = TRUE) ~ "Van",
  grepl(pattern = "Goods", vehicle_type, ignore.case = TRUE) ~ "Goods",
  
  TRUE ~ "Other"
))
# barplot(table(v$vehicle_type2))

## -----------------------------------------------------------------------------
table(v$vehicle_type2)
summary(v$age_of_driver)
summary(v$engine_capacity_cc)
table(v$propulsion_code)
summary(v$age_of_vehicle)

## -----------------------------------------------------------------------------
a = get_stats19(year = 2018, type = "accidents")
va = dplyr::inner_join(v, a)

## -----------------------------------------------------------------------------
dim(v)
dim(va)
names(va)

## ---- out.width="100%"--------------------------------------------------------
xtabs(~vehicle_type2 + accident_severity, data = va) %>% prop.table()
xtabs(~vehicle_type2 + accident_severity, data = va) %>% prop.table() %>% plot()

## -----------------------------------------------------------------------------
vac = va %>% filter(vehicle_type2 == "Car")

## -----------------------------------------------------------------------------
summary(vac$engine_capacity_cc)

## ---- echo=FALSE, eval=FALSE--------------------------------------------------
#  library(tidyverse)
#  
#  max_engine_size = 5000
#  min_engine_size = 300
#  
#  sel_too_big = vac$engine_capacity_cc > max_engine_size
#  sel_too_small = vac$engine_capacity_cc < min_engine_size
#  sum(sel_too_big) / nrow(vac)
#  sum(sel_too_small) / nrow(vac)
#  vac$engine_capacity_cc[sel_too_big | sel_too_small] = NA

## ---- eval=FALSE, echo=FALSE--------------------------------------------------
#  
#  vac = vac %>% filter(val_size)
#  vac %>%
#    mutate(age = formatC(age_band_of_driver, digits = 2, flag = "0")) %>%
#    ggplot() +
#    geom_violin(aes(age, engine_capacity_cc))
#  
#  vac$sev_factor = factor(vac$accident_severity, labels = 3:1)
#  vac$sev_numeric = vac$sev_factor %>% as.character() %>%  as.numeric()
#  summary(vac$sev_factor)
#  summary(vac$sev_numeric)
#  
#  m = lm(sev_numeric ~ engine_capacity_cc + age_of_driver + speed_limit, data = vac)
#  summary(m)

## ---- echo=FALSE, eval=FALSE--------------------------------------------------
#  table_vehicle_type = xtabs(cbind(accident_severity, vehicle_type) ~ accident_severity, data = va)
#  group_totals = va %>%
#    group_by(accident_severity) %>%
#    summarise(n = n())
#  
#  # fails
#  fit = glm(data = va, factor(accident_severity) ~
#              engine_capacity_cc +
#              age_of_driver +
#              engine_capacity_cc +
#              factor(propulsion_code) +
#              age_of_vehicle
#            )
#  # works but result does not have probabilities
#  ?nnet::multinom
#  mod = nnet::multinom(formula = accident_severity ~
#              engine_capacity_cc +
#              age_of_driver +
#              engine_capacity_cc +
#              propulsion_code +
#              age_of_vehicle, data = va)
#  mod
#  summary(mod)
#  mod$nunits
#  class(mod$fitted.values)
#  dim(mod$fitted.values)
#  colnames(mod$fitted.values)
#  summary(mod$fitted.values) # result!
#  
#  probs = as.data.frame(mod$fitted.values)
#  head(probs)
#  head(rowSums(probs))
#  colSums(probs) / group_totals$n
#  nrow(probs)
#  nrow(va)
#  
#  # install.packages("mlogit")
#  install.packages("AER")
#  vignette(package = "mlogit")
#  vignette("c2.formula.data")
#  library(mlogit)
#  data("TravelMode", package = "AER")
#  head(TravelMode)
#  ?TravelMode
#  summary(TravelMode$choice)
#  summary(TravelMode$mode)
#  TM = mlogit.data(TravelMode, choice = "choice", shape = "long",
#                   alt.levels = c("air", "train", "bus", "car"))
#  
#  vamld = mlogit.data(va, choice = "accident_severity", alt.levels = c("Slight", "Serious", "Fatal"), shape = "wide")
#  
#  mlogit(accident_severity ~ speed_limit | 0, vamld[1:999, ])
#  
#  vamld = mlogit::mlogit.data(va, choice = "accident_severity", shape = "wide")
#  head(vamld)
#  
#  
#  # fails
#  # vm = mlogit(accident_severity ~ engine_capacity_cc + speed_limit | 0, vamld[1:999, ])
#  # apply(fitted(vm, outcome = FALSE), 2, mean)
#  # vm

