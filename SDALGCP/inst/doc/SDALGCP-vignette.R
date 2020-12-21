## ----setup, include = FALSE---------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  eval= FALSE
)

## -----------------------------------------------------------------------------
#  require(SDALGCP)

## -----------------------------------------------------------------------------
#  data("PBCshp")

## -----------------------------------------------------------------------------
#  data <- as.data.frame(PBCshp@data)

## -----------------------------------------------------------------------------
#  data("pop_den")

## -----------------------------------------------------------------------------
#  pop_den[is.na(pop_den[])] <- 0

## -----------------------------------------------------------------------------
#  FORM <- X ~ propmale + Income + Employment + Education + Barriers + Crime +
#    Environment +  offset(log(pop))

## -----------------------------------------------------------------------------
#  phi <- seq(500, 1700, length.out = 20)

## ---- results = "hide",  warning = FALSE, message = FALSE---------------------
#  my_est <- SDALGCPMCML(data=data, formula=FORM, my_shp=PBCshp, delta=200, phi=phi, method=1, pop_shp=pop_den,
#                        weighted=TRUE, par0=NULL, control.mcmc=NULL)

## ---- results = "hide",  warning = FALSE, message = FALSE---------------------
#  summary(my_est)
#  #and for confidence interval use
#  confint(my_est)

## ---- results = "hide",  warning = FALSE, message = FALSE---------------------
#  phiCI(my_est, coverage = 0.95, plot = TRUE)

## ---- results = "hide",  warning = FALSE, message = FALSE---------------------
#  Dis_pred <- SDALGCPPred(para_est=my_est,  continuous=FALSE)

## ---- results = "hide",  warning = FALSE, message = FALSE---------------------
#  #to map the incidence
#  plot(Dis_pred, type="incidence", continuous = FALSE)
#  #and its standard error
#  plot(Dis_pred, type="SEincidence", continuous = FALSE)
#  #to map the covariate adjusted relative risk
#  plot(Dis_pred, type="CovAdjRelRisk", continuous = FALSE)
#  #and its standard error
#  plot(Dis_pred, type="SECovAdjRelRisk", continuous = FALSE)
#  #to map the exceedance probability that the covariate-adjusted relative risk is greter than a particular threshold
#  plot(Dis_pred, type="CovAdjRelRisk", continuous = FALSE, thresholds=3.0)

## ---- results = "hide",  warning = FALSE, message = FALSE---------------------
#  Con_pred <- SDALGCPPred(para_est=my_est, cellsize = 300, continuous=TRUE)

## ---- results = "hide",  warning = FALSE, message = FALSE---------------------
#  #to map the covariate adjusted relative risk
#  plot(Con_pred, type="relrisk")
#  #and its standard error
#  plot(Con_pred, type="SErelrisk")
#  #to map the exceedance probability that the relative risk is greter than a particular threshold
#  plot(Con_pred, type="relrisk", thresholds=1.5)

## ---- results = "hide",  warning = FALSE, message = FALSE---------------------
#  my_est <- SDALGCPMCML(data=data, formula=FORM, my_shp=PBCshp, delta=200, phi=phi, method=1,
#                        weighted=FALSE,  plot=FALSE, par0=NULL, control.mcmc=NULL, messages = TRUE, plot_profile = TRUE)

## ---- results = "hide",  warning = FALSE, message = FALSE---------------------
#  require(rgdal)
#  require(sp)
#  require(spacetime)
#  ohiorespMort <- read.csv("https://raw.githubusercontent.com/olatunjijohnson/dataset/master/OhioRespMort.csv")
#  download.file("https://github.com/olatunjijohnson/dataset/raw/master/ohio_shapefile.zip", "ohio_shapefile.zip")
#  unzip("ohio_shapefile.zip")
#  ohio_shp <- rgdal::readOGR("ohio_shapefile/","tl_2010_39_county00")
#  ohio_shp <- sp::spTransform(ohio_shp, sp::CRS("+init=epsg:32617"))

## ---- results = "hide",  warning = FALSE, message = FALSE---------------------
#  m <- length(ohio_shp)
#  TT <- 21
#  Y <- ohiorespMort$y
#  X <- ohiorespMort$year
#  pop <- ohiorespMort$n
#  E <- ohiorespMort$E
#  data <- data.frame(Y=Y, X=X, pop=pop, E=E)
#  formula <- Y ~  X + offset(log(E))
#  phi <- seq(10, 300, length.out = 10)
#  control.mcmc <- controlmcmcSDA(n.sim=10000, burnin=2000, thin=80, h=1.65/((m*TT)^(1/6)), c1.h=0.01, c2.h=0.0001)
#  time <- as.POSIXct(paste(1968:1988, "-01-01", sep = ""), tz = "")
#  st_data <- spacetime::STFDF(sp = ohio_shp, time = time, data = data)

## ---- results = "hide",  warning = FALSE, message = FALSE---------------------
#  spacetime::stplot(st_data[,,"Y"])

## ---- results = "hide",  warning = FALSE, message = FALSE---------------------
#  model.fit <- SDALGCPMCML_ST(formula=formula, st_data = st_data,  delta=800,
#                              phi=phi, method=2, pop_shp=NULL,  kappa=0.5,
#                              weighted=FALSE, par0=NULL, control.mcmc=control.mcmc,
#                              plot=TRUE, plot_profile=TRUE, rho=NULL,
#                              giveup=50, messages=TRUE)
#  summary(model.fit)

## ---- results = "hide",  warning = FALSE, message = FALSE---------------------
#  dis_pred <- SDALGCPPred_ST(para_est = model.fit, continuous = FALSE)

## ---- results = "hide",  warning = FALSE, message = FALSE---------------------
#  plot(dis_pred, type="CovAdjRelRisk", main="Relative Risk", continuous=FALSE)
#  plot(dis_pred,  type="incidence", main="Incidence", continuous=FALSE)

## ---- results = "hide",  warning = FALSE, message = FALSE---------------------
#  con_pred <- SDALGCPPred_ST(para_est = model.fit, cellsize = 2500, continuous=TRUE, n.window = 1)

## ---- results = "hide",  warning = FALSE, message = FALSE---------------------
#  plot(con_pred, type="relrisk", continuous=TRUE)

