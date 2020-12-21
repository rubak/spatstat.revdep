### R code from vignette source 'pointpatterns.Rnw'

###################################################
### code chunk number 1: preliminaries
###################################################
options(prompt = "R> ", continue = "+  ", width = 70, useFancyQuotes = FALSE)
set.seed(123)


###################################################
### code chunk number 2: r_setup
###################################################
library("GET")
library("spatstat")
library("ggplot2")
theme_set(theme_bw(base_size = 9))


###################################################
### code chunk number 3: workflow1
###################################################
X <- spruces
X


###################################################
### code chunk number 4: workflow2
###################################################
env <- envelope(X, nsim=999, savefuns=TRUE,
                simulate=expression(runifpoint(ex=X)),
                verbose=FALSE)


###################################################
### code chunk number 5: workflow3
###################################################
res <- global_envelope_test(env)


###################################################
### code chunk number 6: workflow4
###################################################
plot(res)


###################################################
### code chunk number 7: spruces
###################################################
X <- unmark(spruces)
plot(X)


###################################################
### code chunk number 8: spruces_csr
###################################################
nsim <- 1999 # Number of simulations
env <- envelope(X, fun="Lest", nsim=nsim,
                savefuns=TRUE, # save the functions
                correction="translate", # edge correction for L
                transform = expression(.-r), # centering
                simulate=expression(runifpoint(ex=X)), # Simulate CSR
                verbose=FALSE)
# The rank envelope test (ERL)
res <- global_envelope_test(env, type="erl")
# Plot the result
plot(res)


###################################################
### code chunk number 9: spruces_csr_continues
###################################################
cset <- crop_curves(env, r_min=1, r_max=7)
# Do the rank envelope test (erl)
res <- global_envelope_test(cset, type="erl")
plot(res, ylab=expression(italic(hat(L)(r)-r)))


###################################################
### code chunk number 10: spruces_mpp
###################################################
mpp <- spruces
plot(mpp)


###################################################
### code chunk number 11: spruces_randomlabeling
###################################################
nsim <- 1999 # Number of simulations
env <- envelope(mpp, fun=Kmark, nsim = nsim, f=function(m1, m2) { m1*m2 },
        correction="translate", returnL=TRUE,
        simulate=expression(rlabel(mpp, permute=TRUE)), # Permute the marks
        savefuns=TRUE, # Save the functions
        verbose=FALSE)


###################################################
### code chunk number 12: spruces_randomlabeling2
###################################################
# Crop curves to desired r-interval
curve_set <- crop_curves(env, r_min=1.5, r_max=9.5)
# Center the functions, i.e. take \hat{L}_mm(r)-the mean of simulated functions.
curve_set <- residual(curve_set)
# The global envelope test
res <- global_envelope_test(curve_set)
plot(res, ylab=expression(italic(L[mm](r)-L(r))))


###################################################
### code chunk number 13: combined1
###################################################
data(saplings)
X <- as.ppp(saplings, W = square(75))

nsim <- 499 # Number of simulations

# Specify distances for different test functions
n <- 500 # the number of r-values
rmin <- 0; rmax <- 20; rstep <- (rmax-rmin)/n
rminJ <- 0; rmaxJ <- 8; rstepJ <- (rmaxJ-rminJ)/n
r <- seq(0, rmax, by=rstep)    # r-distances for Lest
rJ <- seq(0, rmaxJ, by=rstepJ) # r-distances for Fest, Gest, Jest


###################################################
### code chunk number 14: combined2
###################################################
env_L <- envelope(X, nsim=nsim,
          simulate=expression(runifpoint(ex=X)),
          fun="Lest", correction="translate",
          transform=expression(.-r), # Take the L(r)-r function instead of L(r)
          r=r,                       # Specify the distance vector
          savefuns=TRUE,             # Save the estimated functions
          savepatterns=TRUE,         # Save the simulated patterns
          verbose=FALSE)


###################################################
### code chunk number 15: combined3
###################################################
simulations <- attr(env_L, "simpatterns")


###################################################
### code chunk number 16: combined4
###################################################
env_F <- envelope(X, nsim=nsim,
                simulate=simulations,
                fun="Fest", correction="Kaplan", r=rJ,
                savefuns=TRUE, verbose=FALSE)
env_G <- envelope(X, nsim=nsim,
                simulate=simulations,
                fun="Gest", correction="km", r=rJ,
                savefuns=TRUE, verbose=FALSE)
env_J <- envelope(X, nsim=nsim,
                simulate=simulations,
                fun="Jest", correction="none", r=rJ,
                savefuns=TRUE, verbose=FALSE)


###################################################
### code chunk number 17: combined5
###################################################
curve_set_L <- crop_curves(env_L, r_min=rmin, r_max=rmax)
curve_set_F <- crop_curves(env_F, r_min=rminJ, r_max=rmaxJ)
curve_set_G <- crop_curves(env_G, r_min=rminJ, r_max=rmaxJ)
curve_set_J <- crop_curves(env_J, r_min=rminJ, r_max=rmaxJ)


###################################################
### code chunk number 18: combined6
###################################################
res <- global_envelope_test(curve_sets=list(curve_set_L, curve_set_F,
                                          curve_set_G, curve_set_J))
plot(res, labels=c("L(r)-r", "F(r)", "G(r)", "J(r)"))


###################################################
### code chunk number 19: goodness-of-fit1
###################################################
X <- unmark(spruces)
# Minimum distance between points in the pattern
min(nndist(X))
# Fit a model
fittedmodel <- ppm(X, interaction=Hardcore(hc=1)) # Hardcore process


###################################################
### code chunk number 20: goodness-of-fit2
###################################################
#env <- envelope(fittedmodel, fun="Jest", nsim=999, savefuns=TRUE,
#                correction="none", r=seq(0, 4, length=500))


###################################################
### code chunk number 21: goodness-of-fit3
###################################################
simulations <- NULL
nsim <- 999 # Number of simulations
for(j in 1:nsim) {
   simulations[[j]] <- rHardcore(beta=exp(fittedmodel$coef[1]),
                                 R=fittedmodel$interaction$par$hc,
                                 W=X$window)
}
env <- envelope(X, simulate=simulations, fun="Jest",
                nsim=length(simulations),
                savefuns=TRUE, correction="none",
                r=seq(0, 4, length=500),
                verbose=FALSE)
curve_set <- crop_curves(env, r_min=1, r_max=3.5)
res <- global_envelope_test(curve_set, type="erl")
plot(res, ylab=expression(italic(J(r))))


###################################################
### code chunk number 22: saplings
###################################################
data(saplings)
saplings <- as.ppp(saplings, W = square(75))
plot(saplings)


###################################################
### code chunk number 23: saplings_adjusted_init
###################################################
rmin <- 0.3; rmax <- 10; rstep <- (rmax-rmin)/500
r <- seq(0, rmax, by=rstep)
nsim <- 19 # Increase nsim for serious analysis!


###################################################
### code chunk number 24: saplings_adjusted_fit
###################################################
M1 <- kppm(saplings~1, clusters = "MatClust", statistic="K")


###################################################
### code chunk number 25: saplings_adjusted_test
###################################################
adjenvL <- GET.composite(X = M1, nsim = nsim,
            testfuns = list(L = list(fun="Lest", correction="translate",
                 transform = expression(.-r), r=r)), # passed to envelope
            type = "area", r_min=rmin, r_max=rmax, verbose=FALSE)
plot(adjenvL)


###################################################
### code chunk number 26: spatialF
###################################################
data(bei)


###################################################
### code chunk number 27: spatialF_bei
###################################################
fullmodel <- ~ grad
reducedmodel <- ~ 1
fitppm <- function(X, model, covariates) {
  ppm(X, model, covariates=covariates)
}

nsim <- 19 # Increase nsim for serious analysis!
res <- GET.spatialF(bei, fullmodel, reducedmodel,
                    fitppm, bei.extra, nsim)

plot(res$F)
plot(res$S)


###################################################
### code chunk number 28: spatialF_forestfires
###################################################
# Example of forest fires
data("clmfires")
# Choose the locations of the lightnings in years 2004-2007:
pp.lightning <- unmark(subset(clmfires, cause == "lightning" &
                 date >= "2004-01-01" & date < "2008-01-01"))

covariates <- clmfires.extra$clmcov100
covariates$forest <- 
  covariates$landuse == "conifer" | covariates$landuse == "denseforest" |
  covariates$landuse == "mixedforest"

fullmodel <- ~ elevation + landuse
reducedmodel <- ~ landuse
nsim <- 19 # Increase nsim for serious analysis!
res <- GET.spatialF(pp.lightning, fullmodel, reducedmodel,
                    fitppm, covariates, nsim)
plot(res$F)
plot(res$S)


###################################################
### code chunk number 29: spatialF_extra
###################################################
# fitfun for the log Gaussian Cox Process with exponential covariance function
fitLGCPexp <- function(X, model, covariates) {
  kppm(X, model, clusters="LGCP", model="exponential", covariates=covariates)
}
# fitfun for the hardcore process with hardcore radius 0.01
fitHardcore <- function(X, model, covariates) {
  ppm(X, model, interaction=Hardcore(0.01), covariates = covariates)
}


###################################################
### code chunk number 30: saplings_data
###################################################
data(saplings)
saplings <- as.ppp(saplings, W = square(75))


###################################################
### code chunk number 31: saplings_init
###################################################
nr <- 500
rmin <- 0.3; rminJ <- 0.3
rmax <- 10; rmaxJ <- 6
rstep <- (rmax-rmin)/nr; rstepJ <- (rmaxJ-rminJ)/nr
r <- seq(0, rmax, by=rstep)
rJ <- seq(0, rmaxJ, by=rstepJ)


###################################################
### code chunk number 32: saplings_csr_L1
###################################################
nsim <- 999 # Number of simulations
env <- envelope(saplings, nsim=nsim,
  simulate=expression(runifpoint(saplings$n, win=saplings$window)), # Simulate CSR
  fun="Lest", correction="translate", # estimator of L with translational edge corr.
  transform = expression(.-r),        # Take the L(r)-r function instead of L(r)
  r=r,                                # Specify the distance vector
  savefuns=TRUE,                      # Save the estimated functions
  verbose=FALSE)


###################################################
### code chunk number 33: saplings_csr_L2
###################################################
curve_set <- crop_curves(env, r_min = rmin, r_max = rmax)


###################################################
### code chunk number 34: saplings_csr_L3
###################################################
res <- global_envelope_test(curve_set, type="erl")
plot(res, ylab=expression(italic(hat(L)(r)-r)))


###################################################
### code chunk number 35: saplings_matclust
###################################################
fitted_model <- kppm(saplings~1, clusters = "MatClust", statistic="K")


###################################################
### code chunk number 36: saplings_matclust2
###################################################
nsim <- 19 # 19 just for experimenting with the code!!
#nsim <- 499 # 499 is ok for type = 'qdir' (takes > 1 h)
adjenvL <- GET.composite(X = fitted_model,
                      fun="Lest", correction="translate",
                      transform = expression(.-r), r=r,
                      type = "qdir", nsim = nsim, nsimsub = nsim,
                      r_min=rmin, r_max=rmax, verbose=FALSE)
# Plot the test result
plot(adjenvL, ylab=expression(italic(L(r)-r)))


###################################################
### code chunk number 37: saplings_matclust3
###################################################
adjenvJ <- GET.composite(X = fitted_model,
                        fun="Jest", correction="none", r=rJ,
                        type = "qdir", nsim = nsim, nsimsub = nsim,
                        r_min=rminJ, r_max=rmaxJ, verbose=FALSE)
# Plot the test result
plot(adjenvJ, ylab=expression(italic(J(r))))


###################################################
### code chunk number 38: saplings_matclust4
###################################################
adjenvLJ <- GET.composite(X = fitted_model,
      testfuns = list(L = list(fun="Lest", correction="translate",
                          transform = expression(.-r), r=r),
                      J = list(fun="Jest", correction="none", r=rJ)),
      type = "erl", nsim = nsim, nsimsub = nsim,
      r_min=c(rmin, rminJ), r_max=c(rmax, rmaxJ),
      save.cons.envelope=TRUE, verbose=FALSE)
plot(adjenvLJ)


