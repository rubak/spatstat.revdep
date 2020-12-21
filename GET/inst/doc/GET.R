### R code from vignette source 'GET.Rnw'

###################################################
### code chunk number 1: preliminaries
###################################################
options(prompt = "R> ", continue = "+  ", width = 70, useFancyQuotes = FALSE)
library("ggplot2")
library("patchwork")
library("GET")
theme_set(theme_bw(base_size = 9))


###################################################
### code chunk number 2: girls_data
###################################################
library("fda")
years <- paste(1:18)
curves <- growth[['hgtf']][years,]
cset1 <- create_curve_set(list(r = as.numeric(years),
  obs = curves))
cset2 <- create_curve_set(list(r = as.numeric(years[-1]),
  obs = curves[-1,] - curves[-nrow(curves),]))


###################################################
### code chunk number 3: girls_forder
###################################################
A1 <- forder(cset1, measure = 'area'); order(A1)[1:10]
A2 <- forder(cset2, measure = 'area'); order(A2)[1:10]


###################################################
### code chunk number 4: girls_combined_forder
###################################################
csets <- list(Height = cset1, Change = cset2)
A <- forder(csets, measure = 'area'); order(A)[1:10]


###################################################
### code chunk number 5: girls_curves
###################################################
col3 <- c("#21908CFF", "#440154FF", "#5DC863FF")
p1 <- plot(cset1, idx = order(A)[1:3], col_idx = col3) +
  labs(x = "Age (years)", y = "Height")
p2 <- plot(cset2, idx = order(A)[1:3], col_idx = col3) +
  labs(x = "Age (years)", y = "Change")
p1 + p2 + plot_layout(guides = "collect")


###################################################
### code chunk number 6: girls_fboxplot
###################################################
res <- fBoxplot(csets, type = 'area', factor = 1.5)
plot(res) + labs(x = "Age (years)", y = "Value")


###################################################
### code chunk number 7: girls_fboxplot_plot
###################################################
plot(res) + labs(x = "Age (years)", y = "Value")


###################################################
### code chunk number 8: adulttrees
###################################################
data("adult_trees")
ggplot(adult_trees) + geom_point(aes(x=x, y=y)) +
  xlim(c(0, 75)) + ylim(c(0, 75)) + coord_cartesian(expand = FALSE)


###################################################
### code chunk number 9: adulttrees_seed
###################################################
set.seed(190503)


###################################################
### code chunk number 10: adulttrees_sims
###################################################
library("spatstat.core")
data("adult_trees")
X <- as.ppp(adult_trees, W = square(75))
nsim <- 999
obs.L <- Lest(X, correction = "translate")
r <- obs.L[['r']]
obs <- obs.L[['trans']] - r
sim <- matrix(nrow = length(r), ncol = nsim)
for(i in 1:nsim) {
  sim.X <- runifpoint(ex = X)
  sim[, i] <- Lest(sim.X, correction = "translate", r = r)[['trans']] - r
}


###################################################
### code chunk number 11: adulttrees_curve_set
###################################################
cset <- create_curve_set(list(r = r, obs = obs, sim_m = sim))


###################################################
### code chunk number 12: adulttrees_CSR
###################################################
res <- global_envelope_test(cset, type = "erl")
plot(res) + ylab(expression(italic(hat(L)(r)-r)))


###################################################
### code chunk number 13: adulttrees_CSR_plot
###################################################
plot(res, ylab = expression(italic(hat(L)(r)-r)))


###################################################
### code chunk number 14: adulttrees_seed2
###################################################
set.seed(190503)


###################################################
### code chunk number 15: adulttrees_envelope
###################################################
env <- envelope(X, nsim = 999, fun = "Lest", correction = "translate",
  transform = expression(.-r), simulate = expression(runifpoint(ex=X)),
  savefuns = TRUE, verbose = FALSE)
res <- global_envelope_test(env, type = "erl")


###################################################
### code chunk number 16: adjtest_data
###################################################
library("fda.usc")
data("poblenou")
dat <- poblenou[['nox']][['data']][,'H10']
n <- length(dat)


###################################################
### code chunk number 17: adjtest_seed
###################################################
set.seed(200127)


###################################################
### code chunk number 18: adjtest_est
###################################################
mu <- mean(dat)
sigma <- sd(dat)


###################################################
### code chunk number 19: adjtest_sim1
###################################################
nsim <- nsimsub <- 199 # The number of simulations
r <- seq(min(dat), max(dat), length=100)
obs <- stats::ecdf(dat)(r)
sim <- sapply(1:nsimsub, function(i) {
  x <- rnorm(n, mean = mu, sd = sigma)
  stats::ecdf(x)(r)
})
cset <- create_curve_set(list(r = r, obs = obs, sim_m = sim))


###################################################
### code chunk number 20: adjtest_sim2
###################################################
cset.ls <- list()
for(rep in 1:nsim) {
  x <- rnorm(n, mean = mu, sd = sigma)
  mu2 <- mean(x)
  sigma2 <- sd(x)
  obs2 <- stats::ecdf(x)(r)
  sim2 <- sapply(1:nsimsub, function(i) {
    x2 <- rnorm(n, mean = mu2, sd = sigma2)
    stats::ecdf(x2)(r)
  })
  cset.ls[[rep]] <- create_curve_set(list(r = r, obs = obs2,
    sim_m = sim2))
}


###################################################
### code chunk number 21: adjtest_get
###################################################
res <- GET.composite(X=cset, X.ls=cset.ls, type='erl')
plot(res) + labs(x = "NOx", y = "Ecdf")


###################################################
### code chunk number 22: adjtest_log
###################################################
p1 <- plot(res) + labs(x = "NOx", y = "Ecdf")
# For log(nox) values:
dat <- log(poblenou[['nox']][['data']][,'H10'])
set.seed(200127)
# 1. Fit the model
mu <- mean(dat)
sigma <- sd(dat)
# 2. Simulate a sample from the fitted null model and
#    compute the test vectors for data (obs) and each simulation (sim)
nsim <- nsimsub <- 199
r <- seq(min(dat), max(dat), length=100)
obs <- stats::ecdf(dat)(r)
sim <- sapply(1:nsimsub, function(i) {
  x <- rnorm(n, mean = mu, sd = sigma)
  stats::ecdf(x)(r)
})
cset <- create_curve_set(list(r = r, obs = obs, sim_m = sim))
# 3. Simulate another sample from the fitted null model.
# 4. Fit the null model to each of the patterns,
#    simulate a sample from the null model,
#    and compute the test vectors for all.
cset.ls <- list()
for(rep in 1:nsim) {
  x <- rnorm(n, mean = mu, sd = sigma)
  mu2 <- mean(x)
  sigma2 <- sd(x)
  obs2 <- stats::ecdf(x)(r)
  sim2 <- sapply(1:nsimsub, function(i) {
    x2 <- rnorm(n, mean = mu2, sd = sigma2)
    stats::ecdf(x2)(r)
  })
  cset.ls[[rep]] <- create_curve_set(list(r = r, obs = obs2,
    sim_m = sim2))
}
# Perform the adjusted test
res <- GET.composite(X=cset, X.ls=cset.ls, type='erl')
p2 <- plot(res) + labs(x = "NOx", y = "Ecdf")


###################################################
### code chunk number 23: normalitytest
###################################################
combined <- p1 + p2 & theme(legend.position = "bottom")
combined + plot_layout(guides = "collect")


###################################################
### code chunk number 24: necdfs_10_seed
###################################################
set.seed(190503)


###################################################
### code chunk number 25: necdfs_10
###################################################
fm10.l <- list(Girls=growth$hgtf["10",], Boys=growth$hgtm["10",])
res10 <- GET.necdf(fm10.l, nsim=1999)
plot(res10)


###################################################
### code chunk number 26: necdfs_14
###################################################
set.seed(190503)
fm14.l <- list(Girls=growth$hgtf["14",], Boys=growth$hgtm["14",])
res14 <- GET.necdf(fm14.l, nsim=1999)


###################################################
### code chunk number 27: necdfs
###################################################
nG <- length(fm10.l$Girls); nB <- length(fm10.l$Boys)
df <- data.frame(Height = c(fm10.l$Girls, fm14.l$Girls, fm10.l$Boys, fm14.l$Boys),
  Gender = c(rep("Girls", times=2*nG), rep("Boys", times=2*nB)),
  Age = c(rep("Age 10", times=nG), rep("Age 14", times=nG),
          rep("Age 10", times=nB), rep("Age 14", times=nB)))
ggplot(df, aes(Height, colour = Gender)) + facet_wrap("Age") + stat_ecdf() +
  scale_color_manual(values = c("#440154FF", "#5DC863FF")) +
  xlab(substitute(paste(italic(i), " (", j, ")", sep=""),
    list(i="x", j="Height in cm")))


###################################################
### code chunk number 28: necdfs_means_GET
###################################################
p1 <- plot(res10)
p2 <- plot(res14)
combined <- p1 + p2 & theme(legend.position = "bottom")
combined + plot_layout(guides = "collect")


###################################################
### code chunk number 29: graphfanova_data
###################################################
library("fda.usc")
data("poblenou")
fest <- poblenou$df$day.festive; week <- as.integer(poblenou$df$day.week)
Type <- vector(length=length(fest))
Type[fest == 1 | week >= 6] <- "Free"
Type[fest == 0 & week %in% 1:4] <- "MonThu"
Type[fest == 0 & week == 5] <- "Fri"
Type <- factor(Type, levels = c("MonThu", "Fri", "Free"))


###################################################
### code chunk number 30: NOx_data
###################################################
maxs <- apply(poblenou[['nox']]$data, MARGIN = 1, FUN = max)
df <- do.call(rbind, lapply(1:24, FUN = function(x) {
  data.frame(Hour = x, NOx = poblenou[['nox']]$data[,x], Type = Type,
    Date = rownames(poblenou[['nox']]$data), Max = maxs)
}))
ggplot(df) + geom_line(aes(x = Hour, y = NOx, group = Date, col = Max)) +
  facet_wrap(vars(Type)) + scale_color_continuous(type = "viridis")


###################################################
### code chunk number 31: graphfanova_cset
###################################################
cset <- create_curve_set(list(r=0:23,
  obs=t(log(poblenou[['nox']][['data']]))))


###################################################
### code chunk number 32: graphfanova_seed
###################################################
set.seed(190503)


###################################################
### code chunk number 33: graphfanova_logNOx
###################################################
res.c <- graph.fanova(nsim = 2999, curve_set = cset, groups = Type,
  variances = "unequal", contrasts = TRUE)
plot(res.c)


###################################################
### code chunk number 34: abide_9002_23_subj1and27
###################################################
data("abide_9002_23")
plot(abide_9002_23$curve_set, idx=c(1, 27))


###################################################
### code chunk number 35: flm_graph_seed
###################################################
set.seed(202001)


###################################################
### code chunk number 36: flm_graph
###################################################
res <- graph.flm(nsim = 999, formula.full = Y ~ Group + Sex + Age,
  formula.reduced = Y ~ Sex + Age,
  curve_sets = list(Y = abide_9002_23[['curve_set']]),
  factors = abide_9002_23[['factors']], contrasts = TRUE,
  GET.args = list(type = "area"))


###################################################
### code chunk number 37: flm_graph_plot
###################################################
plot(res)


###################################################
### code chunk number 38: flm_graph_abide
###################################################
plot(res)


###################################################
### code chunk number 39: flm_frank_seed
###################################################
set.seed(202001)


###################################################
### code chunk number 40: flm_frank
###################################################
res.F <- frank.flm(nsim = 999, formula.full = Y ~ Group + Age + Sex,
  formula.reduced = Y ~ Age + Sex,
  curve_sets = list(Y = abide_9002_23[['curve_set']]),
  factors = abide_9002_23[['factors']], GET.args = list(type = "area"))
plot(res.F)


###################################################
### code chunk number 41: flm_frank_abide
###################################################
plot(res.F)


###################################################
### code chunk number 42: polynomial
###################################################
set.seed(190504)
# Polynomial regression from Narisetty & Nair (2016)
# Simulate regression data according to the cubic model
# f(x) = 0.8x - 1.8x^2 + 1.05x^3 for x in [0,1]
par <- c(0,0.8,-1.8,1.05) # Parameters of the true polynomial model
res <- 100 # Resolution
x <- seq(0, 1, by=1/res); x2=x^2; x3=x^3;
f <- par[1] + par[2]*x + par[3]*x^2 + par[4]*x^3 # The true function
d <- f + rnorm(length(x), 0, 0.04) # Data

# Estimate polynomial regression model
reg <- lm(d ~ x + x2 + x3)
ftheta <- reg$fitted.values
resid0 <- reg$residuals
s0 <- sd(resid0)

# Bootstrap regression
B <- 2000 # Number of bootstrap samples
ftheta1 <- array(0, c(B,length(x)))
s1 <- array(0,B)
for(i in 1:B) {
  u <- sample(resid0, size=length(resid0), replace=TRUE)
  reg1 <- lm((ftheta+u) ~ x + x2 + x3)
  ftheta1[i,] <- reg1$fitted.values
  s1[i] <- sd(reg1$residuals)
}

# Centering and scaling
meanftheta <- apply(ftheta1, 2, mean)
m <- array(0, c(B,length(x)))
for(i in 1:B) { m[i,] <- (ftheta1[i,]-meanftheta)/s1[i] }

# Central region computation
boot.cset <- create_curve_set(list(r=1:length(x), obs=ftheta+s0*t(m)))
cr <- central_region(boot.cset, coverage=0.95, type="erl")
plot(cr) + labs(x = expression(italic(x)), y = expression(italic(f(x)))) +
  geom_point(data = data.frame(id = 1:length(d), points = d),
    aes(x = id, y = points)) + # data points
  geom_line(data = data.frame(id = 1:length(d), points = f),
    aes(x = id, y = points)) # true function


