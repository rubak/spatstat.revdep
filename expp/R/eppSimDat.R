#' Type I error rate simulations
#' 
#' A helper function to perform Type I error rate simulations.
#' 
#' All default values match the values found in one of our study populations
#' ('Westerholz').
#' 
#' @param N Number of breeding pairs; default value is 10
#' @param meanClutch Mean clutch size (integer); clutch size it is assumed to
#' be Poisson distributed; default is 10
#' @param eppRate Proportion of extra-pair young in population; default is 0.10
#' @param eppMax Maximum number of extra-pair young by male; default is 12
#' @param eppMales Proportion of extra-pair males in population; default is
#' 0.35
#' @param nLags \code{maxlag} parameter to pass to
#' \link[expp]{DirichletPolygons}
#' @return An object of class 
#' \link[expp]{epp} The \code{data.frame} of the
#' \code{EPP} slot contains two variable (\code{trait_MALE} \code{trait_FEMALE}
#' ) simulated independent from the \code{epp} variable.
#' @examples
#' 
#' d = eppSimDat()
#' plot(d)
#' 
#' 
#' \donttest{
#' # Type I error rate simulation
#' 
#' require(lme4)
#' pval_glmer = vector(mode = "numeric", length = 0)
#' pval_glm = vector(mode = "numeric", length = 0)
#' 
#' # For meaningful results increase i to e.g. 500 and N in eppSimDat to e.g. 120
#' for(i in 1:5) { 
#'   x = as.data.frame(eppSimDat(N = 25, meanClutch = 10, eppRate = 0.10, eppMax = 12, 
#'       eppMales = 0.35, nLags = 3))
#'   
#'   fm1glmer = glmer(epp ~ rank + trait_MALE + trait_FEMALE + (1 | male) + (1 | female) , 
#'   data = x, family = binomial, nAGQ =  0)
#'   fm0glmer = update(fm1glmer, epp ~ 1 + (1 | male) + (1 | female) )
#'   pval_glmer[i] = anova(fm0glmer, fm1glmer)$"Pr(>Chisq)"[2]
#'   
#'   fm1glm = glm(epp ~ rank + trait_MALE + trait_FEMALE  , data = x, family = binomial)
#'   fm0glm = update(fm1glm, epp ~ 1 )
#'   pval_glm[i] = anova(fm0glm, fm1glm, test = "Chisq")$"Pr(>Chi)"[2]
#'   
#'   print(i)
#'  }
#' 
#' # Type I error rate of glmer models
#' table(pval_glmer<0.05)[2]/length(pval_glmer)
#' 
#' 
#' # Type I error rate of the equivalent glm models
#' table(pval_glm<0.05)[2]/length(pval_glm)
#' 
#' 
#' }
#' 
#' 
#' @export eppSimDat
eppSimDat <-function(N = 10, meanClutch = 10, eppRate = 0.10, eppMax = 12, eppMales = 0.35, nLags = 3) {
  
  # breeding
  d = data.frame(x = rnorm(N), y = rnorm(N), id = 1:N, male = paste0("m", 1:N), female = paste0("f", 1:N), stringsAsFactors = F )
  d$clutch = rpois(N, meanClutch)
  d$trait = rnorm(N)
  
  # epp
  d$epy = round(d$clutch * rnorm(N, mean = eppRate, sd = eppRate)    )
  d$epy = ifelse(d$epy > d$clutch, d$clutch, ifelse(d$epy < 0, 0, d$epy) )
  
  d$epmale = FALSE
  d[ sample(1:N, round(N*eppMales) ), "epmale"] = TRUE
  
  # epp pairs
  females = rep(d$female,  times = d$epy)
  males = d[d$epmale, "male"]
  males = sapply(females, function(f)    sample( setdiff(males, d[d$female == f, "male"] ), 1) )
  
  eppPairs = data.frame(female = females, male = males)
  
  e= eppMatrix(eppPairs, pairs = ~ male + female)
  
  # the epp object
  d = SpatialPointsBreeding(d[, c("x", "y", "male", "female", "id", "trait")], id = 'id')
  
  epp(d, DirichletPolygons(d) , e, maxlag = nLags)
  
  
  }

