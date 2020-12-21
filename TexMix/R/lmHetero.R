#' @title Function: Multiplicately weighted regression model
#'
#' @description
#' \code{lmHetero} accounts for heteroscedasticity in regression models
#'
#' @details
#' This function estimates the parameters of a regression model whose normally
#' distributed disturbances have a variance that multiplicatively depends on a
#' set of strictly positive weights variables. That is,
#'
#' \deqn{\sigma^2_i = \exp(\gamma_0 + \gamma_1 \cdot \log(z_{i1}) + ...)}
#'
#' The weights variables z must be entered in their logarithmic forms. The
#' paramater \eqn{\exp(\gamma_0)} expresses the global variance.
#'
#' @param formula Formula object (perhaps, multiple parts formula
#'  \code{y~x1+x2+...| log(z1)+log(z2)+...}) linking the dependent
#' variable to its set of associated independent variables x and a second
#' expression separated by \code{|} modelling the local variances with the
#' variables z. The formula \code{y~x1+x2+...} assumes homoscedasticity.
#' \strong{Important:} the weights variables z must be entered in
#' their log-transformed form! This is the preferred specification of the function call.
#'
#' @param hetero \strong{Valid only if} omitting the second \code{formula} term
#' \code{| log(z1)+log(z2)+...}. The default parameter
#' \code{hetero= ~1} assumes a constant variance (homoscedasticity), whereas
#' \code{hetero=~log(z1)+log(z2)+...} explicitly models heteroscedasticity.
#'
#'
#' @param hetero Optional formula specification without "\code{| log(z1)+log(z2)+...}":
#' A second parameter modeling the heteroscedasticity with a right-handed formula
#' by a set of variables z. That is, "\code{hetero=~log(z1)+log(z2)+...}".
#' Omitting the second parameter assumes "\code{hetero=~1}". This is only included
#' for backward compatibility).
#'
#' @param data An optional data frame containing the variables in the model.
#'  By default the variables are taken from the environment of the formula.
#'
#' @param subset An optional vector specifying a subset of observations to
#' be used in fitting the model.
#'
#' @param na.action A function that indicates what should happen when the data
#' contain \code{NA}s. The default is set by the "\code{na.action}" option
#'
#' @param contrasts An optional list. See the "\code{contrasts.arg}" of
#' \code{\link{model.matrix.default}}.
#'
#' @param iter Logical indicating whether the interation history should be
#' displayed. The default setting if "\code{FALSE}".
#'
#' @param ... Currently not in use.
#'
#' @export
#' @return a list with 10 elements:
#' \item{CALL}{function call}
#' \item{sigma2}{global variance estimate \eqn{exp(gamma_0)}}
#' \item{gamma}{vector of estimated gamma coefficients}
#' \item{namesGamma}{vector of variable names expressed by Z}
#' \item{beta}{vector of estimated weight adjusted regression parameters}
#' \item{weights}{vector of weights \eqn{1/\sigma^2_i} estimates for each
#' observation. It can be used in the call \code{lm(..., weights=weights)}
#' to adjust for heteroscedasticity}
#' \item{covBeta}{covariance matrix of the estimated regression coefficients}
#' \item{covGamma}{covariance matrix of the estimated gamma coefficients}
#' \item{logLikeH1}{log-likelihood of the heteroscedastic adjusted regression
#' model}
#' \item{logLikeH0}{log-likelihood of the unadjusted regression model}
#'
#' @source The maximum likelihood estimation procedure for multiplicately
#' weighted regression is given in Greene W. H. (2000). Econometric Analysis.
#' 4th edition. Upper Saddle River: Prentice Hall. pp 516-521
#' (Note: page numbers will differ for other editions)
#'
#' @author Michael Tiefelsdorf (\email{tiefelsdorf@@utdallas.edu}) & Yongwan Chun
#'
#' @examples
#' library(sp)
#' data(tractShp)
#' validTractShp <- tractShp[!is.na(tractShp$BUYPOW), ]         # Remove 2 tracts with NA's
#' ## Population at risk
#' totPop <- validTractShp$MALE+validTractShp$FEMALE
#'
#' ## H0 model (homoscedasticity)
#' mod.lm <- mod.lmH <- lmHetero(PERCAPINC~PCTNOHINS+PCTMINOR+PCTUNIVDEG+PCTWHITE,
#'                               data=validTractShp)
#' summary(mod.lm)
#'
#' ## Preferred heteroscedasticiy function call
#' mod.lmH <- lmHetero(PERCAPINC~PCTNOHINS+PCTMINOR+PCTUNIVDEG+PCTWHITE|log(totPop),
#'                     data=validTractShp)
#' summary(mod.lmH)
#'
#' ## Alternative equivalent heteroscedasticiy function call
#' mod.lmH <- lmHetero(PERCAPINC~PCTNOHINS+PCTMINOR+PCTUNIVDEG+PCTWHITE, hetero=~log(totPop),
#'                     data=validTractShp)
#' summary(mod.lmH)
#'
#' ## Use estimated weights as input for weighted lm-model.
#' ## This also to perform further model diagnostics.
#' mod.lmW <- lm(PERCAPINC~PCTNOHINS+PCTMINOR+PCTUNIVDEG+PCTWHITE, weights=mod.lmH$weights,
#'               data=validTractShp)
#' summary(mod.lmW)
#' hist(weighted.residuals(mod.lmW))
#'

lmHetero <- function(formula, hetero=~1, data, subset, na.action,
                     contrasts = NULL, iter=FALSE, ...) {
  ################################################################################################
  ##
  ## Purpose: Calculate multiplicately weighted regression models with respect to externally
  ##          given exogenous variables
  ##
  ## Source: The ML estimation procedure for multiplicately weighted regression is given
  ##         in Greene W. H. (2000). Econometric Analysis. 4th edition. Upper Saddle River:
  ##         Prentice Hall. pp 516-521 (Note: page numbers will differ for other editions)
  ##
  ## Objective: estimate gamma values for mulitiplicative heteroscedasticity models
  ##
  ## Syntax:
  ## [1] lmHetero(y~x1+x2, hetero= ~log(z1)+log(z2), data=mydata)
  ## [2] lmHetero(y~x1+x2 | log(z1)+log(z2), data=mydata)
  ## Note: An expression for Z must be present.
  ## y : Dependent variable
  ## X : Independent variable(s) with added intercept
  ## Z : Weight variable(s) with intercept added.
  ##
  ##     !!! Important user input: weighte variables Z must be enter log-transformed !!!
  ##
  ## Authors: Michael Tiefelsdorf (tiefelsdorf@utd.edu) & Yongwan Chun
  ###############################################################################################

  ## Parseing the function call
  #require(Formula)
  cl <- match.call()
  if (missing(data)) data <- environment(formula)
  mf <- match.call(expand.dots = FALSE)
  m <- match(c("formula", "data", "subset", "na.action"), names(mf), 0)
  mf <- mf[c(1, m)]
  mf$drop.unused.levels <- TRUE
  if (!missing(hetero)) {
    formula <- Formula::as.Formula(formula, hetero)
    cl$hetero <- NULL
    cl$formula <- formula(formula)
  }
  else {
    formula <- Formula::as.Formula(formula)
  }
  stopifnot(length(formula)[1] == 1L, length(formula)[2] %in% 1:2)
  has_dot <- function(formula) inherits(try(stats::terms(formula), silent = TRUE), "try-error")
  if (has_dot(formula)) {
    f1 <- formula(formula, rhs = 1)
    f2 <- formula(formula, lhs = 0, rhs = 2)
    if (!has_dot(f1) & has_dot(f2))
      formula <- Formula::as.Formula(f1, stats::update(formula(formula, lhs = 0, rhs = 1), f2))
  }
  mf$formula <- formula
  mf[[1]] <- as.name("model.frame")
  mf <- eval(mf, parent.frame())
  y <- stats::model.response(mf, "numeric")
  mt <- stats::terms(formula, data = data, ...)
  mtX <- stats::terms(formula, data = data, rhs = 1, ...)
  X <- stats::model.matrix(mtX, mf, contrasts)
  if (length(formula)[2] < 2L) {
    #stop("Set of weights variables is missing in function call!")
    #mtZ <- NULL
    #Z <- NULL
    zNames <- "(Intercept)"
    Z <- matrix(1, nrow=length(y), ncol=1)
  }
  else {
    mtZ <- stats::delete.response(stats::terms(formula, data = data, rhs = 2, ...))
    Z <- stats::model.matrix(mtZ, mf, contrasts)
    zNames <- colnames(Z)
    if (!all(exp(Z) > 0)) stop("All weight variables must be positive!")
    #Z[,-1] <- log(Z[,-1])
  }
  ## done::parsing and starting the estimation

  ##
  ## Calculate start values for Beta and Gamma
  ##
  nofreg <- length(y)
  Beta <- qr.solve(crossprod(X),crossprod(X,y))
  res <- y - X %*% Beta
  res <- log(res^2)

  Gamma <- qr.solve(crossprod(Z),crossprod(Z,res))
  Gamma[1] <- Gamma[1]+1.2704    # Harvey (1976) correction
  if (iter==TRUE) { cat("Beta:",Beta);cat("\n");cat("Gamma:",Gamma);cat("\n") }

  qrcz <- qr(crossprod(Z))

  ## Start interative estimation procedure
  MaxDiff <- Inf
  while (MaxDiff > 0.001) {
    Sigma2 <- exp(Z %*% Gamma)
    W <- diag(1 / Sigma2[,1])
    ## Calculate the values at the next step
    #  newBeta <- solve(t(x) %*% W %*% x, t(x) %*% W %*% y)
    xW <- crossprod(X,W)
    newBeta <- qr.solve(xW %*% X, xW %*% y)
    res2 <- (y - X %*% newBeta)^2
    newGamma <- Gamma + qr.coef(qrcz,crossprod(Z, res2 / Sigma2 - 1))
    MaxDiff <- max(abs(newGamma - Gamma))
    if (iter==TRUE) {cat("Beta:",newBeta);cat("\n");cat("Gamma:",newGamma);cat("\n")}
    Beta <- newBeta
    Gamma <- newGamma
  }  # end while

  Sigma2 <- exp(Z %*% newGamma)
  cBeta <- qr.solve(xW %*% X)      # The global covariance matrix is block diagonal
  cGamma <- 2*qr.solve(crossprod(Z))
  logLikeH1 <- -(nofreg/2)*log(2*pi) - 0.5*sum(log(Sigma2[,1])) - 0.5* sum(res2/Sigma2[,1])
  logLikeH0 <- stats::logLik(stats::lm(y~X))[1]
  ##
  ## Report results. As ML-estimates the values are slightly biased
  ##
  rval <- list()
  rval$CALL <- cl
  rval$sigma2 <- exp(newGamma[1])
  rval$gamma <- newGamma
  rval$namesGamma <- zNames
  rval$beta <- newBeta
  rval$weights <- 1/Sigma2[,1]
  rval$covBeta <- cBeta
  rval$covGamma <- cGamma
  rval$logLikeH1 <- logLikeH1
  rval$logLikeH0 <- logLikeH0
  #   rval$formula <- formula(formula)
  #   rval$terms <- list(regressors = mtX, hetero = mtZ, full = mt)
  #   rval$na.action <- attr(mf, "na.action")
  #   rval$levels <- .getXlevels(mt, mf)
  #   rval$contrasts <- list(regressors = attr(X, "contrasts"),
  #                        hetero = attr(Z, "contrasts"))
  class(rval) <- c("lmHetero","list")
  return(rval)
} # end::lmHetero

