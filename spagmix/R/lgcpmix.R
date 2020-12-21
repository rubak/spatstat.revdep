lgcpmix <- function(lambda,covmodel="exp",covpars=NULL){
  if(is.im(lambda)){
    if(any(lambda<0)) stop("'lambda' must be a non-negative real-valued spatial intensity function")
  } else if(inherits(lambda,"bivden")){
    lambda <- lambda$z
  } else {
    stop("'lambda' must be an object of class 'im' (spatstat) or 'bivden' (sparr)")
  }

  if(is.null(covpars)){
    warning("'covpars' is empty; setting var = 1 and scale = 1")
    covpars <- list(var=1,scale=1)
  }
  if(!is.list(covpars)) stop("'covpars' must be a named list; members matching the required contents for 'covmodel'")
  if(is.null(covpars$var)) covpars$var <- 1

  modgen <- spatstat::getRandomFieldsModelGen(covmodel)
  rfmodel <- do.call(modgen,covpars) + RandomFields::RMtrend(mean=-covpars$var/2)
  if(!inherits(rfmodel,"RMmodel")) stop("Problem generating RandomFields covariance model object",call.=FALSE)

  # offset <- 0.5*covpars$var
  lamvec <- log(as.vector(t(as.matrix(lambda))))
  # lamvec <- lamvec - offset
  lamvec[is.na(lamvec)] <- -Inf

  spc <- RandomFields::RFoptions()$general$spConform
  if(spc) RandomFields::RFoptions(spConform=FALSE)
  z <- RandomFields::RFsimulate(rfmodel,lambda$xcol,lambda$yrow,grid=TRUE)
  if(spc) RandomFields::RFoptions(spConform=TRUE)

  w <- as.mask(lambda)
  logLambda <- lamvec + z
  result <- matrix(t(exp(logLambda)),nrow=dim(lambda)[1],ncol=dim(lambda)[2])
  result <- as.im(result,W=w)[w,drop=FALSE]

  return(result)
}
