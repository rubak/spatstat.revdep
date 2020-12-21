
#####################################################################
# modificada a partir de "Kphylo.sp.envelope2"; Esta version es la que se llamaba  llamada "envelope.fv.idar4.R"

# toma como entrada el ppp original y como simulados solo los simulados de la especie focal (de simulador2: en data.frame)
# mippp: ppp multivariao
# mippp.sp.sim <- list of simulated focal patterns from simulador2.R
# mippp.sp <- observed univariate focal pp
# cross.idar: computing idar between different unities/age/size classes? IF TRUE, remove focal species from nmultivariate pp (mippp).
#      When computing "typical" isar, cross.idar=FALSE: the observed is replaced by the simulated
#      When computing cross-isar (btween different size-age-ckasses) the focal is not within the multivarite and 
#       there's nothin to replace and this is what happens when cross.idar=TRUE


# idar: name of the individual phylogenetic -area relationship to compute (isar, ipsear, ipsvar, ifdar, etc)
# tree: ape phylogenetic tree or variance/covarance matrix
# traits: data.frame of functional traits (species as columns, traits in rows)

# Es obligatorio suministar namesmark: mitable requiere un ppp multivariado con solo una columna de marcas y sin mnamesmark no se le puede suministrar

# fast envelope for isar analysis
# requires a list of simulated focal patterns


# En esta version 3, la matriz de distancias de gower se calcula SOLO UNA VEZ sobre toda la data.frame de traits
# Y se pasa a idar2 (en vez de traits) para ser subseteada segun sea necesario
# Eso hace que las distancias de gower entre especies sean siempre las mismas y no dependan de el tamano de la
# comunidad local (lo que ocurre ya que el Sijk usa el rango Rij para calcular las distancias individuales) 

envelope4idar <- function(mippp, mippp.sp.sim, mippp.sp, mimark=NULL, namesmark=NULL, r,
                                        idar="isar",
                                        buffer="adapt", bfw=NULL,nsim=NULL, nrank=1,
                                        tree = NULL,
                                        traits = NULL,
                                        cross.idar=FALSE,
                                        savefuns=TRUE,
                                       correct.phylo="exclude",
                                       correct.trait.na=FALSE,
                                       correct.trait="mean"
                                        ){


    if (!is.null(namesmark)) 
        mippp$marks <- factor((mippp$marks[namesmark][[1]]))
    if (!is.null(mimark)) 
        if (mimark %in% levels(mippp$marks) == FALSE) {
            stop(paste(mimark, " can't be recognized as a mark\n\n\n
                      have you indicated in which column of thedataframe are the species\n
                      marks? (argument 'namesmark'\n\n"))
        }
  
      # phylogenetic check and correction
      if(!is.null(tree)) tree<- checktree(tree=tree, mippp=mippp, idar=idar, correct.phylo=correct.phylo) ## NUEVA FUNCION de MANIPULACION DE TREES
      
      # Trait check and correction
       if (!is.null(traits)) traits <-checktraits(traits=traits, mippp=mippp, idar=idar, correct.trait.na=correct.trait.na, correct.trait=correct.trait)
       #if (!is.null(traits) & idar!="icwmar") traits <-checktraits(traits=traits, idar=idar, correct.trait.na=correct.trait.na, correct.trait=correct.trait)

   #*************************
   # BUFFER preprocessing
   # ************************* 
    #if the buffer is fixed, compute it only once and pass it to each simulation
    if (is.numeric(buffer) & is.null(bfw))  bfw <- owin(mippp$window$xrange + c(buffer, -buffer), 
                                                                        mippp$window$yrange + c(buffer, -buffer))

   #*************************
   # cross.idar preprocessing
   # *************************
    if (cross.idar == FALSE) {
        # we split the original ppp in focal (mippp.sp) and the rest (mippp)
        # in order to recompose  it with the simulated focals in idar2
        mippp.sp <- mippp[mippp$marks == mimark]
        # remove the focal species from the original mippp
        mippp <- mippp[mippp$marks != mimark]
    }

   #*************************
   # Computing observed and simulated idars
   # *************************
   # compute observed idar
    obs <- idar2(mippp.sp = mippp.sp, mippp = mippp, mimark = mimark, 
        buffer = buffer, bfw = bfw, r = r, cross.idar = cross.idar, 
        idar = idar, tree = tree, traits = traits)

    # compute simulated  idars   
    if (is.null(nsim)) nsim <- length(mippp.sp.sim)
    simvals = NULL
    cat("computing envelopes   ")
    for (i in 1:nsim) {
        progressreport(i, nsim)
        simvals <- cbind(simvals, idar2(mippp.sp = mippp.sp.sim[[i]], 
                          mippp = mippp, mimark = mimark, buffer = buffer, 
                          bfw = bfw, r = r, cross.idar = cross.idar, idar = idar, 
                         tree = tree, traits = traits))
    }
   #*************************
   # end of computing observed and simulated idars
   # *************************
   
    #*************************
    # Formating  and returning results
    #*************************
    # extract values for fv format
    if (nrank == 1) {
        lohi <- apply(simvals, 1, range)
    }
    else {
        lohi <- apply(simvals, 1, function(x, n) {
            sort(x)[n]
        }, n = c(nrank, nsim - nrank + 1))
    }
    lo <- lohi[1, ]
    hi <- lohi[2, ]
    lo.name <- paste("lower pointwise envelope of %s from simulations")
    hi.name <- paste("upper pointwise envelope of %s from simulations")
    m <- apply(simvals, 1, mean, na.rm = TRUE) # Beware that we ignore NA results !!!
      
    results <- data.frame(r = r, obs = obs, mmean = m, lo = lo, 
        hi = hi)

    if (idar == "isar") 
        ylab = substitute(ISAR(r), NULL)
    if (idar == "ipsrar") 
        ylab = substitute(IPSRAR(r), NULL)
    if (idar == "ipsvar") 
        ylab = substitute(IPSVAR(r), NULL)
    if (idar == "ipsear") 
        ylab = substitute(IPSEAR(r), NULL)
    if (idar == "ipscar") 
        ylab = substitute(IPSCAR(r), NULL)
    if (idar == "ifdar") 
        ylab = substitute(IFDAR(r), NULL)
    if (idar == "imntdar") 
        ylab = substitute(IMNTDAR(r), NULL)
    if (idar == "icwmar") 
        ylab = substitute(ICWMAR(r), NULL)
   if (idar == "icwmar.O") 
        ylab = substitute(ICWMAR.O(r), NULL)	
    if (idar == "iraodar")                                                
        ylab = substitute(IRaoDAR(r), NULL) 
 if (idar == "iraodar.O")                                                
        ylab = substitute(IRaoDAR.O(r), NULL) 

    result <- fv(results, argu = "r", ylab = ylab, valu = "obs", 
        fmla = . ~ r, alim = c(min(r), max(r)), labl = c("r", 
            "%s[obs](r)", "bar(%s)(r)", "%s[lo](r)", "%s[hi](r)"), 
        desc = c("distance argument r", "observed value of %s for data pattern", 
            "sample mean of %s from simulations", "lower pointwise envelope of %s from simulations", 
            "upper pointwise envelope of %s from simulations"), 
        fname = idar)
    class(result) <- c("envelope", "fv", "data.frame")
    attr(result, "dotnames") <- c("obs", "mmean", "hi", "lo")
    attr(result, "einfo") <- list(global = FALSE, nrank = nrank, 
        nsim = nsim, Nsim = nsim, dual = FALSE, nsim2 = nsim, 
        valname = "isar", VARIANCE = FALSE, csr = FALSE)
    if (savefuns == TRUE) 
        attr(result, "simfuns") <- data.frame(r = r, simvals)
    return(result)
}

#-------------------------------------------------------------------



####################################################################################################
####################################################################################################
## idar2
## derivada de  "Kphylo.sp2", para emplear en las envelopes. La entrada es mas flexible

idar2 <- function(mippp.sp, mippp, mimark, idar="isar", buffer,bfw, r, cross.idar=FALSE,  tree = NULL, traits = NULL) {
        
       #mippp.sp puede ser el univariado observado o viene de simulador2. contiene solo la especie focal (observada o simulada)
       # mippp tiene todas las otras especies excepto la focal.
       # viene de un paso externo (solo tiene que hacerse una vez)
       #    mippp <- mippp[mippp$marks!=mimark]
       # 
       # Si el bufer es fijo, tampoco hay que calcularlo cada vez: se suministra de un paso anterior:
       #    bfw<- owin( mippp$window$xrange + c(buffer,-buffer), mippp$window$yrange + c(buffer,-buffer))

    if (cross.idar == FALSE) {
       # We paste the focal ppp (observed or simulated) to the multivariate ppp    
        mippp$x <- c(mippp$x, mippp.sp$x)
        mippp$y <- c(mippp$y, mippp.sp$y)
        mippp$marks <- (c(as.character(mippp$marks), as.character(mippp.sp$marks)))
    }

    # If buffer is fixed, we trim now the focal pattern     
    if (!is.null(bfw)) 
        mippp.sp <- mippp.sp[inside.owin(mippp.sp, w = bfw)]

   # get local comunities around each focal tree for all  r's
    # result will be a list with length=  length(r) and each element in the list a data.frame with dimensions [mippp.sp$n, mippp$n]

    ######TODO: allow selection for functions mitable2, mitable3 etc  with sum of point sizes etc (from prototipomitable2m_3.R)
    cosamt <- mitable(mippp.sp, mippp, r) 
    
    # compute abundance in "rings" (for cwm.O computation)
    # BEWARE: this must be made before buffer correction
    # TODO: implement in frotran ???
    if( idar %in% c("iraodar.O", "icwmar.O")){
       cosamt.O <- cosamt
       for ( i in 2:length(cosamt.O)) cosamt.O[[i]] <- cosamt[[i]]-cosamt[[i-1]]
       # apply buffer correction if indicated
       if (buffer == "adapt") {
            bdp <- bdist.points(mippp.sp)
            for (i in 1:length(r))cosamt.O[[i]] <- cosamt.O[[i]][bdp >= r[i], ]
        }
    }
    
    
    #If indicated, correct edge effect with an adaptive buffer
    if(!idar %in% c("iraodar.O", "icwmar.O")){
       if (buffer == "adapt") {
           bdp <- bdist.points(mippp.sp)
          for (i in 1:length(r)) cosamt[[i]] <- cosamt[[i]][bdp >= r[i], ]
       }
    }

    #Define function to compute ISAR --------------------------------------------------------------------------------
    # TODO: this could be an external function, as psv, pse, FD, etc
    isar.r <- function(x) {
        # bivariate emptiness probability
        Ptj <- function(x) sum(x == 0)/length(x)
        #compute isar    
        result <- sum(apply(x, 2, function(x) 1 - Ptj(x)))
        return(result)
    }
    # ---- end function isar.r -------------------------------------------------------------------------------------------------------

   # Definition of function cwm.r
   # OJO: traits here is just a vector with the values of only one trait
   cwm.r<- function(x, traits){
            if(!is.null(dim(traits))) stop ("to compute icwmar, 'traits' mut be a vector with the values of juts one trait")
            return(sum(x*traits, na.rm=TRUE)/sum(x))
   }

    # Choose analysis
    if (idar == "icwmar") {
       traits.ok<- traits[names(traits)%in%colnames(cosamt[[1]])] # added 20/04/2016
       micwmar <-sapply(cosamt, function(tr) {
                         if (is.null(dim(tr)))   tr <- rbind(tr, tr) # If there is only one point, repeat its local assemblage to build a matrix and still be able to use "apply"
                         cwm.ind<-apply(tr, 1, function(x) cwm.r(x, traits=traits.ok))
                         cwm.ind[is.na(cwm.ind)] <- 0 # Set 0 fot those points without neighbours (i.e., no traits in their surroundings)
                        return(mean(cwm.ind))
            })
       result<-micwmar 
    }
    
     if (idar == "icwmar.O") {
           traits.ok<- traits[names(traits)%in%colnames(cosamt.O[[1]])] # added 20/04/2016
           micwmar <-sapply(cosamt.O, function(tr) { # use cosamt.O instead of cosamt
	                # If there is only one point, repeat its local assemblage 
	                #to build a matrix and still be able to use "apply"
                         if (is.null(dim(tr)))   tr <- rbind(tr, tr) 
                         cwm.ind<-apply(tr, 1, function(x) cwm.r(x, traits=traits.ok))
			 # Set 0 fot those points without neighbours (i.e., no traits in their surroundings)
                         cwm.ind[is.na(cwm.ind)] <- 0 
                        return(mean(cwm.ind))
            })
       result<-micwmar 
    }
    
    
    if (idar == "iraodar") {
       miraodar <-sapply(cosamt, function(x) {
             if (is.null(dim(x)))   x <- rbind(x, x)
	     if(dim(x)[1]<1) return(NA)
             res<- raoDmod(comm=x, phy=tree)
	     return(res)
           })
          result <- miraodar
    }
    
    if (idar == "iraodar.O") {
       miraodar <-sapply(cosamt.O, function(x) {
             if (is.null(dim(x)))   x <- rbind(x, x)
	     if(dim(x)[1]<1) return(NA)
             res<- raoDmod(comm=x, phy=tree)
	     return(res)
           })
          result <- miraodar
    }
       

    if (idar == "ifdar") {
       #Compute IFDAR(r)  (i.e., IFdisp) for each radius r
        mifdar <- sapply(cosamt, function(x) {
                       if (is.null(dim(x)))   x <- rbind(x, x)
                       res <- fdis(x, traits = traits) 
                       return(res)
        })
        result <- mifdar
    }
    
   if (idar == "isar") {
        #Compute ISAR(r) for each radius r
        misar <- sapply(cosamt, function(x) {
                      if (is.null(dim(x)))  x <- rbind(x, x)
                      res <- isar.r(x)
                      return(res)
        })
        result <- misar
    }

    if (idar == "ipsear") {
       # first for each individual
        mipsear <- lapply(cosamt, function(x) {
            if (is.null(dim(x)))  x <- rbind(x, x)
            if (sum(colSums(x) > 0) > 1)   res <- pse(x, tree = tree)$PSE
            if (sum(colSums(x) > 0) <= 1)  res <- rep(NA, dim(x)[1]) # repit NA for all individuals of he focal species
            return(res)
        })
        # then average ipse for each radius
        result <- sapply(mipsear, mean, na.rm = T)# TODO: optional na.rm. At short distances lot of NA's and maybe unreliable means
    }

    if (idar == "ipsvar") {
        # first for each individual
        mipsvar <- lapply(cosamt, function(x) {
            if (is.null(dim(x)))  x <- rbind(x, x)
            if (sum(colSums(x) > 0) > 1)   res <- psv(x, tree = tree, compute.var = F)$PSVs
            if (sum(colSums(x) > 0) <= 1) res <- rep(NA, dim(x)[1]) # repit NA for all individuals of he focal species
            return(res)
        })
        # then average ipsv for each radius
        result <- sapply(mipsvar, mean, na.rm = T)# TODO: optional na.rm. At short distances lot of NA's and maybe unreliable means
    }

    if (idar == "ipsrar") {
         # first for each individual
        mipsrar <- lapply(cosamt, function(x) {
            if (is.null(dim(x)))   x <- rbind(x, x)
            if (sum(colSums(x) > 0) > 1)   res <- psr(x, tree = tree, compute.var = F)$PSR
            if (sum(colSums(x) > 0) <= 1)  res <- rep(NA, dim(x)[1]) # repit NA for all individuals of he focal species
            return(res)
        })
        # then average ipsr for each radius
        result <- sapply(mipsrar, mean, na.rm = T)# TODO: optional na.rm. At short distances lot of NA's and maybe unreliable means
    }

    if (idar == "ipscar") {
        # first for each individual
        mipscar <- lapply(cosamt, function(x) {
            if (is.null(dim(x))) 
                x <- rbind(x, x)
            if (sum(colSums(x) > 0) > 1) 
                res <- psc(x, tree = tree)$PSCs
            if (sum(colSums(x) > 0) <= 1) 
                res <- rep(NA, dim(x)[1])
            return(res)
        })
        # then average ipsc for each radius
        result <- sapply(mipscar, mean, na.rm = T)# TODO: optional na.rm. At short distances lot of NA's and maybe unreliable means
    }

    if (idar == "imntdar") {
        # first for each individual
        mimntdar <- lapply(cosamt, function(x) {
            if (is.null(dim(x))) 
                x <- rbind(x, x)
            if (sum(colSums(x) > 0) > 1) 
                res <- mntd(x, dis = tree, abundance.weighted = FALSE)
            if (sum(colSums(x) > 0) <= 1) 
                res <- rep(NA, dim(x)[1])
            return(res)
        })
        # then average imntdar for each radius
        result <- sapply(mimntdar, mean, na.rm = T) # TODO: optional na.rm. At short distances lot of NA's and maybe unreliable means
    }
    return(result)
}


######################################################
# Funcion simulador2
# para generar heterogeneos poisson para calcular las envueltas
#
# la salida no es un ppp sino un dataframe con la especie
# simulada. Esto evita tener que escribir repetidamente 
# cientos de miles de datos que no se simulan
# (se simulan como mucho cientos de individuos de la especie focal)
######################################################

simulador2 <- function(mimark, milambda, nsim=99){
    cosapp <- vector("list", nsim)
    cosapp <- lapply(cosapp, function(x) {
                      x<- rpoispp(milambda)
                      x$marks<-factor(rep(mimark, x$n))
                      return(x)
                      }
              )
return(cosapp)
}
#-------------------------------------------------------------------





