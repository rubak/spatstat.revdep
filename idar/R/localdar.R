
###################################################################################
# funcion para dibujar (y sacar datos) de la distribucion de riqueza a una escala determinada en una rejilla
# mippp: patron marcado de puntos
      
localdar<- 
function (mippp, mippp.sp=NULL, nx=NULL, ny=NULL, mimark=NULL, idar = "isar", buffer=0, bfw=NULL, 
    r, cross.idar = FALSE, tree = NULL, traits = NULL, namesmark=NULL, correct.trait.na=TRUE, correct.trait = "mean", correct.phylo="mean") 
{
     gridok <- FALSE
     bufferex<-FALSE
     bufferect<- FALSE

    # 1) Check that the column of species is clearly set
    
    if(!is.marked(mippp)) stop ("mapIDAR requires a marked point pattern")
    if (!is.null(namesmark))  mippp$marks <- factor((mippp$marks[namesmark][[1]]))
    
    dmark<-dim(marks(mippp))
      
      if(is.null(namesmark) & !is.null(dmark)){
          if(dmark[2]==1) marks(mippp)<-mippp$marks[,1] else stop(
	  "you should indicate which column of the dataframe of marks
	  stores species names (argument 'namesmark')\n\n") 
      }
      # from here on, mippp has only a vector of marks.




   # 2) Check  th existence of a focal pattern, indicated by "mimark" or "mippp.sp" 
      
   #  if (!is.null(mimark))  if (mimark %in% levels(mippp$marks) == FALSE) {
   #         stop(paste(mimark, " can't be recognized as a mark\n\n\n
   #	    have you indicated in which column of thedataframe are the species\n
   #	    marks? (argument 'namesmark'\n\n"))
   #     }

    # if some mark is given in mimark, we will estimate the idar values around each individual
    
    if(!is.null(mippp.sp))  gridok <- FALSE # it is not a "inner" grid
    if (!is.null(mimark)){
 
           if (mimark %in% levels(mippp$marks) == FALSE) {
               stop(paste(mimark, " can't be recognized as a mark\n\n\n
	       have you indicated in which column of thedataframe are the species\n
	       marks? (argument 'namesmark'\n\n"))
             }
       mippp.sp<- unmark(mippp[mippp$marks==mimark])
       gridok<- FALSE
       
     }

   if(is.null(mippp.sp)){
        if(is.null(nx)& is.null(ny)) nx<-ny<-30 #por defecto 30 x 30 pixels
	if(is.null(nx)& !is.null(ny)) nx<-ny
	if(is.null(nx)& !is.null(ny)) ny<-nx
	
        gridxy <-gridcentres(mippp$window, nx, ny)
	# select only points inside the window
	okbig<-inside.owin(gridxy, w=mippp$window)
        mippp.sp <- ppp(x=gridxy$x[okbig], y=gridxy$y[okbig], window=mippp$window)
	# if we use a regular grid, by default cross.idar is TRUE
      #  cross.idar =TRUE  ##!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	gridok<- TRUE
 #print(mippp.sp)  
   }

    # BUFFER MANAGEMENT#--------------------------------------------------------------------------------------------------

     if(!is.null(bfw)){
           bufferex<- TRUE
           ok<- inside.owin(mippp.sp, w = bfw)
	   mippp.sp0 <- mippp.sp # keep the original focal pattern to return resuilts
	   mippp.sp <- mippp.sp[ok]

     }

   if(buffer!=0 & buffer !="adapt" & is.null(bfw)){ 
         if (is.numeric(buffer) & is.null(bfw)){ 
              if(mippp$window$type!="rectangle") stop("numeric buffer only available for rectangular windows")
              bfw <- owin(mippp$window$xrange + c(buffer, -buffer), 
                      mippp$window$yrange + c(buffer, -buffer))
	bufferect<- TRUE
#print(mippp.sp)	   
	}	      
	if (!gridok){
	   ok<- inside.owin(mippp.sp, w = bfw)
	   mippp.sp0 <- mippp.sp # keep the original focal pattern to return resuilts
	   mippp.sp <- mippp.sp[ok]

	}
	if (gridok){
	   #ok<-inside.owin(gridxy, w=bfw)
	   #mippp.sp <- ppp(x=gridxy$x[ok], y=gridxy$y[ok], window=mippp$window)
	   okbig<-inside.owin(gridxy, w=bfw)
	   mippp.sp <- ppp(x=gridxy$x[okbig], y=gridxy$y[okbig], window=mippp$window)
	}
   }
    
    
        #si el idar es cruzado, se elimina el ppp focal del global
    if (cross.idar == TRUE & gridok ==FALSE) {
        if(is.null(mimark)) stop("for crossed maps you should indicate a focal species (argument 'mimark')")
        mippp<- mippp[mippp$marks!=mimark]
    }

       
	
     if (!is.null(tree)) 
        tree <- checktree(tree = tree, mippp = mippp, idar = idar, 
            correct.phylo = correct.phylo)
    
    if (!is.null(traits)) 
        traits <- checktraits(traits = traits, mippp = mippp, 
            idar = idar, correct.trait.na = correct.trait.na, 
            correct.trait = correct.trait)
	    
    

	
	#*********************************************************************************
	 # calculamos las distancias entre los puntos de la grid y llos arboles 
         # el resultado sera una lista de longitud = length(r) y en cada entrada una tabla de dimension mippp.sp$n X mippp$n
        # con el numero de individuos de cada especie a la distancia r de cada individuo de la especie focal 
		
        cosamt <- mitable(mippp.sp, mippp, r)
    
         #*********************************************************************************
    
    
    
    if (idar %in% c("iraodar.O", "icwmar.O")) {
        cosamt.O <- cosamt
        for (i in 2:length(cosamt.O)) cosamt.O[[i]] <- cosamt[[i]] - 
            cosamt[[i - 1]]
        if (buffer == "adapt") {
	    ok<- NULL
            bdp <- bdist.points(mippp.sp)
            for (i in 1:length(r)){
		ok<- cbind(ok, bdp >= r[i])
	        cosamt.O[[i]] <- cosamt.O[[i]][bdp >=  r[i], ]
            }
        }
    }
    if (!idar %in% c("iraodar.O", "icwmar.O")) {
        if (buffer == "adapt") {
	    ok<-NULL
            bdp <- bdist.points(mippp.sp)
            for (i in 1:length(r)){
	         ok<- cbind(ok, bdp >= r[i])
	         cosamt[[i]] <- cosamt[[i]][bdp >= r[i], ]
	     }
        }
    }
    
    
    isar.r <- function(x) {
        Ptj <- function(x) sum(x == 0)/length(x)
        result <- sum(apply(x, 2, function(x) 1 - Ptj(x)))
        return(result)
    }
    cwm.r <- function(x, traits) {
        if (!is.null(dim(traits))) 
            stop("to compute icwmar, 'traits' mut be a vector with the values of juts one trait")
	    
        return(sum(x * traits, na.rm = TRUE)/sum(x))
    }
    
    
    
    
    if (idar == "icwmar") {
    
        micwmar <- sapply(cosamt, function(tr) {
            if (is.null(dim(tr))) 
                tr <- rbind(tr, tr)
      
         #traits<-traits[match( dimnames(tr)[[2]], names(traits))] #!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        #cwm.ind <- apply(tr, 1, function(x) cwm.r(x, traits = traits))
	cwm.ind <- apply(tr, 1, function(x) cwm.r(x, traits = traits[match( dimnames(tr)[[2]], names(traits))])) #!!!!!!!!!!!!!!!!!!!!!!!!!OJO EN FUNCIONES INDIVIDUALES: meyter mejor en checktraits
            cwm.ind[is.na(cwm.ind)] <- 0   # OJO: esto tiene que ser 0 o NA ??? !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            #return(mean(cwm.ind))
	    return(cwm.ind)
        })
        result <- micwmar
    }
    if (idar == "icwmar.O") {
        micwmar <- sapply(cosamt.O, function(tr) {
            if (is.null(dim(tr))) 
                tr <- rbind(tr, tr)
            #cwm.ind <- apply(tr, 1, function(x) cwm.r(x, traits = traits))
	    cwm.ind <- apply(tr, 1, function(x) cwm.r(x, traits = traits[match( dimnames(tr)[[2]], names(traits))])) #!!!!!!!!!!!!!!!!!!!!!!!!!OJO EN FUNCIONES INDIVIDUALES: meyter mejor en checktraits
            cwm.ind[is.na(cwm.ind)] <- 0  # OJO: esto tiene que ser 0 o NA ???  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1
            #return(mean(cwm.ind))
	    return(cwm.ind)
        })
        result <- micwmar
    }
    if (idar == "iraodar") {
        miraodar <- sapply(cosamt, function(x) {
            if (is.null(dim(x))) 
                x <- rbind(x, x)
            if (dim(x)[1] < 1) 
                return(NA)
            res <- raoDmap(comm = x, phy = tree)
            return(res)
        })
        result <- miraodar
    }
    if (idar == "iraodar.O") {
        miraodar <- sapply(cosamt.O, function(x) {
            if (is.null(dim(x))) 
                x <- rbind(x, x)
            if (dim(x)[1] < 1) 
                return(NA)
            res <- raoDmap(comm = x, phy = tree)
            return(res)
        })
        result <- miraodar
    }
    if (idar == "ifdar") {
        mifdar <- sapply(cosamt, function(x) {
            if (is.null(dim(x))) 
                x <- rbind(x, x)
            #res <- fdis(x, traits = traits)
	    res <- fdismap(x, traits = traits)
	    
            return(res)
        })
        result <- mifdar
    }
    if (idar == "isar") {
        misar <- sapply(cosamt, function(x) {
            if (is.null(dim(x))) 
                x <- rbind(x, x)
            #res <- isar.r(x)
            res<-apply(x,1,function(y) sum(y>0))
	    
	    return(res)
        })
        result <- misar
    }
    if (idar == "ipsear") {
        mipsear <- lapply(cosamt, function(x) {
            if (is.null(dim(x))) 
                x <- rbind(x, x)
            if (sum(colSums(x) > 0) > 1) 
                res <- pse(x, tree = tree)$PSE
            if (sum(colSums(x) > 0) <= 1) 
                res <- rep(NA, dim(x)[1])
            return(res)
        })
        #result <- sapply(mipsear, mean, na.rm = T)
	result <- mipsear
    }
    if (idar == "ipsvar") {
        mipsvar <- lapply(cosamt, function(x) {
            if (is.null(dim(x))) 
                x <- rbind(x, x)
            if (sum(colSums(x) > 0) > 1) 
                res <- psv(x, tree = tree, compute.var = F)$PSVs
            if (sum(colSums(x) > 0) <= 1) 
                res <- rep(NA, dim(x)[1])
            return(res)
        })
        #result <- sapply(mipsvar, mean, na.rm = T)
	result <-mipsvar
    }
    if (idar == "ipsrar") {
        mipsrar <- lapply(cosamt, function(x) {
            if (is.null(dim(x))) 
                x <- rbind(x, x)
            if (sum(colSums(x) > 0) > 1) 
                res <- psr(x, tree = tree, compute.var = F)$PSR
            if (sum(colSums(x) > 0) <= 1) 
                res <- rep(NA, dim(x)[1])
            return(res)
        })
        #result <- sapply(mipsrar, mean, na.rm = T)
	result <-mipsrar
    }
    if (idar == "ipscar") {
        mipscar <- lapply(cosamt, function(x) {
            if (is.null(dim(x))) 
                x <- rbind(x, x)
            if (sum(colSums(x) > 0) > 1) 
                res <- psc(x, tree = tree)$PSCs
            if (sum(colSums(x) > 0) <= 1) 
                res <- rep(NA, dim(x)[1])
            return(res)
        })
        #result <- sapply(mipscar, mean, na.rm = T)
	result <-mipscar
    }
    if (idar == "imntdar") {
        mimntdar <- lapply(cosamt, function(x) {
            if (is.null(dim(x))) 
                x <- rbind(x, x)
            if (sum(colSums(x) > 0) > 1) 
                res <- mntd(x, dis = tree, abundance.weighted = FALSE)
            if (sum(colSums(x) > 0) <= 1) 
                res <- rep(NA, dim(x)[1])
            return(res)
        })
        #result <- sapply(mimntdar, mean, na.rm = T)
	result <-mimntdar
    }
    
    
    # reconversion a data.frame de nlos resultados tipo matriz, para que puedan ser 
    if(is.matrix(result)) result<- as.data.frame(result)
    
    
  #print("hasta aqui hemos llegado")  
    	# si se ha suministrado o se espera una rejilla #######################################################################33
	if(gridok) {
#cat("okbig=", sum(okbig), "\n")
#cat("ok=", sum(ok), "\n")
	   # if(!is.null(bfw))ok<- inside.owin(gridxy, w=bfw) #!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! OJOJOJOJOJOJOJOJOJOJ
	     
	    result.im=list()
	    for (i in 1: length(cosamt)){
	        #mapaidar.ppp<- ppp(x=gridxy$x[ok], y=gridxy$y[ok], window=mippp$window, marks=result)
	        resultgrid=rep(NA,nx*ny)
		# resultgrid[ok] <-result[[i]] 
		#if(buffer!="adapt") resultgrid[ok] <-result[[i]] else resultgrid[ok[,i]] <-result[[i]] 
		#if(buffer!="adapt") resultgrid[okbig][ok] <-result[[i]] else resultgrid[okbig][ok[,i]] <-result[[i]] 
		if(buffer=="adapt") resultgrid[okbig][ok[,i]] <-result[[i]] 
		if(buffer!="adapt" & buffer!=0 & !bufferect) resultgrid[okbig][ok] <-result[[i]] 
		if(bufferect) resultgrid[okbig] <-result[[i]]
		if(bufferex) resultgrid[okbig][ok] <-result[[i]] 
		if(buffer==0 & !bufferex) resultgrid[okbig]<-result[[i]] 
	        result.im[[i]]<- as.im(list(x=unique(gridxy$x),y=unique(gridxy$y),z=matrix(resultgrid,nx,ny)))
	       }
	     #  names(result.im)<- r
	}
	
	# si sehasuministrado o se esperan resultados para un patron focal
	if(!gridok) {
        if(bufferex | bufferect) mippp.sp<- mippp.sp0 # the original non-trimed focal pattern
	result.im=list()
#print(ok)	
#print(result[[1]])
#print(result[[2]])
#print(mippp.sp)
#print(buffer)
#print(bufferex)
        for (i in 1: length(cosamt)){
		marks(mippp.sp)<-NA
		#if(buffer!="adapt")    marks(mippp.sp) [ok]<- result[[i]] else  marks(mippp.sp)[ok[,i]] <-result[[i]] 
		if(buffer=="adapt") marks(mippp.sp)[ok[,i]] <-result[[i]] 
		 if(buffer!="adapt" & buffer!=0 )   marks(mippp.sp) [ok]<- result[[i]]   
		 if(buffer==0 & !bufferex)   marks(mippp.sp) <- result[[i]]
		 if(bufferex) marks(mippp.sp)[ok]<- result[[i]]   
	        #marks(mippp.sp)<- result[[i]]
	        result.im[[i]] <- mippp.sp    
            }
	   #names(result.im)<- r
	 

	}

        names(result.im)<- paste(idar, "_",r* mippp$window$units$multiplier,"_", mippp$window$units$plural, sep="")
	return(result.im)
		
}




#################################################################################
 #funcion necesaria para mapear los componentes individuales de raoD (sin tener en cuenta la ponderacion; podria hacerse opara que devolviera algun valor ponderado por sampl.relabund)
 raoDmap<-
function (comm, phy = NULL) 
{
    res <- list()
    if (is.null(phy)) {
        tij <- 1 - diag(x = rep(1, length(comm[1, ])))
    }
    else {
        #if (class(phy) != "phylo") {  # !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! CAMBIADO 04/12/2019 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	if (!inherits (phy, what= "phylo")) {
            #if (!is.matrix(phy) & class(phy) != "dist")  # !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! CAMBIADO 04/12/2019 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	     if (!is.matrix(phy) & !inherits (phy, what=  "dist") )
                stop("Phy must be a distance matrix")
            if (is.matrix(phy)) 
                phy <- as.dist(phy)
            dat <- match.comm.dist(comm, phy)
            comm <- dat$comm
            phy <- dat$dist
            tij <- as.matrix(phy/2)
        }
         #if (class(phy) == "phylo") {  # !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! CAMBIADO 04/12/2019 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	if (inherits (phy, what= "phylo")) {
            if (!is.ultrametric(phy)) 
                stop("Phylogeny must be ultrametric")
            dat <- match.phylo.comm(phy, comm)
            comm <- dat$comm
            phy <- dat$phy
            tij <- cophenetic(phy)/2
        }
    }
    x <- as.matrix(comm)
    S <- length(x[1, ])
    N <- length(x[, 1])
    total <- apply(x, 1, sum)
    samp.relabund <- total/sum(x)
    x.combined <- matrix(apply(x, 2, sum), nrow = 1)/sum(x)
    x <- sweep(x, 1, total, "/")
    D <- vector(length = N)
    names(D) <- rownames(x)
    for (k in 1:N) D[k] <- sum(tij * outer(as.vector(t(x[k, ])), 
        as.vector(t(x[k, ]))))
    res$Dkk <- D
    res$alpha <- sum(res$Dkk * samp.relabund)
    return(res$Dkk)
}





################################################################################33333



# fdis para los mapas: no devuelve la media sino los valores en cada asbol/punto
fdismap<-
function (comm, traits) 
{
    sp.a <- colSums(comm) > 0
    filas <- dim(comm)[1]
    
    if (sum(sp.a) < 1) 
        return(NA)
    if (sum(sp.a) == 1) 
        return(0)
    if (sum(sp.a) > 1) {
        m.ok <- comm[, sp.a]
        sp.trait.ok <- !is.na(match(rownames(traits), colnames(m.ok)))
        if (!is.null(dim(m.ok))) {
            com.G.0 <- rowSums(m.ok) > 0
            cosad <- as.dist(traits[sp.trait.ok, sp.trait.ok, 
                drop = FALSE])
            
	    if (sum(cosad, na.rm = TRUE) == 0) 
                return(0)
            
	    if (sum(sp.trait.ok) >= 2) {
                com.G.02 <- rowSums(m.ok[, labels(cosad)]) > 
                  0
                cosaf <- fdisp(d = cosad, a = m.ok[com.G.02, 
                  labels(cosad), drop = FALSE], tol = 1e-07)$FDis
		
		result<-rep(NA, filas)
		result[com.G.02]<- cosaf
            }
            
	    if (sum(sp.trait.ok) < 2) {
                cosaf <- fdisp(d = cosad, a = t(as.matrix(m.ok[com.G.0, 
                  labels(cosad)])), tol = 1e-07)$FDis
                result<-rep(NA, filas)
		result[com.G.0]<- cosaf
	    }
            
        }
        if (is.null(dim(m.ok))) 
            result <- NA
        return(result)
    }
}


