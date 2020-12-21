  #*************************
     # Phylogenetic manipulations  TODO: phylogenetic manipulations as an external function
     #************************* 
    # tree : a phylo object (ape) or a phylogenetic covariance matrix

checktree<- function(tree,  mippp, idar, correct.phylo){
   
        treeold <- tree
        #if (class(tree) == "phylo")  tree <- vcv.phylo(tree, corr = TRUE) #!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! CHANGED 04/12/2019!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	if (inherits(tree, what= "phylo")) tree <- vcv.phylo(tree, corr = TRUE)
        sp.nok <- which(is.na(match(levels(mippp$marks), dimnames(tree)[[1]])))
        len.sp.nok <- length(sp.nok)
        names.sp.nok <- levels(mippp$marks)[sp.nok]
        # if there is some species missing in the phylogenetic tree
        if (len.sp.nok > 0) {
            # option 1) include missing species in the tree with a constant mean phylogenetic covariance
            if (correct.phylo == "mean") {
                meanphylo <- mean(as.dist(tree))
                dimtree <- dim(tree)
                newtree <- matrix(1, nrow = dimtree[1] + len.sp.nok, 
                  ncol = dimtree[1] + len.sp.nok)
                newtree[1:len.sp.nok, ] <- meanphylo
                newtree[, 1:len.sp.nok] <- meanphylo
                newtree[-(1:len.sp.nok), -(1:len.sp.nok)] <- tree
                diag(newtree) <- 1
                dimnames(newtree) <- lapply(dimnames(tree), function(x) c(names.sp.nok, 
                  x))
                warning(paste("species", names.sp.nok, "was not in the phylotree and has been included with constant mean phylogenetc covariance\n"))
                tree <- newtree
            }
            # option 2) eliminate missing species from the ppp
            if (correct.phylo == "exclude") {
                numbers.nok <- table(mippp[mippp$marks %in% names.sp.nok]$marks)[names.sp.nok]
                mippp <- mippp[!mippp$marks %in% names.sp.nok]
                mippp$marks <- factor(mippp$marks)
                warning(paste("species ", names.sp.nok, ", with ", 
                  numbers.nok, " individuals, was not in the phylotree and has been excluded from the analysis\n", 
                  sep = ""))
            }
        }
        if (idar == "imntdar") { # OJO: poner al principio de la funcion, para que corrija las especies faltantes (pasar cophenetic a "matrix" y al final de todas las correcciones devolver otr vez como "as.dist"
            #if (class(treeold) != "phylo")  #!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! CHANGED 04/12/2019!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	    if (!inherits(treeold, what= "phylo") )   stop("Computing IMNTDAR requires a phylogenetic tree of class 'phylo'\n\n")
            tree <- cophenetic(treeold)
        }
    return(tree)
    }
         

      #*************************
      # Trait manipulations. TODO: trait manipulations as an external function
      # ************************* 

checktraits<- function(traits, mippp, idar, correct.trait.na, correct.trait){

        # If traits come as a vector , do this:
        #if (is.null(dim(traits)) & class(traits) != "dist") {   #!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! CHANGED 04/12/2019!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	if (is.null(dim(traits)) & !inherits(traits, what= "dist")) {
            # first, check that the vector of traits has names simuilar to species names in ppp$marks 
            if (sum(!is.na(match(levels(mippp$marks), names(traits)[[1]]))) < 1)   stop(
            "the vector of traits should be named and should wear the same labels as the names\n of species in ppp$marks")

            # then deal with NA's in the vector of traits
            meantrait <- mean(traits, na.rm = T)
            if (correct.trait.na == TRUE) {
                traits[is.na(traits)] <- meantrait
            }

            # check that all species in the community are in the vector of traits
            sp.nok <- which(is.na(match(levels(mippp$marks), 
                                    names(traits))))
            len.sp.nok <- length(sp.nok)
            names.sp.nok <- levels(mippp$marks)[sp.nok]

            # if there is some species missing in the vector of traits...
            if (len.sp.nok > 0) {

                # option 1) include missing species in the vector with a mean trait value
                if (correct.trait == "mean") {
                  # repeat mean traits as many times as there are missing species
                  traitplus <- NULL
                  for (i in 1:len.sp.nok) traitplus <- c(traitplus, meantrait)
                  names(traitplus) <- names.sp.nok
                  traits <- c(traits, traitplus)
                  traits <- traits[order(names(traits))]
                  warning(paste("species", names.sp.nok,
                  "was not in the vector of traits and\n has been included with mean trait value\n"))
                }

                # option 2) eliminate missing species from the ppp
                if (correct.trait == "exclude") {
                  numbers.nok <- table(mippp[mippp$marks %in% names.sp.nok]$marks)[names.sp.nok]
                  mippp <- mippp[!mippp$marks %in% names.sp.nok]
                  mippp$marks <- factor(mippp$marks)
                  warning(paste("species ", names.sp.nok, ", with ", numbers.nok,
                  " individuals, was not in the data.frame of traits and has\n been excluded from the analysis\n",
                    sep = ""))
                }
            }
            #traits se sustituye por la matriz de distancias de traits solo si vamos a calcular el IFDAR
            if(idar=="ifdar") traits <- as.matrix(gowdis(traits)) 
        }

        # If traits come as a data.frame, do this:
        if (!is.null(dim(traits))) {# for traits as a data.frame
            # first, check that the data.frame of traits has dimnames simuilar to species names in ppp$marks 
            if (sum(!is.na(match(levels(mippp$marks), dimnames(traits)[[1]]))) < 1)  stop(
            "the rownames of the data.frame of traits should  wear the same labels as the names\n of species in ppp$marks")
            # then deal with NA's in the table of traits
            meantrait <- apply(traits, 2, mean, na.rm = T)
            if (correct.trait.na == TRUE) {
                for (i in 1:dim(traits)[2]) traits[is.na(traits[, i]), i] <- meantrait[i]
            }

            # check that all species in the community are in the table of traits
            sp.nok <- which(is.na(match(levels(mippp$marks), dimnames(traits)[[1]])))
            len.sp.nok <- length(sp.nok)
            names.sp.nok <- levels(mippp$marks)[sp.nok]

            # if there is some species missing in the table of traits
            if (len.sp.nok > 0) {

                # option 1) include missing species in the table with a mean trait value
                if (correct.trait == "mean") {
                  # repeat mean traits as many times as there are missing species
                  traitplus <- NULL
                  for (i in 1:len.sp.nok) traitplus <- rbind(traitplus,  meantrait)
                  row.names(traitplus) <- names.sp.nok
                  traits <- rbind(traits, traitplus)
                  traits <- traits[order(rownames(traits)), , drop = FALSE]
                  warning(paste("species", names.sp.nok,
                  "was not in the data.frame of traits\n and has been included with mean trait value\n"))
                }

                # option 2) eliminate missing species from the ppp
                if (correct.trait == "exclude") {
                  numbers.nok <- table(mippp[mippp$marks %in% names.sp.nok]$marks)[names.sp.nok]
                  mippp <- mippp[!mippp$marks %in% names.sp.nok]
                  mippp$marks <- factor(mippp$marks)
                  warning(paste("species ", names.sp.nok, ", with ", numbers.nok, 
                   " individuals, was not in the data.frame of traits and has been\n excluded from the analysis\n",
                     sep = ""))
                }
            }
            #traits se sustituye por la matriz de distancias de traits solo si vamos a calcular el IFDAR
            if(idar=="ifdar") traits <- as.matrix(gowdis(traits)) 
        }
        # if traits comes as a distance semimatrix (dist) transform it to a matrix
         # if (class(traits) == "dist")  traits <- as.matrix(traits) #!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! CHANGED 04/12/2019!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	 if (inherits(traits, what= "dist"))  traits <- as.matrix(traits)
        
       #if traits come as a matrix of traits instead of a data.frame and we are computing ifdar, do this:
        if (is.matrix(traits) & idar=="ifdar"){ 
             if (dim(traits)[1] != dim(traits)[2])  traits <- as.matrix(gowdis(traits))
         }
   return(traits) 
   }
   



