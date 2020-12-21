  #fdis: function to compute IFDAR (IFDisp) --------------------------------------------------------------------
  # x is a community data.frame.
  # traits is a distance matrix (as a matrix, i.e., from as.matrix(dist(traits)) or as.matrix(gowdis(traits))

 fdis <- function(x, traits){
     sp.a <- colSums(x) > 0
    if (sum(sp.a) < 1) 
        return(NA)
    if (sum(sp.a) == 1) 
        return(0)
    if (sum(sp.a) > 1) {
        m.ok <- x[, sp.a]
        sp.trait.ok <- !is.na(match(rownames(traits), colnames(m.ok)))
        if (!is.null(dim(m.ok))) {
            com.G.0 <- rowSums(m.ok) > 0
            cosad <- as.dist(traits[sp.trait.ok, sp.trait.ok, 
                drop = FALSE])
            if (sum(cosad, na.rm = TRUE) == 0) 
                return(0)
            if (sum(sp.trait.ok) >= 2) {
                 com.G.02 <- rowSums(m.ok[,labels(cosad)]) > 0 
                cosaf <- fdisp(d = cosad, a = m.ok[com.G.02, labels(cosad), 
                  drop = FALSE], tol = 1e-07)$FDis
            }
            if (sum(sp.trait.ok) < 2) {
                cosaf <- fdisp(d = cosad, a = t(as.matrix(m.ok[com.G.0, 
                  labels(cosad)])), tol = 1e-07)$FDis
            }
            media <- mean(cosaf, na.rm = T)
            result <- media
        }
        if (is.null(dim(m.ok))) 
            result <- NA
        return(result)
    }
}
   # ---- end function fdis ------------------------------------------------------------------------------------------------------------
