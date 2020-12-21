

###############################################################################
#
#    IFDAR.R
#
#   Individual Functional-Diversity -Area Relationship
#
#  $Version: 0.2 $  $Date: 2014/05/26 23:48:03 $
#
#
#   function to compute Indivdual Functional Diversity Area Relationships (i.e., acumulation of FDis curves
#    around individual species) in a community.
#
#    Value: an object of class fv , bassicaly a data.frame with the following
#    elemnts:
#     r: vector of radii at which IFDAR(r) has been estimated
#     ifdar: vector with the IFDAR(r) values
#     ###npoints: number of focal points employed to estimate isar(r) at each 
#              radius r
#
#    Arguments:
#     mippp:    Multivariate point pattern (ppp of spatstat)
#     mippp.sp: Focal (unmarked) point pattern (ppp of spatstat) or character
#                string indicating one of the marks in the multivariate mippp 
#     mimark:   character string indicating one of the marks in the
#                multivariate mippp
#     namesmark: name of the column with species names in the data.frame of
#                marks (only for multimarked ppp's)
#     traits:    A data.frame of traits with cols=traits, rows = species. rownames should be equal to species names in mippp 
#     r:        vector with the sequence of radii (>0) at which estimate
#                IFDAR(r)
#     buffer:   a number indicating the width around the window of the focal 
#                point pattern that will be excluded from the computations of
#                IFDAR(r) to control de edge effect, or the string "adapt", 
#                indicating that an adaptive border (of width ri) will be 
#                excluded for the computation of each ISAR(ri)
#     bfw:      an owin object indicating a subset of the focal point pattern 
#                that will be employed to compute ISAR(r)
#
#   Details:
#     bfw can be employed to select different habitats and compute IFDAR((R) in 
#      each habitat


 
ifdar<- function(mippp, mippp.sp=NULL, mimark=NULL,  namesmark=NULL, traits=NULL,
                      r=NULL, buffer=0, bfw=NULL, correct.trait.na=FALSE, correct.trait="mean") {
    if (is.null(traits))   stop("you should provide a dataframe of traits  to compute ifdar")
    
    #  if (is.data.frame(traits))  traits <- gowdis(traits)
    #  if (class(traits) == "dist")  traits <- as.matrix(traits)
    #  if (is.matrix(traits)) {
    #         if (dim(traits)[1] != dim(traits)[2])  traits <- as.matrix(gowdis(traits))
    #  }

    # Check mippp
    if (!is.null(namesmark))  mippp$marks <- factor((mippp$marks[namesmark][[1]]))

    # check traits
    idar <-"ifdar"
    traits <-checktraits(traits=traits, mippp=mippp, idar=idar, correct.trait.na=correct.trait.na, correct.trait=correct.trait)

    # Check mippp.sp
    if (!is.null(mippp.sp) & !is.ppp(mippp.sp)) {
        mippp.sp <- NULL
        mimark <- mippp.sp
    }
    if (!is.null(mimark)) 
        if (mimark %in% levels(mippp$marks) == FALSE) {
            stop(paste(mimark, " can't be recognized as a mark\n\n\n
                    have you indicated in which column of thedataframe are the species\n
                    marks? (argument 'namesmark'\n\n"))
        }
    if (is.null(mippp.sp)) 
        mippp.sp <- mippp[mippp$marks == mimark]
    if (buffer != "adapt") {
        if (is.null(bfw)) 
            bfw <- erosion(mippp$window, buffer)
        mippp.sp <- mippp.sp[inside.owin(mippp.sp, w = bfw)]
        npoints <- rep(mippp.sp$n, length(r))
        names(npoints) <- r
    }
    cosamt <- mitable(mippp.sp, mippp, r)
    if (buffer == "adapt") {
        bdp <- bdist.points(mippp.sp)
        for (i in 1:length(r)) cosamt[[i]] <- cosamt[[i]][bdp >= 
            r[i], ]
        npoints <- sapply(cosamt, function(x) dim(x)[1])
    }
    ifdar <- sapply(cosamt, function(x) fdis(x, traits = traits)) 
    result <- data.frame(r = r, ifdar = ifdar)
    result <- fv(result, argu = "r", ylab = substitute(IFDAR(r), 
        NULL), valu = "ifdar", fmla = ifdar ~ r, alim = c(min(r), 
        max(r)), labl = c("r", "%s(r)"), desc = c("radius of circle", 
        "%s"), fname = "IFDAR")
    return(result)
}
###############################################################################