
midar<-
function (mippp, mippp.sp = NULL, mimark = NULL, namesmark = NULL, 
    traits = NULL, tree=NULL, r = NULL, buffer = 0, bfw = NULL, what=NULL) 
{
    if(!is.expression(what)& !is.function(what)) stop("'what' should be a valid R function or expression")
    
    # TODO: include a test that the expression "what" will accept a matrix and will return a unique value
   
    if (!is.null(namesmark)) 
        mippp$marks <- factor((mippp$marks[namesmark][[1]]))
   
   
    if (!is.null(mippp.sp) & !is.ppp(mippp.sp)) {
        mippp.sp <- NULL
        mimark <- mippp.sp
    }
    if (!is.null(mimark)) 
        if (mimark %in% levels(mippp$marks) == FALSE) {
            stop(paste(mimark, " can't be recognized as a mark\n\n\n\n                    have you indicated in which column of thedataframe are the species\n\n                    marks? (argument 'namesmark'\n\n"))
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
   
    if(is.expression(what)) midar <- sapply(cosamt, function(x) eval(what))
    if(is.function(what))   midar <- sapply(cosamt, what)
    
    result <- data.frame(r = r, midar = midar)
    #print(result)
    result <- fv(result, argu = "r", ylab = substitute(MIDAR(r), 
        NULL), valu = "midar", fmla = midar ~ r, alim = c(min(r), 
        max(r)), labl = c("r", "%s(r)"), desc = c("radius of circle", 
        "%s"), fname = "MIDAR")
    return(result)
}
