#
# interactive plot for lpp objects
#

iplot.lpp <- function(x, ..., xname) {
  if(missing(xname))
    xname <- short.deparse(substitute(x))
  stopifnot(is.lpp(x))
  ## predigest
  L <- domain(x)
  v <- vertices(L)
  deg <- vertexdegree(L)
  dv <- textstring(v, txt=paste(deg))
  y <- layered(lines=as.psp(L),
               vertices=v,
               degree=dv,
               points=as.ppp(x))
  iplot(y, ..., xname=xname, visible=c(TRUE, FALSE, FALSE, TRUE))
}

#
# interactive plot for linnet objects
#

iplot.linnet <- function(x, ..., xname) {
  if(missing(xname))
    xname <- short.deparse(substitute(x))
  if(!inherits(x, "linnet"))
    stop("x should be a linnet object")
  ## predigest
  v <- vertices(x)
  deg <- vertexdegree(x)
  dv <- textstring(v, txt=paste(deg))
  y <- layered(lines=as.psp(x),
               vertices=v,
               degree=dv)
  iplot(y, ..., xname=xname, visible=c(TRUE, FALSE, FALSE))
}

