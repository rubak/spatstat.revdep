colormapgraph <- function( ncol = 3, mapval = NULL, notSeenAsBlack = TRUE,
                           txtfont = "sans", pointsize = 10,
                           symbol = "circle", size = 1, inch = 0.18 ) {

  lumth <- 0.4

  if( is.null( mapval ) ) {
    texteval <- "vfenv$nv$pmapsettings"
    mapval <- eval( parse( text = texteval ) )
  }

  total <- nrow( mapval )
  if( notSeenAsBlack ) total <- total + 1
  nrow <- ceiling( total / ncol )

# get coordinates to plot
  coords   <- NULL
  coords$x <- rep( 1:ncol, nrow )
  coords$y <- rep( 1:nrow, ncol )
  coords$y <- coords$y[order( coords$y, decreasing = TRUE )]
  coords   <- as.data.frame( coords )
  coords   <- coords[1:total,]

  xmin     <- min( coords$x ) - 1 / ncol
  xmax     <- max( coords$x ) + 1 / ncol
  ymin     <- min( coords$y ) - 1 / nrow
  ymax     <- max( coords$y ) + 1 / nrow
# get rgb and text to plot
  loccol <- NULL
  txtval   <- NULL
  idx <- 0
  if( notSeenAsBlack ) {
    idx <- 1
    loccol$red[idx]   <- 0
    loccol$green[idx] <- 0
    loccol$blue[idx]  <- 0
    txtval[idx]       <- "NS"
  }
  for( i in 1:nrow( mapval ) ) {
    loccol$red[i+idx]   <- mapval$red[i]
    loccol$green[i+idx] <- mapval$green[i]
    loccol$blue[i+idx]  <- mapval$blue[i]
    txtval[i+idx]       <- as.character( mapval$cutoffs[i] )
  }
  loccol <- as.data.frame( loccol )
  colval <- rgb( loccol )

# opar <- par( no.readonly = TRUE )
  oplt    <- par()$plt
  ops     <- par()$ps
  ofamily <- par()$family
  par( plt = c( 0, 1, 0, 1 ) )
  par( ps = pointsize )
  par( family = txtfont )

# legend
  plot( coords$x, coords$y, type = "n", axes = FALSE, xlab = "", ylab = "", xlim = c( xmin, xmax ), ylim = c( ymin, ymax ) )
  dimensions <- t( matrix( data = rep( size, nrow( coords ) ),nrow = length( size ), ncol = nrow( coords ) ) )
  evaltxt <- paste( "symbols( coords$x, coords$y, " , symbol, " = dimensions, add = TRUE, inches = inch, bg = colval, fg = colval, lwd = 1 )", sep = "" )
  eval( parse( text = evaltxt ) )

  coltxt <- rep( "black", length( coords$x ) )
  coltxt[( 0.2126 * loccol$red + 0.7152 * loccol$green + 0.0722 * loccol$blue ) < lumth] <- "white"
  coltxt[loccol$red < 0.1 & loccol$green < 0.6 & loccol$blue < 0.1] <- "white" # ad-hoc patch to make green scale look good
  text( coords$x, coords$y, txtval, col = coltxt )
  
  par( plt = oplt )
  par( ps = ops )
  par( family = ofamily )
}
