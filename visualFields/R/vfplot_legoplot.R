vfplot_legoplot <- function( vals, patternMap, vftiles, vfhull, loccolout, loccolin, radius = 2,
                             xmin = NULL, xmax = NULL, ymin = NULL, ymax = NULL,
                             txtfont = "sans", pointsize = 10,
                             showaxis = FALSE, colaxis = "black" ) {
  lumth <- 0.4
  
  # if not imposed, calculate limits of plot
  # expand by 5% each axis
  xrange <- max( patternMap$xod ) - min( patternMap$xod )
  yrange <- max( patternMap$yod ) - min( patternMap$yod )
  if( is.null( xmin ) ) xmin <- min( patternMap$xod ) - 0.025 * xrange
  if( is.null( xmax ) ) xmax <- max( patternMap$xod ) + 0.025 * xrange
  if( is.null( ymin ) ) ymin <- min( patternMap$yod ) - 0.025 * yrange
  if( is.null( ymax ) ) ymax <- max( patternMap$yod ) + 0.025 * yrange

  # init
  oplt    <- par()$plt
  ops     <- par()$ps
  ofamily <- par()$family
  par( plt = c( 0, 1, 0, 1 ) )
  par( ps = pointsize )
  par( family = txtfont )

  plot( patternMap$xod, patternMap$yod, type = "n", axes = FALSE, xlab = "", ylab = "", asp = 1, xlim = c( xmin, xmax ), ylim = c( ymin, ymax ) )
  plot( vftiles, clipp = vfhull, add = TRUE, close = TRUE, showpoints = FALSE, border = "darkgray", fillcol = rgb( loccolout ) )
  symbols( patternMap$xod, patternMap$yod, circles = rep( radius, nrow( patternMap ) ), inches = FALSE, fg = rgb( loccolin ), bg = rgb( loccolin ), add = TRUE )
  if( showaxis ) {
    axis( 1, pos = 0, labels = FALSE, lwd.ticks = 0, at = c( xmin, xmax ), col = colaxis )
    axis( 2, pos = 0, labels = FALSE, lwd.ticks = 0, at = c( ymin, ymax ), col = colaxis )
  }
  idx <- which( is.na( vals ) )
  if( length( idx ) > 0 ) {
    patternMap <- patternMap[-idx,]
    vals       <- vals[-idx]
  }

  coltxt <- rep( "grey10", length( patternMap$xod ) )
  coltxt[( 0.2126 * loccolin$red + 0.7152 * loccolin$green + 0.0722 * loccolin$blue ) < lumth] <- "white"
  text( patternMap$xod, patternMap$yod, vals, col = coltxt )
  par( plt    = oplt )
  par( ps     = ops )
  par( family = ofamily )
}