vfplotloc <- function( vals, patternMap, loccol, vftiles, vfhull,
                       xmin, xmax, ymin, ymax,
                       txtfont = "sans", pointsize = 10,
                       showaxis = FALSE, colaxis = "black" ) {
  lumth <- 0.4

# init
  oplt    <- par()$plt
  ops     <- par()$ps
  ofamily <- par()$family
  par( plt = c( 0, 1, 0, 1 ) )
  par( ps = pointsize )
  par( family = txtfont )

  plot( patternMap$xod, patternMap$yod, type = "n", axes = FALSE, xlab = "", ylab = "", asp = 1, xlim = c( xmin, xmax ), ylim = c( ymin, ymax ) )
  plot( vftiles, clipp = vfhull, add = TRUE, close = TRUE, showpoints = FALSE, border = "darkgray", fillcol = rgb( loccol ) )
  if( showaxis ) {
    axis( 1, pos = 0, labels = FALSE, lwd.ticks = 0, at = c( xmin, xmax ), col = colaxis )
    axis( 2, pos = 0, labels = FALSE, lwd.ticks = 0, at = c( ymin, ymax ), col = colaxis )
  }
  coltxt <- rep( "grey10", length( patternMap$xod ) )
  coltxt[( 0.2126 * loccol$red + 0.7152 * loccol$green + 0.0722 * loccol$blue ) < lumth] <- "white"
  coltxt[loccol$red < 0.1 & loccol$green < 0.6 & loccol$blue < 0.1] <- "white" # ad-hoc patch to make green scale look good
  text( patternMap$xod, patternMap$yod, vals, col = coltxt )
  par( plt    = oplt )
  par( ps     = ops )
  par( family = ofamily )
}