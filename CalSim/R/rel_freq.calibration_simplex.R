#Function returning the relative frequency in each bin.

rel_freq = function(x) {
  UseMethod("rel_freq")
}

rel_freq.calibration_simplex = function(x) {
  return(x$freq/x$n_obs)
}
