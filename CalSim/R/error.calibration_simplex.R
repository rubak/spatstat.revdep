#Function to compute the miscalibration error in each bin used to construct the
#calibration simplex.

error = function(x,true_error) {
  UseMethod("error")
}

error.calibration_simplex = function(x,
                                     true_error = TRUE) {

  error = matrix(rep(NA,x$n_bins*3),ncol = 3)
  # error = data.frame(c1 = rep(NA,x$n_bins),
  #                    c3 = rep(NA,x$n_bins))
  if(true_error) {

    error[x$freq > 0,] = x$cond_rel_freq[x$freq > 0,] - x$cond_ave_prob[x$freq > 0,]
    # error$c3[x$freq > 0] = x$cond_rel_freq_3[x$freq > 0] - x$cond_p3_ave[x$freq > 0]
    # error$c1[x$freq > 0] = x$cond_rel_freq_1[x$freq > 0] - x$cond_p1_ave[x$freq > 0]
  }
  else {
    rounded_forecasts = make_forecasts(x)

    error[x$freq > 0,] = x$cond_rel_freq[x$freq > 0,] - rounded_forecasts[x$freq > 0,]
    # error$c3[x$freq > 0] = x$cond_rel_freq_3[x$freq > 0] - rounded_forecasts$p3[x$freq > 0]
    # error$c1[x$freq > 0] = x$cond_rel_freq_1[x$freq > 0] - rounded_forecasts$p1[x$freq > 0]
  }
  return(error)
}
