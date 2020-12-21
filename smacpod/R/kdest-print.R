#' Print a \code{kdenv} object
#'
#' Print an \code{kdenv} object produced by
#' \code{\link[smacpod]{kdest}}.
#'
#' @param x An object produced by the
#'   \code{\link[smacpod]{kdest}} function.
#' @param ... Not currently implemented.
#' @param extra A logical value indicating whether extra
#'   information related to the internal
#'   \code{\link[spatstat.core]{fv}} object should be printed.
#'   The default is \code{FALSE}.
#' @return Information about the \code{kdest}
#' @author Joshua French
#' @export
print.kdenv = function(x, ..., extra = FALSE) {
  cat("Difference in K functions\n")
  cat("case label: ", x$case_label, "\n")
  cat("control label: ", x$control_label, "\n")
  if (x$nsim > 0) {
    cat("number of simulations: ", x$nsim, "\n")
    cat("level: ", x$level, "\n")
  }
  if (extra) {
    x$out  
  }
}
