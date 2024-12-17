#' Auxiliary Function for Controlling Optimization
#' @description Auxiliary function for controlling the optimization.
#'
#' @param method character; optimization method passed to \code{\link{optim}}. Default: "\code{BFGS}".
#' @param maxit integer; maximum number of iterations. Default: \code{100}.
#' @param start list; list of initial values for the location, scale and weight parameters.
#' If a parameter is set to \code{NULL} (default), corresponding initial values are internally determined.
#' @param hessian logical; if \code{TRUE} the numerical Hessian matrix from the \code{\link{optim}} output is used for the covariance matrix estimation.
#' If \code{FALSE} (default) no covariance matrix is determined.
#' @param trace logical; if \code{TRUE}, the optimization is traced, otherwise \code{FALSE} (default) not.
#' @param lower integer; lower bound for the "\code{L-BFGS-B}", "\code{Brent}" methods.
#' @param upper integer; upper bound for the "\code{L-BFGS-B}", "\code{Brent}" methods.
#'
#' @return An object of class \code{mixnhreg.optim}.
#'
#' @rdname control_optim
#' @export
control_optim <- function(method = "BFGS",
                          maxit = 100,
                          start = list(location = NULL, scale = NULL, weight = NULL),
                          hessian = FALSE,
                          trace = FALSE,
                          lower = -Inf,
                          upper = Inf) {

  out <- list(method = method,
              maxit = maxit,
              start = start,
              hessian = hessian,
              trace = trace,
              lower = lower,
              upper = upper)
  class(out) <- "mixnhreg.optim"

  return(out)

}
