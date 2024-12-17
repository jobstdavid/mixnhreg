#' Auxiliary Function for Controlling the Gradient-Boosting Estimation
#' @description Auxiliary function for controlling the non-cylcic gradient-boosting estimation (see references).
#'
#' @param maxit integer; giving the number of initial boosting iterations. Default: \code{maxit = 6000}.
#' @param nu double; defining the step size or shrinkage parameter. Default: \code{nu = 0.1}.
#' @param mstop character; stopping criterion "\code{max}", "\code{aic}" (default), "\code{bic}", "\code{cv}".
#' @param nfolds integer; number of folds in cross validation.
#' @param foldid vector; an optional vector of values between 1 and \code{nfolds} identifying the fold for each observation.
#' @param cores integer; the number of cores used for cross validation. Default: \code{cores = 1L}.
#'
#' @return An object of class \code{mixnhreg.boost}.
#'
#' @rdname control_boost
#' @export
control_boost <- function(maxit = 6000,
                          nu = 0.1,
                          mstop = "aic",
                          nfolds = 10,
                          foldid = NULL,
                          cores = 1L) {

  out <- list(maxit = maxit,
              nu = nu,
              mstop = mstop,
              nfolds = nfolds,
              foldid = foldid,
              cores = cores)
  class(out) <- "mixnhreg.boost"

  return(out)

}
