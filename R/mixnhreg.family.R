#' Family objects in \pkg{mixnhreg}
#'
#' @description In the \pkg{mixnhreg} package, family objects are utilized to define the necessary information for fitting a mixture model via \code{\link{mixnhreg}}.
#' This encompasses specifics like the distribution/family name, link functions for the parameters, density function, negative log-likelihood/crps function, and derivatives of the negative log-likelihood/crps with respect to the predictors.
#' Additionally, these family objects are required for predictions via \code{\link{predict.mixnhreg}}.
#'
#' @details
#'  The following lists the minimum requirements on a \pkg{mixnhreg} family object to be used with \pkg{mixnhreg}:
#'  \itemize{
#'    \item The family object is expected to return a \code{\link{list}} of class "\code{mixnhreg.family}".
#'    \item The object must contain the family name ("\code{family}") and its abbreviation ("\code{dist}") as a character string.
#'    \item The object must contain link functions ("\code{links}") for the location and scale parameters as character string.
#'    \item The family object must contain either the negative log-likelihood ("\code{nll}") function or the continuous ranked probability score ("\code{crps}") function
#'    incl. its corresponding gradient function ("\code{grad_nll}", "\code{grad_crps}") with respect to the predictors.
#'  }
#' Remaining arguments in the family object, such as the probability density function ("\code{density}"), cumulative distribution function ("\code{probability}") or
#' quantile function ("\code{quantile}") are not necessary for the model estimation but might be used for predictions via \code{\link{predict.mixnhreg}}.
#'
#' @param ... additional arguments passed to the family object.
#'
#' @return An object of class \code{mixnhreg.family}.
#'
#' @name mixnhreg.family
#' @rdname mixnhreg.family
#'
#' @export
mixnorm <- function(...) {

  f <- list("family" = "Mixture Normal Distribution",
            "dist" = "mixnorm",
            "links" = c(location = "id", scale = "log", weight = "softmax"),
            "crps" = function(x, location, scale, weight) {
              crps_mixnorm(x, location = location, scale = scale, weight = weight)
            },
            "nll" = function(x, location, scale, weight) {
              logs_mixnorm(x, location = location, scale = scale, weight = weight)
            },
            "grad_crps" = function(x, location, scale, weight) {
              grad_crps_mixnorm(x, location = location, scale = scale, weight = weight)
            },
            "grad_nll" = function(x, location, scale, weight) {
              grad_logs_mixnorm(x, location = location, scale = scale, weight = weight)
            },
            "density" = function(x, location, scale, weight) {
              dmixnorm(x, location = location, scale = scale, weight = weight)
            },
            "probability" = function(q, location, scale, weight) {
              pmixnorm(q, location = location, scale = scale, weight = weight)
            },
            "quantile" = function(p, location, scale, weight) {
              qmixnorm(p, location = location, scale = scale, weight = weight)
            }
  )


  class(f) <- "mixnhreg.family"

  return(f)

}
