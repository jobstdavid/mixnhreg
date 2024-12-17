#' @title Scoring Rules for Mixture Normal Distribution
#' @description Logarithmic score (logs) a.k.a. negative log-likelihood and continuous ranked probability score (crps) for the mixture normal distribution.
#' @name scores_mixnorm
#'
#' @param x vector of quantiles.
#' @param location matrix of location parameters.
#' @param scale matrix of scale parameters.
#' @param weight matrix of weight parameters.
#'
#' @references
#'
#' Grimit E.P., Gneiting T., Berrocal V., Johnson N.A. (2006). The continuous ranked probability score for circular variables and its application to mesoscale forecast ensemble verification.
#' Quarterly Journal of the Royal Meteorological Society, 132, pp. 2925–2942. doi: <https://doi.org/10.1256/qj.05.235>.
#'
#' @rdname scores_mixnorm
#' @export
logs_mixnorm <- function(x, location, scale, weight) {

  if (length(unique(c(ncol(location), ncol(scale), ncol(weight)))) != 1) {
    stop("Number of columns in 'location', 'scale', 'weight' are not equal!")
  }
  d <- ncol(location)
  arg <- data.frame(x = x, location = location, scale = scale, weight = weight)

  logs_mixnorm_cpp(x = as.numeric(arg$x),
                  location = as.matrix(arg[, 2:(d+1)]),
                  scale = as.matrix(arg[, (d+2):(2*d+1)]),
                  weight = as.matrix(arg[, (2*d+2):(3*d+1)]))
}

#' @rdname scores_mixnorm
#' @export
crps_mixnorm <- function(x, location, scale, weight) {

  if (length(unique(c(ncol(location), ncol(scale), ncol(weight)))) != 1) {
    stop("Number of columns in 'location', 'scale', 'weight' are not equal!")
  }
  d <- ncol(location)
  arg <- data.frame(x = x, location = location, scale = scale, weight = weight)

  crps_mixnorm_cpp(x = as.numeric(arg$x),
                   location = as.matrix(arg[, 2:(d+1)]),
                   scale = as.matrix(arg[, (d+2):(2*d+1)]),
                   weight = as.matrix(arg[, (2*d+2):(3*d+1)]))
}

grad_logs_mixnorm <- function(x, location, scale, weight) {

  if (length(unique(c(ncol(location), ncol(scale), ncol(weight)))) != 1) {
    stop("Number of columns in 'location', 'scale', 'weight' are not equal!")
  }
  d <- ncol(location)
  arg <- data.frame(x = x, location = location, scale = scale, weight = weight)

  grad_logs_mixnorm_cpp(x = as.numeric(arg$x),
                        location = as.matrix(arg[, 2:(d+1)]),
                        scale = as.matrix(arg[, (d+2):(2*d+1)]),
                        weight = as.matrix(arg[, (2*d+2):(3*d+1)]))
}

grad_crps_mixnorm <- function(x, location, scale, weight) {

  if (length(unique(c(ncol(location), ncol(scale), ncol(weight)))) != 1) {
    stop("Number of columns in 'location', 'scale', 'weight' are not equal!")
  }
  d <- ncol(location)
  arg <- data.frame(x = x, location = location, scale = scale, weight = weight)

  grad_crps_mixnorm_cpp(x = as.numeric(arg$x),
                        location = as.matrix(arg[, 2:(d+1)]),
                        scale = as.matrix(arg[, (d+2):(2*d+1)]),
                        weight = as.matrix(arg[, (2*d+2):(3*d+1)]))
}
