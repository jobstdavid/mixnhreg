#' @title The Mixture Normal Distribution
#' @description Density, distribution function, quantile function and random generation for the mixture normal distribution
#' as well as the corresponding \code{family.mixnhreg} object.
#' @name Mixture_Normal
#'
#' @param x,q vector of quantiles.
#' @param location matrix of location parameters.
#' @param scale matrix of scale parameters.
#' @param weight matrix of weight parameters.
#' @param log,log.p logical; if \code{TRUE}, probabilities p are given as log(p).
#'
#' @rdname Mixture_Normal
#' @export
dmixnorm <- function(x, location, scale, weight, log = FALSE) {

  if (length(unique(c(ncol(location), ncol(scale), ncol(weight)))) != 1) {
    stop("Number of columns in 'location', 'scale', 'weight' are not equal!")
  }
  d <- ncol(location)
  arg <- data.frame(x = x, location = location, scale = scale, weight = weight)

  dmixnorm_cpp(x = as.numeric(arg$x),
               location = as.matrix(arg[, 2:(d+1)]),
               scale = as.matrix(arg[, (d+2):(2*d+1)]),
               weight = as.matrix(arg[, (2*d+2):(3*d+1)]),
               log_p = log)
}


#' @param lower.tail logical; if \code{TRUE} (default), probabilities are \eqn{P[X \leq x]} otherwise, \eqn{P[X > x]}.
#' @rdname Mixture_Normal
#' @export
pmixnorm <- function(q, location, scale, weight, lower.tail = TRUE, log.p = FALSE) {

  if (length(unique(c(ncol(location), ncol(scale), ncol(weight)))) != 1) {
    stop("Number of columns in 'location', 'scale', 'weight' are not equal!")
  }
  d <- ncol(location)
  arg <- data.frame(q = q, location = location, scale = scale, weight = weight)

  pmixnorm_cpp(q = as.numeric(arg$q),
               location = as.matrix(arg[, 2:(d+1)]),
               scale = as.matrix(arg[, (d+2):(2*d+1)]),
               weight = as.matrix(arg[, (2*d+2):(3*d+1)]),
               log_p = log.p,
               lower_tail = lower.tail)
}

#' @param p vector of probabilities.
#' @rdname Mixture_Normal
#' @export
qmixnorm <- function(p, location, scale, weight, lower.tail = TRUE, log.p = FALSE) {

  if (length(unique(c(ncol(location), ncol(scale), ncol(weight)))) != 1) {
    stop("Number of columns in 'location', 'scale', 'weight' are not equal!")
  }
  d <- ncol(location)
  arg <- data.frame(p = p, location = location, scale = scale, weight = weight)

  qmixnorm_cpp(p = as.numeric(arg$p),
               location = as.matrix(arg[, 2:(d+1)]),
               scale = as.matrix(arg[, (d+2):(2*d+1)]),
               weight = as.matrix(arg[, (2*d+2):(3*d+1)]),
               log_p = log.p,
               lower_tail = lower.tail)
}

#' @param n number of observations. If \code{length(n) > 1}, the length is taken to be the number required.
#'
#' @rdname Mixture_Normal
#' @importFrom stats runif
#' @export
rmixnorm <- function(n, location, scale, weight) {

  if (length(unique(c(ncol(location), ncol(scale), ncol(weight)))) != 1) {
    stop("Number of columns in 'location', 'scale', 'weight' are not equal!")
  }
  d <- ncol(location)
  arg <- data.frame(p = runif(n), location = location, scale = scale, weight = weight)

  rmixnorm_cpp(p = as.numeric(arg$p),
               location = as.matrix(arg[, 2:(d+1)]),
               scale = as.matrix(arg[, (d+2):(2*d+1)]),
               weight = as.matrix(arg[, (2*d+2):(3*d+1)]))

}








