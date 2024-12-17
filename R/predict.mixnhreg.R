#' Prediction for \pkg{mixnhreg}
#' @description This is the main function of \pkg{mixnhreg} to get predictions from an \code{mixnhreg} object.
#'
#' @param object an object of class \code{mixnhreg}.
#' @param newdata a data frame containing the predictors at which predictions should be made.
#' @param type type of prediction:
#' \itemize{
#'  \item "\code{location}", "\code{scale}", "\code{weight}" returns the respective distribution parameters.
#'  \item "\code{parameter}" returns a list of all distribution parameters.
#'  \item "\code{density}" evaluates the probability density function (PDF) at "\code{response}".
#'  \item "\code{probability}" evaluates the cumulative distribution function (CDF) at "\code{response}".
#'  \item "\code{quantile}" evaluates the quantile function (QF) at "\code{response}".
#'  \item "\code{crps}" returns the CRPS of the distribution evaluated at "\code{response}".
#'  \item "\code{nll}" returns the negative log-likelihood of the distribution evaluated at "\code{response}".
#' }
#' @param response a vector of values needed to evaluate the distribution for a specified \code{type}:
#' "\code{density}", "\code{probability}", "\code{quantile}", "\code{crps}", "\code{nll}".
#' @param ... unused.
#'
#' @return a list or data frame of predictions.
#'
#' @examples
#' # load data
#' data("station")
#'
#' # fit mixture normal distribution with two components via BFGS
#' (fit_optim <- mixnhreg(formula = obs ~ sin1 + cos1 + temp_mean | temp_ctrl,
#'                        scale.formula = ~ sin1 + cos1 + log(temp_sd) | 1,
#'                        weight.formula = ~ sin1 + cos1 | 1,
#'                        data = station,
#'                        control = control_optim()))
#'
#' # predict location parameter
#' location <- predict(fit_optim, type = "location")
#'
#' # predict all distribution parameter
#' parameter <- predict(fit_optim, type = "parameter")
#'
#' # predict negative log-likelihood
#' nll <- predict(fit_optim, newdata = station, type = "nll", response = station$obs)
#'
#'
#' @rdname predict.mixnhreg
#' @export

predict.mixnhreg <- function(object,
                             newdata,
                             type = "parameter",
                             response = 0.5,
                             ...) {

  # number of components
  formula <- object$formula
  scale.formula <- object$scale.formula
  weight.formula <- object$weight.formula
  f <- Formula(formula)
  K <- length(f)[2]

  # get model frame
  fake_formula <- reformulate(response = as.character(formula[2]), termlabels = c(attr(terms(as.Formula(formula)), "term.labels"),
                                                                                  attr(terms(as.Formula(scale.formula)), "term.labels"),
                                                                                  attr(terms(as.Formula(weight.formula)), "term.labels")))
  if (missing(newdata)) {
    mf <- object$model_frame
  } else {
    mf <- model.frame(fake_formula,
                      newdata,
                      na.action = na.pass)
  }
  N <- nrow(mf)

  # predictor matrices
  X <- c(lapply(1:K, function(j) {
    model.matrix(formula(as.Formula(formula), lhs = F, rhs = j), mf)}),
    lapply(1:K, function(j) {
      model.matrix(formula(as.Formula(scale.formula), lhs = F, rhs = j), mf)}),
    lapply(1:K, function(j) {
      model.matrix(formula(as.Formula(weight.formula), lhs = F, rhs = j), mf)}))

  # coefficients
  x <- object$coef

  # get link functions
  links <- as.character(object$family$links)

  # get predictor
  if (inherits(object, "mixnhreg.optim")) {
    parameter <- get_parameter_optim_cpp(X, unlist(x), links, N, K)
  } else if (inherits(object, "mixnhreg.boost")) {
    parameter <- get_parameter_boost_cpp(X, x, links, N, K)
  }
  names(parameter) <- c("location", "scale", "weight")

  # combine arguments
  args <- append(list(response), parameter)
  names(args) <- NULL

  out <-  switch(type,
                 "location" = data.frame(location = parameter$location),
                 "scale" = data.frame(scale = parameter$scale),
                 "weight" = data.frame(weight = parameter$weight),
                 "parameter" = parameter,
                 "density" = do.call(object$family$density, args),
                 "probability" = do.call(object$family$probability, args),
                 "quantile" = do.call(object$family$quantile, args),
                 "crps" = do.call(object$family$crps, args),
                 "nll" = do.call(object$family$nll, args))

  return(out)

}
