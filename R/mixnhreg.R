#' Mixtures of Non-Homogeneous Linear Regression Models
#' @description This is the main function of \pkg{mixnhreg} for estimating mixtures of non-homogeneous linear regression models.
#'
#' @param formula a formula object, with the response on the left of the \code{~} operator followed by the mixture components separated by an \code{|} operator,
#' where the predictors related to the location parameter of each mixture component are split up by the \code{+} operator.
#' @param scale.formula a formula object starting with the \code{~} operator followed by the mixture components separated by an \code{|} operator,
#' where the predictors related to the scale parameter of each mixture component are split up by the \code{+} operator.
#' @param weight.formula a formula object starting with the \code{~} operator followed by the mixture components separated by an \code{|} operator,
#' where the predictors related to the weight parameter of each mixture component are split up by the \code{+} operator.
#' @param data a data frame containing the variables occurring in the formulas.
#' @param na.action a function which indicates what should happen when the data contains \code{NA}s. Default: \code{na.omit}.
#' @param family a \code{mixnhreg.family} object, specifying details of the modeled distribution. Default: \code{\link{mixnorm}()}.
#' @param loss loss function used for optimization, which is either negative log-likelihood "\code{nll}" (default) or
#' continuous ranked probability score "\code{crps}".
#' @param control either \code{\link{control_optim}()} (default) for estimation via \code{\link{optim}} or \code{\link{control_boost}()} for estimation via gradient-boosting.
#' @param ... unused.
#'
#' @return An object of class \code{mixnhreg}.
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
#' # fit mixture normal distribution with two components via gradient-boosting
#' (fit_boost <- mixnhreg(formula = obs ~ sin1 + cos1 + temp_mean | temp_ctrl,
#'                        scale.formula = ~ sin1 + cos1 + log(temp_sd) | 1,
#'                        weight.formula = ~ sin1 + cos1 | 1,
#'                        data = station,
#'                        control = control_boost(mstop = "cv")))
#'
#' # model summary
#' summary(fit_optim)
#' summary(fit_boost)
#'
#' # parameter predictions
#' par_optim <- predict(fit_optim, type = "parameter")
#' par_boost <- predict(fit_boost, type = "parameter")
#'
#' # get observations and mean forecasts
#' obs <- na.omit(station$obs)
#' mean_optim <- rowSums(par_optim$location * par_optim$weight)
#' mean_boost <- rowSums(par_boost$location * par_boost$weight)
#'
#' # residual plot
#' plot(obs - mean_optim,
#'      xlab = "Index",
#'      ylab = "Observation - Mean",
#'      pch = 19,
#'      col = adjustcolor("red", alpha = 0.5))
#' points(obs - mean_boost,
#'        pch = 19,
#'        col = adjustcolor("steelblue", alpha = 0.25))
#' legend("bottomright",
#'        legend = c("optim", "boost"),
#'        pch = 19,
#'        col = c("red", "steelblue"))
#' grid()
#'
#'
#' @references
#' Hepp, T., J. Zierk, and E. Bergherr (2023). "Component-wise boosting for mixture distributional regression models".
#' In: Proceedings of the 37th International Workshop on Statistical Modelling. url: https://iwsm2023.statistik.tu- dortmund.de/.
#'
#' Jobst, D. (2024). "Gradient-Boosted Mixture Regression Models for Postprocessing Ensemble Weather Forecasts". doi: <https://doi.org/10.48550/arXiv.2412.09583>.
#'
#' Messner, J. W., G. J. Mayr, and A. Zeileis (2017). "Nonhomogeneous Boosting for Predictor Selection in Ensemble Postprocessing".
#' In: Monthly Weather Review 145.1, pp. 137–147. doi: <https://doi.org/10.1175/mwr-d-16-0088.1>.
#'
#' Thomas, J. et al. (2017). "Gradient boosting for distributional regression: faster tuning and improved variable selection via noncyclical updates".
#' In: Statistics and Computing 28.3, pp. 673–687. doi: <https://doi.org/10.1007/s11222-017-9754-6>.
#'
#' @rdname mixnhreg
#' @importFrom Formula Formula as.Formula
#' @importFrom stats formula model.frame model.matrix model.response reformulate as.formula update.formula terms lm optim sigma na.omit na.pass var
#' @export
mixnhreg <- function(formula,
                     scale.formula,
                     weight.formula,
                     data,
                     na.action = na.omit,
                     family = mixnorm(),
                     loss = "nll",
                     control = control_optim(),
                     ...) {

  # number of components
  f <- Formula(formula)
  K <- length(f)[2]

  # get model frame
  fake_formula <- reformulate(response = as.character(f[2]), termlabels = c(attr(terms(f), "term.labels"),
                                                                            attr(terms(as.Formula(scale.formula)), "term.labels"),
                                                                            attr(terms(as.Formula(weight.formula)), "term.labels")))

  mf <- model.frame(fake_formula,
                    data,
                    na.action = na.action)
  N <- nrow(mf)

  # response
  y <- as.numeric(model.response(mf))
  # predictor matrices
  X <- c(lapply(1:K, function(j) {
    model.matrix(formula(as.Formula(formula), lhs = F, rhs = j), mf)}),
    lapply(1:K, function(j) {
      model.matrix(formula(as.Formula(scale.formula), lhs = F, rhs = j), mf)}),
    lapply(1:K, function(j) {
      model.matrix(formula(as.Formula(weight.formula), lhs = F, rhs = j), mf)}))


  if (inherits(control, "mixnhreg.optim")) {

    # determine initial values for x as list
    # initial coefficients for location parameter
    if (is.null(control$start$location)) {
      start <- lapply(1:K, function(j) {
        init <- as.numeric(coef(lm(formula = formula(f, lhs = T, rhs = j),
                                  data = data,
                                  na.action = na.action)))
        names(init) <- colnames(X[[j]])
        init
      })
    } else {
      start <- control$start$location
    }
    # initial coefficients for scale parameter
    if (is.null(control$start$scale)) {
      start.scale <- lapply(1:K, function(j) {
        if ("(Intercept)" %in% colnames(X[[K+j]])) {
          scale <- sigma(lm(formula = formula(f, lhs = T, rhs = j),
                            data = data,
                            na.action = na.action))
          init <- c(as.vector(get_parameter_cpp(as.matrix(scale), family$links[2], FALSE)), rep(0, ncol(X[[K+j]])-1))
        } else {
          init <- rep(0, ncol(X[[K+j]]))
        }
        names(init) <- colnames(X[[K+j]])
        init
      })
      start <- c(start, start.scale)
    } else {
      start <- c(start, control$start$scale)
    }
    # initial coefficients for weight parameter
    if (is.null(control$start$weight)) {
      start.weight <- lapply(1:K, function(j) {
        init <- rep(0, ncol(X[[2*K+j]]))
        names(init) <- colnames(X[[2*K+j]])
        init
      })
      start <- c(start, start.weight)
    } else {
      start <- c(start, control$start$weight)
    }

    # get link functions
    links <- as.character(family$links)

    # get loss function
    loss_fun <- switch(loss,
                       "crps" = family$crps,
                       "nll" = family$nll)
    # get gradient of loss function
    grad_fun <- switch(loss,
                       "crps" = family$grad_crps,
                       "nll" = family$grad_nll)

    # optimization function
    fn <- function(x, y, X, links, loss_fun, K, N, ...) {

      # get parameter
      parameter <- get_parameter_optim_cpp(X, x, links, N, K)

      sum(loss_fun(x = y, location = parameter[[1]], scale = parameter[[2]], weight = parameter[[3]]))

    }
    # gradient of optimization function
    gr <- function(x, y, X, links, grad_fun, K, N, ...) {

      # get parameter
      parameter <- get_parameter_optim_cpp(X, x, links, N, K)

      g <- do.call("cbind", grad_fun(x = y, location = parameter[[1]], scale = parameter[[2]], weight = parameter[[3]]))

      grad_optim_cpp(X, g)

    }

    fit <- suppressWarnings(optim(par = unlist(start),
                                  fn = fn,
                                  gr = gr,
                                  y = y,
                                  X = X,
                                  links = links,
                                  loss_fun = loss_fun,
                                  grad_fun = grad_fun,
                                  K = K,
                                  N = N,
                                  method = control$method,
                                  lower = control$lower,
                                  upper = control$upper,
                                  hessian = control$hessian,
                                  control = list(trace = control$trace,
                                                 maxit = control$maxit)))


    coef <- update_list(start, fit$par, K)
    names(coef) <- rep(c("location", "scale", "weight"), each = K)

    # prepare output object
    out <- list(formula = formula,
                scale.formula = scale.formula,
                weight.formula = weight.formula,
                model_frame = mf,
                coef = coef,
                family = family,
                loss = loss,
                optim = fit,
                control = control)
    class(out) <- c("mixnhreg.optim", "mixnhreg")

    # add fit statistics
    loglik <- -sum(predict.mixnhreg(out, type = "nll", response = out$model_frame[, 1]))
    df <- length(unlist(out$coef))
    stats <- list(n = N,
                  df = df,
                  CRPS = sum(predict.mixnhreg(out, type = "crps", response = out$model_frame[, 1])),
                  logLik = loglik,
                  AIC = -2*loglik + 2*df,
                  BIC = -2*loglik + df*log(N))
    out <- append(out,
                  list(stats = stats),
                  after = 9)
    class(out) <- c("mixnhreg.optim", "mixnhreg")


  }

  if (inherits(control, "mixnhreg.boost")) {

    # standardize predictor matrices
    X <- lapply(1:length(X), function(k) standardize.matrix(X[[k]]))
    # standardize response
    y <- standardize.matrix(y, center = mean(y), scale = sqrt(var(y)*(N-1)/N))

    # save center and scale
    standardize <- lapply(1:length(X), function(k) {
      list(center = attr(X[[k]], "standardize:center"), scale = attr(X[[k]], "standardize:scale"))
    })
    standardize <- append(standardize,
                          list(list(center = c("(Intercept)" = attr(y, "standardize:center")), scale = c("(Intercept)" = attr(y, "standardize:scale")))),
                          after = 0)

    # coefficients
    x <- c(lapply(1:K, function(j) {rep(0, ncol(X[[j]]))}),
           lapply(1:K, function(j) {rep(0, ncol(X[[K+j]]))}),
           lapply(1:K, function(j) {rep(0, ncol(X[[2*K+j]]))}))

    # get link functions
    links <- as.character(family$links)

    # implemented family discrimination for faster estimation
    if (family$dist == "mixnorm") {

      # get loss function
      loss_fun <- switch(loss,
                         "crps" = crps_mixnorm_cpp,
                         "nll" = logs_mixnorm_cpp)
      # get gradient of loss function
      grad_fun <- switch(loss,
                         "crps" = grad_crps_mixnorm_cpp,
                         "nll" = grad_logs_mixnorm_cpp)
      # get negative logLik function
      nll_fun <- logs_mixnorm_cpp

    } else {

      # get loss function
      loss_fun <- switch(loss,
                         "crps" = family$crps,
                         "nll" = family$nll)
      # get gradient of loss function
      grad_fun <- switch(loss,
                         "crps" = family$grad_crps,
                         "nll" = family$grad_nll)

      # get negative logLik function
      nll_fun <- family$nll

    }

    # initial gradient-boosting
    fit <- boost_noncylcic_cpp(y,
                               X,
                               x,
                               links,
                               loss_fun,
                               grad_fun,
                               nll_fun,
                               control$nu,
                               control$maxit,
                               N,
                               K)

    # find optimal number of boosting iterations
    fit_mopt <- get_mopt(fit,
                         control,
                         y,
                         X,
                         x,
                         links,
                         loss_fun,
                         grad_fun,
                         nll_fun,
                         standardize,
                         N,
                         K)

    coef <- update_list(x, fit_mopt$coef, K)
    coef_path <- fit_mopt$coef_path
    nll_path <- fit_mopt$nll_path
    df_path <- fit_mopt$df_path
    iterations <- fit_mopt$iterations
    names(coef) <- names(coef_path) <- rep(c("location", "scale", "weight"), each = K)


    # prepare output object
    out <- list(formula = formula,
                scale.formula = scale.formula,
                weight.formula = weight.formula,
                model_frame = mf,
                coef = coef,
                family = family,
                loss = loss,
                paths = list(coef_path = coef_path,
                             nll_path = nll_path,
                             df_path = df_path),
                iterations = iterations,
                control = control,
                standardize = standardize)
    class(out) <- c("mixnhreg.boost", "mixnhreg")

    # add fit statistics
    loglik <- -sum(predict.mixnhreg(out, type = "nll", response = out$model_frame[, 1]))
    df <- sum(unlist(out$coef) != 0)
    stats <- list(n = N,
                  df = df,
                  CRPS = sum(predict.mixnhreg(out, type = "crps", response = out$model_frame[, 1])),
                  logLik = loglik,
                  AIC = -2*loglik + 2*df,
                  BIC = -2*loglik + df*log(N))
    out <- append(out,
                  list(stats = stats),
                  after = 11)
    class(out) <- c("mixnhreg.boost", "mixnhreg")

  }

  return(out)


}

#' update list entries
#' @noRd
update_list <- function(L, x, K) {

  out <- c()
  iter <- l <- 1
  u <- 0

  for (j in 1:3) {
    for (k in 1:K) {
      u <- u + length(L[[iter]])
      out[[iter]] <- x[l:u]
      l <- u+1
      iter <- iter + 1
    }
  }

  return(out)

}

# (re)standardize vectors and model matrices
#' @noRd
standardize.matrix <- function(x, center = NULL, scale = NULL, restandardize = FALSE) {
  if(is.null(center)) center <- attr(x, "standardize:center")
  if(is.null(scale)) scale <- attr(x, "standardize:scale")
  x <- as.matrix(x)
  if(is.null(center)) {
    center <- colMeans(x)
  }
  center[grep("(Intercept)", colnames(x))] <- 0
  mcenter <- matrix(rep(center, nrow(x)), nrow(x), ncol(x), byrow = TRUE)

  if(is.null(scale)) {
    scale <- sqrt(colSums((x-mcenter)^2)/nrow(x))
  }

  scale[grep("(Intercept)", colnames(x))] <- 1
  mscale  <- matrix(rep(scale , nrow(x)), nrow(x), ncol(x), byrow = TRUE)

  if(!restandardize) {
    x <- (x - mcenter)/mscale
  } else {
    x <- x*mscale + mcenter
  }
  attr(x, "standardize:center") <- center
  attr(x, "standardize:scale") <- scale
  x

}

# (re)standardize coefficients
#' @noRd
standardize.coefficients <- function(coef, center, scale, restandardize = FALSE) {

  interceptind <- grep("(Intercept)", names(coef))
  if(length(interceptind > 0)) {
    intercept <- coef[interceptind]
    center <- center[-interceptind]
    scale <- scale[-interceptind]
    coef <- coef[-interceptind]
  } else {
    intercept <- 0
    names(intercept) <- "(Intercept)"
  }
  if(restandardize) {
    intercept <- intercept - sum(coef*center/scale)
    coef <- coef/scale
  } else {
    coef <- coef*scale
    intercept <- intercept + sum(coef*center/scale)
  }
  c(intercept, coef)
}

#' @noRd
#' @importFrom parallel mclapply parLapply makeCluster stopCluster
get_mopt <- function(fit, control, y, X, x, links, loss_fun, grad_fun, nll_fun, standardize, N, K) {

  if (control$mstop == "aic") {
    aic_path <- 2*fit$nll_path + 2*rowSums(fit$coef_path != 0)
    iterations <- which.min(aic_path)
  } else if (control$mstop == "bic") {
    bic_path <- 2*fit$nll_path + rowSums(fit$coef_path != 0)*log(N)
    iterations <- which.min(bic_path)
  } else if (control$mstop == "max") {
    iterations <- control$maxit
  } else if (control$mstop == "cv") {

    # initial settings
    foldid <- control$foldid
    nfolds <- control$nfolds
    maxit <- control$maxit
    nu <- control$nu
    mc.cores <- control$cores

    if(is.null(foldid)) {
      foldid <- sample(1:nfolds, size = N, replace = TRUE)
    } else {
      nfolds <- length(unique(foldid))
    }

    if (as.character(Sys.info()["sysname"]) != "Windows") {

      # start cv
      cv_loss <- mclapply(1:nfolds, function(i) {

        # indices (rows) training and test data set
        train <- foldid != unique(foldid)[i]
        test <- !train

        y_cv_train <- y[train]
        X_cv_train <- c(lapply(1:K, function(j) {X[[j]][train, , drop = FALSE]}),
                        lapply(1:K, function(j) {X[[K+j]][train, , drop = FALSE]}),
                        lapply(1:K, function(j) {X[[2*K+j]][train, , drop = FALSE]}))
        N_train <- length(y_cv_train)

        cv_fit <- boost_noncylcic_cpp(y_cv_train,
                                      X_cv_train,
                                      x,
                                      links,
                                      loss_fun,
                                      grad_fun,
                                      nll_fun,
                                      nu,
                                      maxit,
                                      N_train,
                                      K)

        y_cv_test <- y[test]
        X_cv_test <- c(lapply(1:K, function(j) {X[[j]][test, , drop = FALSE]}),
                       lapply(1:K, function(j) {X[[K+j]][test, , drop = FALSE]}),
                       lapply(1:K, function(j) {X[[2*K+j]][test, , drop = FALSE]}))
        N_test <- length(y_cv_test)
        coef_path_cv <- cv_fit$coef_path

        # loss matrix for ith-fold
        get_loss_cv_cpp(y_cv_test,
                        X_cv_test,
                        coef_path_cv,
                        links,
                        loss_fun,
                        maxit,
                        N_test,
                        K)



      }, mc.cores = mc.cores)

    } else {

      cl <- makeCluster(mc.cores)
      # start cv
      cv_loss <- parLapply(cl, 1:nfolds, function(i) {

        # indices (rows) training and test data set
        train <- foldid != unique(foldid)[i]
        test <- !train

        y_cv_train <- y[train]
        X_cv_train <- c(lapply(1:K, function(j) {X[[j]][train, , drop = FALSE]}),
                        lapply(1:K, function(j) {X[[K+j]][train, , drop = FALSE]}),
                        lapply(1:K, function(j) {X[[2*K+j]][train, , drop = FALSE]}))
        N_train <- length(y_cv_train)

        cv_fit <- boost_noncylcic_cpp(y_cv_train,
                                      X_cv_train,
                                      x,
                                      links,
                                      loss_fun,
                                      grad_fun,
                                      nll_fun,
                                      nu,
                                      maxit,
                                      N_train,
                                      K)

        y_cv_test <- y[test]
        X_cv_test <- c(lapply(1:K, function(j) {X[[j]][test, , drop = FALSE]}),
                       lapply(1:K, function(j) {X[[K+j]][test, , drop = FALSE]}),
                       lapply(1:K, function(j) {X[[2*K+j]][test, , drop = FALSE]}))
        N_test <- length(y_cv_test)
        coef_path_cv <- cv_fit$coef_path

        # loss matrix for ith-fold
        get_loss_cv_cpp(y_cv_test,
                        X_cv_test,
                        coef_path_cv,
                        links,
                        loss_fun,
                        maxit,
                        N_test,
                        K)



      })
      stopCluster(cl)

    }
    iterations <- which.min(colSums(do.call(rbind, cv_loss)))
  }


  # adapt coefficient path and coefficients
  u <- 0
  l <- iter <- 1
  coef_path <- c()
  coef <- c()
  for (j in 1:3) {
    for (k in 1:K) {
      u <- u + ncol(X[[iter]])
      coef_path[[iter]] <- as.data.frame(fit$coef_path[, seq(l, u), drop = FALSE])
      colnames(coef_path[[iter]]) <- colnames(X[[iter]])

      # restandardize coefficient path
      z <- apply(coef_path[[iter]], 1, standardize.coefficients, center = standardize[[iter+1]]$center, scale = standardize[[iter+1]]$scale, restandardize = TRUE)
      if (!is.null(dim(z))) {
        z <- t(z)
      } else {
        z <- matrix(z, ncol = 1, dimnames = list(c(), colnames(coef_path[[iter]])))
      }

      if (j == 1) {
        z <- z*standardize[[1]]$scale
        z[, grep("(Intercept)", colnames(z))] <- z[, grep("(Intercept)", colnames(z)), drop = FALSE] + standardize[[1]]$center
      } else if (j == 2) {
        z[, grep("(Intercept)", colnames(z))] <- z[, grep("(Intercept)", colnames(z)), drop = FALSE] + as.vector(get_parameter_cpp(as.matrix(standardize[[1]]$scale), links[2], FALSE))
      }

      coef_path[[iter]] <- z
      coef_tmp <- z[iterations, ]
      coef <- c(coef, coef_tmp)
      l <- u + 1
      iter <- iter + 1
    }
  }

  out <- list(coef = coef,
              coef_path = coef_path,
              df_path = rowSums(fit$coef_path != 0),
              nll_path = fit$nll_path,
              iterations = iterations)

  return(out)


}
