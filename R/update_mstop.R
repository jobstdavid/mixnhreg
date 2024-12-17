#' Auxiliary Function for Updating the Stopping Criterion in the Gradient-Boosting Estimation
#' @description Auxiliary function for updating the stopping criterion "\code{mstop}" and consequently an \code{mixnhreg.boost} object
#' based on the initial settings provided for the "\code{control}" argument in \code{\link{mixnhreg}}.
#'
#' @param object an object of class \code{mixnhreg.boost}.
#' @param mstop character; new stopping criterion "\code{max}", "\code{aic}" (default), "\code{bic}", "\code{cv}".
#'
#' @return An object of class \code{mixnhreg}.
#'
#' @examples
#' # load data
#' data("station")
#'
#' # fit mixture normal distribution with two components via gradient-boosting
#' (fit_boost1 <- mixnhreg(formula = obs ~ sin1 + cos1 + temp_mean | temp_ctrl,
#'                         scale.formula = ~ sin1 + cos1 + log(temp_sd) | 1,
#'                         weight.formula = ~ sin1 + cos1 | 1,
#'                         data = station,
#'                         control = control_boost(mstop = "cv")))
#'
#' # update stopping criterion to "bic"
#' (fit_boost2 <- update_mstop(fit_boost1, mstop = "bic"))
#'
#' @rdname update_mstop
#' @export
update_mstop <- function(object, mstop = "aic") {

  f <- Formula(object$formula)
  K <- length(f)[2]
  loss <- object$loss

  # paths
  path <- object$paths$coef_path
  nll_path <- object$paths$nll_path
  df_path <- object$paths$df_path


  if (mstop == "aic") {
    aic_path <- 2*nll_path + 2*df_path
    iterations <- which.min(aic_path)
  } else if (mstop == "bic") {
    bic_path <- 2*nll_path + df_path*log(object$stats$n)
    iterations <- which.min(bic_path)
  } else if (mstop == "max") {
    iterations <- object$control$maxit
  } else if (mstop == "cv") {

    # initial settings
    mf <- object$model_frame
    y <- mf[, 1]
    X <- c(lapply(1:K, function(j) {
      model.matrix(formula(as.Formula(object$formula), lhs = F, rhs = j), mf)}),
      lapply(1:K, function(j) {
        model.matrix(formula(as.Formula(object$scale.formula), lhs = F, rhs = j), mf)}),
      lapply(1:K, function(j) {
        model.matrix(formula(as.Formula(object$weight.formula), lhs = F, rhs = j), mf)}))

    # standardization
    standardize <- object$standardize
    # standardize predictor matrices
    X <- lapply(1:length(X), function(k) standardize.matrix(X[[k]], center = standardize[[k+1]]$center, scale = standardize[[k+1]]$scale))
    # standardize response
    y <- standardize.matrix(y, center = standardize[[1]]$center, scale = standardize[[1]]$scale)


    # coefficients
    x <- c(lapply(1:K, function(j) {rep(0, ncol(X[[j]]))}),
           lapply(1:K, function(j) {rep(0, ncol(X[[K+j]]))}),
           lapply(1:K, function(j) {rep(0, ncol(X[[2*K+j]]))}))
    N <- object$stats$n

    # get link functions
    links <- as.character(object$family$links)

    # implemented family discrimination for faster estimation
    if (object$family$dist == "mixnorm") {

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
                         "crps" = object$family$crps,
                         "nll" = object$family$nll)
      # get gradient of loss function
      grad_fun <- switch(loss,
                         "crps" = object$family$grad_crps,
                         "nll" = object$family$grad_nll)

      # get negative logLik function
      nll_fun <- object$family$nll

    }

    foldid <- object$control$foldid
    nfolds <- object$control$nfolds
    maxit <- object$control$maxit
    nu <- object$control$nu
    mc.cores <- object$control$cores

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
        X_cv_train <- c(lapply(1:K, function(j) {X[[j]][train, ]}),
                        lapply(1:K, function(j) {X[[K+j]][train, ]}),
                        lapply(1:K, function(j) {X[[2*K+j]][train, ]}))
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
        X_cv_test <- c(lapply(1:K, function(j) {X[[j]][test, ]}),
                       lapply(1:K, function(j) {X[[K+j]][test, ]}),
                       lapply(1:K, function(j) {X[[2*K+j]][test, ]}))
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
        X_cv_train <- c(lapply(1:K, function(j) {X[[j]][train, ]}),
                        lapply(1:K, function(j) {X[[K+j]][train, ]}),
                        lapply(1:K, function(j) {X[[2*K+j]][train, ]}))
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
        X_cv_test <- c(lapply(1:K, function(j) {X[[j]][test, ]}),
                       lapply(1:K, function(j) {X[[K+j]][test, ]}),
                       lapply(1:K, function(j) {X[[2*K+j]][test, ]}))
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
      stopCluster(cl)

    }

    iterations <- which.min(colSums(do.call(rbind, cv_loss)))

  }

  # update coefficients
  coef <- lapply(1:(3*K), function(k) object$paths$coef_path[[k]][iterations, ])
  names(coef) <- rep(c("location", "scale", "weight"), each = K)
  object$coef <- coef

  # update fit statistics
  object$stats$logLik <- -sum(predict.mixnhreg(object, type = "nll", response = object$model_frame[, 1]))
  object$stats$CRPS <- sum(predict.mixnhreg(object, type = "crps", response = object$model_frame[, 1]))
  object$stats$df <- sum(unlist(object$coef) != 0)
  object$stats$AIC <- -2*object$stats$logLik + 2*object$stats$df
  object$stats$BIC <- -2*object$stats$logLik + object$stats$df*log(object$stats$n)

  # update mstop and iterations
  object$iterations <- iterations
  object$control$mstop <- mstop

  return(object)

}
