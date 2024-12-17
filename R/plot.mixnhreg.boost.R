#' Plot Standardized Coefficient and log-likelihood Paths
#' @description This function plots standardized coefficient and log-likelihood paths of an \code{mixnhreg.boost} object.
#'
#' @param x an object of class \code{mixnhreg.boost}.
#' @param type character; specifies for which distribution parameters, i.e. "\code{location}" (default), "\code{scale}", "\code{weight}" the paths should be created.
#' @param loglik logical; if \code{TRUE}, the log-likelihood contribution path is plotted instead of the standardized coefficient path. Default: \code{FALSE}.
#' @param mstop logical; if \code{TRUE} (default) a vertical line for the optimal number of boosting iterations is plotted, otherwise \code{FALSE} not.
#' @param ... unused.
#'
#'
#' @examples
#' # load data
#' data("station")
#'
#' # fit mixture normal distribution with two components via gradient-boosting
#' (fit_boost <- mixnhreg(formula = obs ~ sin1 + cos1 + temp_mean | temp_ctrl,
#'                        scale.formula = ~ sin1 + cos1 + log(temp_sd) | 1,
#'                        weight.formula = ~ sin1 + cos1 | 1,
#'                        data = station,
#'                        control = control_boost(mstop = "cv")))
#'
#' # standardized coefficient path for location parameter
#' plot(fit_boost, type = "location")
#'
#' # log-likelihood contribution path for location parameter
#' plot(fit_boost, type = "location", loglik = TRUE)
#'
#' @rdname plot.mixnhreg.boost
#' @importFrom ggplot2 ggplot geom_step aes labs scale_y_continuous scale_x_continuous dup_axis geom_vline guide_axis theme theme_bw
#' @importFrom scales hue_pal
#' @export
plot.mixnhreg.boost <- function(x, type = "location", loglik = FALSE, mstop = TRUE, ...) {


  y <- variable <- component <- object <- NULL # for CRAN checks

  # this function depends on the following R-packages:
  # ggplot2, (scales, ggh4x)

  if (x$iterations == 1) {
    stop("Number of boosting iterations is to small!")
  }

  # number of components
  f <- Formula(x$formula)
  K <- length(f)[2]
  maxit <- x$control$maxit+1
  col <- hue_pal()(K)

  # coefficient path
  path <- x$paths$coef_path
  nll_path <- x$paths$nll_path
  standardize <- x$standardize
  links <- x$family$links

  # adapt coefficient path and coefficients
  iter <- 1
  for (j in 1:3) {
    for (k in 1:K) {

      # restandardize coefficient path
      z <- t(apply(path[[iter]], 1, standardize.coefficients, center = standardize[[iter+1]]$center, scale = standardize[[iter+1]]$scale, restandardize = FALSE))
      if (j == 1) {
        z[, grep("(Intercept)", colnames(z))] <- z[, grep("(Intercept)", colnames(z)), drop = FALSE] - standardize[[1]]$center
        z <- z/standardize[[1]]$scale
      } else if (j == 2) {
        z[, grep("(Intercept)", colnames(z))] <- z[, grep("(Intercept)", colnames(z)), drop = FALSE] - as.vector(get_parameter_cpp(as.matrix(standardize[[1]]$scale), links[2], FALSE))
      }

      path[[iter]] <- z
      iter <- iter + 1
    }
  }

  if (type == "location") {
    path <- do.call("rbind", lapply(1:K, function(k) {
      path_k <- path[[k]]
      if (!loglik) {
        index <- as.numeric(which(apply(path_k, 2, var) != 0))
        path_k <- path_k[, index, drop = F]
      } else {
        path_k <- as.data.frame(apply(path_k, 2, diff) != 0)
        path_k[["(Intercept)"]][rowSums(path_k) == 2] <- FALSE
        path_k <- rbind(0, diff(-nll_path)*path_k)
        path_k <- apply(path_k, 2, cumsum)
      }
      data.frame(x = rep(1:maxit, ncol(path_k)),
                 y = as.vector(path_k),
                 component = factor(k, levels = 1:K),
                 variable = rep(paste0(colnames(path_k), ".", k), each = maxit),
                 label = rep(colnames(path_k), each = maxit))
    }))
    title <- "Location parameter"
  } else if (type == "scale") {
    path <- do.call("rbind", lapply(1:K, function(k) {
      path_k <- path[[K+k]]
      if (!loglik) {
        index <- as.numeric(which(apply(path_k, 2, var) != 0))
        path_k <- path_k[, index, drop = F]
      } else {
        path_k <- as.data.frame(apply(path_k, 2, diff) != 0)
        path_k[["(Intercept)"]][rowSums(path_k) == 2] <- FALSE
        path_k <- rbind(0, diff(-nll_path)*path_k)
        path_k <- apply(path_k, 2, cumsum)
      }
      data.frame(x = rep(1:maxit, ncol(path_k)),
                 y = as.vector(path_k),
                 component = factor(k, levels = 1:K),
                 variable = rep(paste0(colnames(path_k), ".", k), each = maxit),
                 label = rep(colnames(path_k), each = maxit))
    }))
    title <- "Scale parameter"
  } else if (type == "weight") {
    path <- do.call("rbind", lapply(1:K, function(k) {
      path_k <- path[[2*K+k]]
      if (!loglik) {
        index <- as.numeric(which(apply(path_k, 2, var) != 0))
        path_k <- path_k[, index, drop = F]
      } else {
        path_k <- as.data.frame(apply(path_k, 2, diff) != 0)
        path_k[["(Intercept)"]][rowSums(path_k) == 2] <- FALSE
        path_k <- rbind(0, diff(-nll_path)*path_k)
        path_k <- apply(path_k, 2, cumsum)
      }
      data.frame(x = rep(1:maxit, ncol(path_k)),
                 y = as.vector(path_k),
                 component = factor(k, levels = 1:K),
                 variable = rep(paste0(colnames(path_k), ".", k), each = maxit),
                 label = rep(colnames(path_k), each = maxit))
    }))
    title <- "Weight parameter"
  } else {
    stop("This 'type' is not implemented!")
  }

  p <- ggplot(data = path) +
    geom_step(aes(x = x, y = y, group = variable, color = component), linewidth = 1) +
    labs(x = "Boosting iterations", y = ifelse(!loglik, "Standardized coefficient", "log-likelihood contribution"), color = "Component", title = title) +
    scale_y_continuous(sec.axis = dup_axis(
      name = "",
      breaks = subset(path, x == maxit)$y,
      labels = subset(path, x == maxit)$label,
      guide = guide_axis(check.overlap = T))
    ) +
    theme_bw() +
    theme(legend.position = "bottom")

  if (mstop) {
    p <- p + geom_vline(xintercept = x$iterations, linetype = "dashed", linewidth = 1) +
      scale_x_continuous(sec.axis = dup_axis(
        name = "",
        breaks = x$iterations,
        labels = x$control$mstop))
  }

  return(p)

}


