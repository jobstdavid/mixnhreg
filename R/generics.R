#' Generic functions
#'
#' print model output
#' @noRd
#' @export
print.mixnhreg <- function(x, digits = max(3, getOption("digits") - 3), ...) {

  K <- length(x$coef)/3

  cat(paste("Family: ", x$family$family, " with ", K, " components", "\n", sep = ""))

  cat(paste("CRPS: ", round(x$stats$CRPS, digits),
            ", logLik: ", round(x$stats$logLik, digits),
            ", AIC: ", round(x$stats$AIC, digits),
            ", BIC: ", round(x$stats$BIC, digits),
            sep = ""), "\n")

  if(inherits(x, "mixnhreg.boost")) {

    cat(paste("Boosting iterations: ", x$iterations, ", stopping criterion: ", x$control$mstop, sep = ""))

  }

}
#' print model output
#' @noRd
#' @export
summary.mixnhreg <- function(object, digits = max(3, getOption("digits") - 3), ...) {

  K <- length(object$coef)/3

  if(inherits(object, "mixnhreg.optim")) {

    for (k in 1:K) {
      cat(paste("Coefficients component ", k, " (location model with ", as.character(object$family$links[1]), "-link):\n", sep = ""))
      print.default(format(object$coef[[k]], digits = digits), print.gap = 2, quote = FALSE)
      cat("\n")
    }

    for (k in 1:K) {
      cat(paste("Coefficients component ", k, " (scale model with ", as.character(object$family$links[2]), "-link):\n", sep = ""))
      print.default(format(object$coef[[K+k]], digits = digits), print.gap = 2, quote = FALSE)
      cat("\n")
    }

    for (k in 1:K) {
      cat(paste("Coefficients component ", k, " (weight model with ", as.character(object$family$links[3]), "-link):\n", sep = ""))
      print.default(format(object$coef[[2*K+k]], digits = digits), print.gap = 2, quote = FALSE)
      cat("\n")
    }

  } else if (inherits(object, "mixnhreg.boost")) {

    for (k in 1:K) {
      coefs <- object$coef[[k]]
      coefs <- coefs[coefs != 0]
      if (length(coefs) != 0) {
        cat(paste("Coefficients component ", k, " (location model with ", as.character(object$family$links[1]), "-link):\n", sep = ""))
        print.default(format(coefs, digits = digits), print.gap = 2, quote = FALSE)
        cat("\n")
      } else {
        cat(paste("Coefficients component ", k, " (location model with ", as.character(object$family$links[1]), "-link):\n", sep = ""))
        cat("No non-zero coefficients!")
        cat("\n\n")
      }
    }

    for (k in 1:K) {
      coefs <- object$coef[[K+k]]
      coefs <- coefs[coefs != 0]
      if (length(coefs) != 0) {
        cat(paste("Coefficients component ", k, " (scale model with ", as.character(object$family$links[2]), "-link):\n", sep = ""))
        print.default(format(coefs, digits = digits), print.gap = 2, quote = FALSE)
        cat("\n")
      } else {
        cat(paste("Coefficients component ", k, " (scale model with ", as.character(object$family$links[2]), "-link):\n", sep = ""))
        cat("No non-zero coefficients!")
        cat("\n\n")
      }
    }

    for (k in 1:K) {
      coefs <- object$coef[[2*K+k]]
      coefs <- coefs[coefs != 0]
      if (length(coefs) != 0) {
        cat(paste("Coefficients component ", k, " (weight model with ", as.character(object$family$links[3]), "-link):\n", sep = ""))
        print.default(format(coefs, digits = digits), print.gap = 2, quote = FALSE)
        cat("\n")
      } else {
        cat(paste("Coefficients component ", k, " (weight model with ", as.character(object$family$links[3]), "-link):\n", sep = ""))
        cat("No non-zero coefficients!")
        cat("\n\n")
      }
    }

  }

  print.mixnhreg(object, digits = digits)

}
#' print model log-likelihood
#' @noRd
#' @export
logLik.mixnhreg <- function(object, ...) {

  structure(object$stats$logLik, df = object$stats$df, class = "logLik")

}
#' print AIC
#' @noRd
#' @export
AIC.mixnhreg <- function(object, ...) {

  object$stats$AIC

}
#' print BIC
#' @noRd
#' @importFrom stats BIC
#' @export
BIC.mixnhreg <- function(object, ...) {

  object$stats$BIC

}
