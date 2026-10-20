#' Extract coefficients and covariance of a bilinear model
#'
#' @param object A `"bilinear"` object.
#' @param ... Not used.
#' @return `coef()` returns the named vector of estimates; `vcov()` returns the
#'   approximate covariance matrix \eqn{\hat\sigma^2 (J'J)^{-1}}.
#' @name bilinear-accessors
NULL

#' @rdname bilinear-accessors
#' @export
coef.bilinear <- function(object, ...) object$coef

#' @rdname bilinear-accessors
#' @export
vcov.bilinear <- function(object, ...) object$vcov

# Model equation as text, used by print.summary.bilinear()
bl_equation <- function(object) {
  m <- object$model
  lhs <- paste(c("X[t]", if (length(m$ar)) paste0("a", m$ar, " X[t-", m$ar, "]")), collapse = " + ")
  rhs <- c(if (m$intercept) "alpha",
           if (length(m$ma)) paste0("c", m$ma, " e[t-", m$ma, "]"),
           if (nrow(m$bl)) paste0("b(", m$bl[, 1], ",", m$bl[, 2], ") X[t-", m$bl[, 1], "] e[t-", m$bl[, 2], "]"),
           "e[t]")
  paste(lhs, "=", paste(rhs, collapse = " + "))
}

#' Print a bilinear model
#'
#' @param x A `"bilinear"` object.
#' @param digits Number of significant digits.
#' @param ... Not used.
#' @return `x`, invisibly.
#' @export
print.bilinear <- function(x, digits = 6, ...) {
  cat("\nCall:\n")
  print(x$call)
  cat("\nCoefficients:\n")
  tab <- cbind(Estimate = x$coef, `Std. Error` = x$se)
  print(signif(tab, digits))
  cat(sprintf("\nsigma^2 = %s   log-likelihood = %s   AIC = %s   NAIC = %s\n",
              format(signif(x$sigma2, digits)), format(round(x$loglik, 2)),
              format(round(x$aic, 2)), format(round(x$naic, 4))))
  cat(sprintf("n = %d (N = %d, gamma = %d)   converged: %s (%d iterations)   invertible: %s\n\n",
              x$n, x$N, x$gamma, if (x$converged) "yes" else "no", x$iter,
              if (x$invertibility$invertible) "yes" else "no"))
  invisible(x)
}

#' Summarise a bilinear model
#'
#' Coefficient table with z tests, fit statistics and residual diagnostics
#' (see [bl_diagnostics()]).
#'
#' @param object A `"bilinear"` object.
#' @param lag Number of lags used in the Ljung-Box test and in the variance of
#'   the skewness and kurtosis statistics.
#' @param x A `"summary.bilinear"` object.
#' @param digits Number of significant digits.
#' @param ... Not used.
#' @return An object of class `"summary.bilinear"`.
#' @export
summary.bilinear <- function(object, lag = 10, ...) {
  tab <- cbind(Estimate = object$coef, `Std. Error` = object$se,
               `z value` = object$zvalue, `Pr(>|z|)` = object$pvalue)
  structure(list(object = object, coefficients = tab,
                 residual_quantiles = stats::quantile(object$residuals, na.rm = TRUE),
                 diagnostics = bl_diagnostics(object, lag = lag)),
            class = "summary.bilinear")
}

#' @rdname summary.bilinear
#' @export
print.summary.bilinear <- function(x, digits = 6, ...) {
  o <- x$object
  cat("\nBilinear time series model (conditional least squares)\n\nCall:\n")
  print(o$call)
  cat("\nModel:\n", bl_equation(o), "\n", sep = "")
  cat("\nResiduals:\n")
  rq <- x$residual_quantiles
  names(rq) <- c("Min", "1Q", "Median", "3Q", "Max")
  print(signif(rq, 4))
  cat("\nCoefficients:\n")
  stats::printCoefmat(x$coefficients, digits = digits, signif.stars = TRUE,
                      P.values = TRUE, has.Pvalue = TRUE)
  if (nrow(o$model$bl)) cat("Note: standard errors of b may be underestimated; consider bootstrap.\n")
  cat(sprintf("\nsigma^2 = %s   log-likelihood = %s   AIC = %s   NAIC = %s\n",
              format(signif(o$sigma2, digits)), format(round(o$loglik, 2)),
              format(round(o$aic, 2)), format(round(o$naic, 4))))
  cat(sprintf("n = %d (N = %d, gamma = %d), %d parameters\n", o$n, o$N, o$gamma, o$npar))
  cat("\nResidual diagnostics:\n")
  d <- x$diagnostics
  out <- data.frame(Statistic = formatC(d$Statistic, digits = 4, format = "f"),
                    Distribution = d$Distribution,
                    `p-value` = format.pval(d$p.value, digits = 3, eps = 1e-4),
                    row.names = rownames(d), check.names = FALSE)
  print(out)
  cat(sprintf("\nConvergence: %s (%d iterations)\nInvertibility: %s\n\n",
              if (o$converged) "yes" else "no", o$iter,
              if (o$invertibility$invertible) "yes" else "no"))
  invisible(x)
}
