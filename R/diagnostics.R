#' Residual diagnostics for a bilinear model
#'
#' Computes the tests reported by [summary.bilinear()]:
#' \itemize{
#'   \item the independence test \eqn{W} on the squared residuals,
#'     \eqn{H_0: \rho_{e^2}(1) = 0}, approximately \eqn{N(0,1)};
#'   \item Lomnicki's skewness and kurtosis tests \eqn{u_1} (\eqn{H_0: \gamma_1 = 0})
#'     and \eqn{u_2} (\eqn{H_0: \gamma_2 = 0}), each approximately \eqn{N(0,1)},
#'     using \eqn{Var(G_1) = (6/N)\sum\rho^2} and \eqn{Var(G_2) = (24/N)\sum\rho^4};
#'   \item the joint normality statistic \eqn{u_1^2 + u_2^2}, approximately \eqn{\chi^2_2};
#'   \item the Ljung-Box test on the residuals.
#' }
#'
#' @param object A `"bilinear"` object.
#' @param lag Number of lags for the Ljung-Box test and the autocorrelation sums.
#' @return A data frame with columns `Statistic`, `Distribution` and `p.value`.
#' @examples
#' fit <- bilinear(log10(lynx), ar = c(1, 2), bl = list(c(1, 1)))
#' bl_diagnostics(fit)
#' @export
bl_diagnostics <- function(object, lag = 10) {
  e <- as.numeric(stats::na.omit(as.numeric(object$residuals)))
  N <- length(e)
  # independence test W on squared residuals
  y <- e^2
  yb <- mean(y)
  R1 <- sum((y[-N] - yb) * (y[-1] - yb)) / N
  R011 <- sum((y[-N] - yb)^2 * (y[-1] - yb)^2) / N
  W <- sqrt(N) * R1 / sqrt(R011)
  # Lomnicki skewness and kurtosis tests
  eb <- mean(e)
  m2 <- mean((e - eb)^2); m3 <- mean((e - eb)^3); m4 <- mean((e - eb)^4)
  G1 <- m3 / m2^1.5
  G2 <- m4 / m2^2 - 3
  L <- min(lag, N - 1)
  rho <- sapply(seq_len(L), function(s) sum((e[1:(N - s)] - eb) * (e[(1 + s):N] - eb)) / N / m2)
  u1 <- G1 / sqrt(6 / N * (1 + 2 * sum(rho^2)))
  u2 <- G2 / sqrt(24 / N * (1 + 2 * sum(rho^4)))
  # Ljung-Box on the residuals
  fitdf <- length(object$model$ar) + length(object$model$ma)
  lb_lag <- max(lag, fitdf + 1)
  lb <- stats::Box.test(e, lag = lb_lag, type = "Ljung-Box", fitdf = fitdf)
  data.frame(
    Statistic = c(W, u1, u2, u1^2 + u2^2, unname(lb$statistic)),
    Distribution = c("N(0,1)", "N(0,1)", "N(0,1)", "Chi-sq(2)", sprintf("Chi-sq(%d)", lb_lag - fitdf)),
    p.value = c(2 * stats::pnorm(-abs(c(W, u1, u2))),
                stats::pchisq(u1^2 + u2^2, 2, lower.tail = FALSE), lb$p.value),
    row.names = c("Independence W", sprintf("Skewness u1 (G1 = %.3f)", G1),
                  sprintf("Kurtosis u2 (G2 = %.3f)", G2), "Normality u1^2 + u2^2",
                  sprintf("Ljung-Box on e (lag %d)", lb_lag)),
    stringsAsFactors = FALSE)
}
