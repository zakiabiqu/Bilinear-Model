#' Forecast from a bilinear model
#'
#' Computes \eqn{E(X_{N+h} | X_1, \dots, X_N)} for \eqn{h = 1, \dots,} `n.ahead`.
#' The one-step forecast is exact. For \eqn{h \ge 2}, the conditional expectation
#' involves products of future values and future errors, so it is evaluated by
#' Monte Carlo: `nsim` future paths are simulated from the fitted model and
#' averaged. Prediction intervals are the empirical quantiles of the simulated
#' paths, because forecast distributions of bilinear models are usually not
#' symmetric.
#'
#' @param object A `"bilinear"` object.
#' @param n.ahead Forecast horizon.
#' @param nsim Number of simulated paths for \eqn{h \ge 2}.
#' @param innov Source of future errors: `"normal"` draws from
#'   \eqn{N(0, \hat\sigma^2)}; `"bootstrap"` resamples the residuals (useful when
#'   residual normality is rejected).
#' @param level Confidence levels (in percent) of the prediction intervals.
#' @param newx Optional out-of-sample observations following the data. If given,
#'   rolling one-step forecasts are computed for them with the parameters held
#'   fixed, and returned in `onestep`.
#' @param seed Optional random seed for reproducible Monte Carlo results.
#' @param x A `"bilinear_forecast"` object.
#' @param digits Number of decimal places.
#' @param ... Not used.
#'
#' @return An object of class `"bilinear_forecast"`, a list with `pred`, `se`,
#'   `lower`, `upper` (one column per level), `level`, `n.ahead`, `nsim`,
#'   `innov` and `onestep` (a table of actual values, forecasts and errors, and
#'   their mean squared error, or `NULL`).
#'
#' @examples
#' sun <- window(sunspot.year, 1700, 1920)
#' fit <- bilinear(sun, ar = c(1, 2, 9), gamma = 10,
#'                 bl = list(c(2, 1), c(8, 1), c(1, 3)))
#' predict(fit, n.ahead = 5, nsim = 2000, seed = 1)
#'
#' fc <- predict(fit, newx = window(sunspot.year, 1921, 1955))
#' fc$onestep$mse
#' @export
predict.bilinear <- function(object, n.ahead = 1, nsim = 10000, innov = c("normal", "bootstrap"),
                             level = c(80, 95), newx = NULL, seed = NULL, ...) {
  innov <- match.arg(innov)
  if (!is.null(seed)) set.seed(seed)
  m <- object$model
  th <- object$coef
  pa <- length(m$ar); pc <- length(m$ma); pb <- nrow(m$bl)
  a <- th[seq_len(pa)]
  cc <- th[pa + seq_len(pc)]
  b <- th[pa + pc + seq_len(pb)]
  al <- if (m$intercept) th[length(th)] else 0
  x <- as.numeric(object$x)
  N <- length(x)
  e <- as.numeric(object$residuals)
  e[is.na(e)] <- 0
  # conditional mean of X_t given the past
  mu_at <- function(xx, ee, t) {
    v <- al
    if (pa) v <- v - sum(a * xx[t - m$ar])
    if (pc) v <- v + sum(cc * ee[t - m$ma])
    if (pb) v <- v + sum(b * xx[t - m$bl[, 1]] * ee[t - m$bl[, 2]])
    v
  }
  # simulate nsim future paths
  h <- n.ahead
  X <- matrix(rep(c(x, rep(NA_real_, h)), each = nsim), nrow = nsim)
  E <- matrix(rep(c(e, rep(NA_real_, h)), each = nsim), nrow = nsim)
  pool <- as.numeric(stats::na.omit(as.numeric(object$residuals)))
  for (k in seq_len(h)) {
    t <- N + k
    eps <- if (innov == "normal") stats::rnorm(nsim, 0, sqrt(object$sigma2)) else sample(pool, nsim, replace = TRUE)
    v <- rep(al, nsim)
    if (pa) for (i in seq_len(pa)) v <- v - a[i] * X[, t - m$ar[i]]
    if (pc) for (i in seq_len(pc)) v <- v + cc[i] * E[, t - m$ma[i]]
    if (pb) for (i in seq_len(pb)) v <- v + b[i] * X[, t - m$bl[i, 1]] * E[, t - m$bl[i, 2]]
    X[, t] <- v + eps
    E[, t] <- eps
  }
  sims <- X[, N + seq_len(h), drop = FALSE]
  pred <- colMeans(sims)
  pred[1] <- mu_at(x, e, N + 1)
  se <- apply(sims, 2, stats::sd)
  se[1] <- if (innov == "normal") sqrt(object$sigma2) else stats::sd(pool)
  probs <- sort(c((1 - level / 100) / 2, 1 - (1 - level / 100) / 2))
  qs <- apply(sims, 2, stats::quantile, probs = probs)
  qs <- matrix(qs, ncol = h)
  nl <- length(level)
  lower <- t(qs[rev(seq_len(nl)), , drop = FALSE])
  upper <- t(qs[nl + seq_len(nl), , drop = FALSE])
  colnames(lower) <- colnames(upper) <- paste0(level, "%")
  xo <- object$x
  if (stats::is.ts(xo)) {
    st <- stats::tsp(xo)[2] + 1 / stats::frequency(xo)
    pred <- stats::ts(pred, start = st, frequency = stats::frequency(xo))
    se <- stats::ts(se, start = st, frequency = stats::frequency(xo))
    lower <- stats::ts(lower, start = st, frequency = stats::frequency(xo))
    upper <- stats::ts(upper, start = st, frequency = stats::frequency(xo))
  }
  # rolling one-step forecasts on new data
  onestep <- NULL
  if (!is.null(newx)) {
    xn <- c(x, as.numeric(newx))
    en <- c(e, rep(0, length(newx)))
    pr <- numeric(length(newx))
    for (k in seq_along(newx)) {
      t <- N + k
      pr[k] <- mu_at(xn, en, t)
      en[t] <- xn[t] - pr[k]
    }
    tab <- data.frame(actual = as.numeric(newx), forecast = pr, error = as.numeric(newx) - pr)
    onestep <- list(table = tab, mse = mean(tab$error^2))
  }
  structure(list(pred = pred, se = se, lower = lower, upper = upper, level = level,
                 n.ahead = h, nsim = nsim, innov = innov, onestep = onestep),
            class = "bilinear_forecast")
}

#' @rdname predict.bilinear
#' @export
print.bilinear_forecast <- function(x, digits = 4, ...) {
  tab <- cbind(Forecast = as.numeric(x$pred), SE = as.numeric(x$se))
  for (j in seq_along(x$level)) {
    tab <- cbind(tab, as.numeric(x$lower[, j]), as.numeric(x$upper[, j]))
    colnames(tab)[ncol(tab) - 1:0] <- paste(c("Lo", "Hi"), x$level[j])
  }
  rownames(tab) <- paste0("h=", seq_len(x$n.ahead))
  cat(sprintf("\nForecasts (h = 1 exact; h >= 2 Monte Carlo, %d paths, %s innovations)\n\n", x$nsim, x$innov))
  print(round(tab, digits))
  if (!is.null(x$onestep)) cat(sprintf("\nOne-step forecasts on newx: %d values, MSE = %s\n",
                                       nrow(x$onestep$table), format(signif(x$onestep$mse, 6))))
  cat("\n")
  invisible(x)
}
