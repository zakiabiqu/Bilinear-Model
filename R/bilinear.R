#' Fit a subset bilinear time series model
#'
#' Estimates the parameters of a (subset) bilinear time series model
#' \deqn{X_t + \sum_{k \in K} a_k X_{t-k} = \alpha + \sum_{l \in C} c_l e_{t-l}
#'       + \sum_{(i,j) \in B} b_{ij} X_{t-i} e_{t-j} + e_t}{X[t] + sum a[k] X[t-k] = alpha + sum c[l] e[t-l] + sum b[i,j] X[t-i] e[t-j] + e[t]}
#' by conditional least squares (equivalently, conditional Gaussian maximum
#' likelihood), following Subba Rao and Gabr (1984).
#'
#' The residuals \eqn{e_t} are rebuilt recursively with \eqn{e_t = 0} for
#' \eqn{t \le \gamma}, and the sum of squares \eqn{Q = \sum_{t > \gamma} e_t^2} is
#' minimised with Gauss-Newton iterations. The derivatives
#' \eqn{\partial e_t / \partial \theta} are computed exactly by recursion, and every
#' step is protected by step-halving.
#'
#' When `start` is not given, starting values are built in stages: the AR part
#' and the intercept come from ordinary least squares, then MA terms and
#' bilinear terms are added one at a time, in the order given, each with a small
#' grid of starting values. A grid point other than zero is only accepted if it
#' lowers \eqn{Q} *and* gives an invertible model. Because of this, the order of
#' the pairs in `bl` can change the local minimum that is reached.
#'
#' Note the sign convention: AR coefficients sit on the left-hand side, so a
#' model written as \eqn{X_t = 0.4 X_{t-1} + \dots} has `a1 = -0.4`.
#'
#' @param x A numeric vector or `ts` object without missing values.
#' @param ar Integer vector of AR lags (the set \eqn{K}), e.g. `c(1, 2, 9)`.
#'   `NULL` for no AR part.
#' @param ma Integer vector of MA lags (the set \eqn{C}). `NULL` for no MA part.
#' @param bl A list of lag pairs `c(i, j)`, one per bilinear term
#'   \eqn{b_{ij} X_{t-i} e_{t-j}}, e.g. `list(c(2, 1), c(8, 1))`. A single pair
#'   may be given as `c(i, j)`. `NULL` for no bilinear part.
#' @param intercept Logical; include the constant \eqn{\alpha}?
#' @param gamma Number of initial observations used as fixed starting values.
#'   Defaults to the largest lag in the model. Use the same value for every
#'   model you want to compare by AIC, so that all of them use the same sample.
#' @param start Optional numeric vector of starting values, in the order
#'   `a`, `c`, `b`, `alpha`, on the original scale of `x`.
#' @param b_grid Starting values tried for each new bilinear coefficient, in
#'   units of `1 / sd(x)`.
#' @param c_grid Starting values tried for each new MA coefficient.
#' @param scale Divide `x` by this constant before estimation (useful when
#'   `x` is large, which makes `b` very small). Results are returned on the
#'   original scale.
#' @param tol Convergence tolerance on the relative change of \eqn{Q}.
#' @param maxit Maximum number of Gauss-Newton iterations.
#' @param trace Logical; print \eqn{Q} at every iteration of the final fit.
#'
#' @return A list of class `"bilinear"` with components:
#' \item{coef}{a named vector of the estimated AR, MA and bilinear coefficients
#'   and the intercept, which can be extracted by the [coef()] method.}
#' \item{se}{the approximate standard errors of `coef`.}
#' \item{zvalue}{the z statistics, `coef / se`.}
#' \item{pvalue}{the two-sided p-values of the z tests of \eqn{H_0: \theta = 0}.}
#' \item{vcov}{the approximate variance matrix of `coef`,
#'   \eqn{\hat\sigma^2 (J'J)^{-1}}, which can be extracted by the [vcov()] method.}
#' \item{sigma2}{the conditional MLE of the innovations variance, `rss / n`.}
#' \item{rss}{the residual sum of squares \eqn{Q}.}
#' \item{loglik}{the maximized conditional Gaussian log-likelihood.}
#' \item{aic}{the AIC value, \eqn{n \log\hat\sigma^2 + 2q}.}
#' \item{naic}{the normalized AIC, `aic / n`.}
#' \item{n}{the number of residuals used in the fitting, `N - gamma`.}
#' \item{N}{the number of observations in `x`.}
#' \item{gamma}{the number of initial observations not used in the fitting.}
#' \item{npar}{the number of estimated parameters \eqn{q} (including the intercept).}
#' \item{residuals}{the fitted innovations, `NA` for the first `gamma` observations.}
#' \item{fitted.values}{the one-step fitted values `x - residuals`.}
#' \item{converged}{logical: did the Gauss-Newton iteration converge?}
#' \item{iter}{the number of iterations used.}
#' \item{history}{a matrix with the parameter values and \eqn{Q} at each iteration.}
#' \item{invertibility}{a list with `log10delta`, the log10 size at the end of
#'   the series of a unit perturbation of the starting residuals, and
#'   `invertible`, which is `TRUE` when `log10delta < 0`.}
#' \item{jacobian}{the matrix \eqn{J} of derivatives \eqn{\partial e_t / \partial\theta}
#'   at the estimates.}
#' \item{model}{a list giving the model specification: `ar`, `ma`, `bl`,
#'   `intercept` and `scale`.}
#' \item{x}{the series used in the fitting.}
#' \item{call}{the matched call.}
#'
#' @references
#' Subba Rao, T. and Gabr, M. M. (1984). *An Introduction to Bispectral
#' Analysis and Bilinear Time Series Models*. Lecture Notes in Statistics 24.
#' Springer. \doi{10.1007/978-1-4684-6318-7}
#'
#' @seealso [summary.bilinear()], [predict.bilinear()]
#'
#' @examples
#' sun <- window(sunspot.year, 1700, 1920)
#' fit <- bilinear(sun, ar = c(1, 2, 9), gamma = 10,
#'                 bl = list(c(2, 1), c(8, 1), c(1, 3), c(4, 3),
#'                           c(1, 6), c(2, 4), c(3, 2)))
#' fit
#' summary(fit)
#'
#' x <- c(-0.6, -0.2, -1.2, -0.5, 0.0, 0.1, 1.6, 1.9, 0.8, -0.5, 0.4)
#' small <- bilinear(x, ar = 1, bl = list(c(1, 1)), intercept = FALSE,
#'                   start = c(-0.5546, 0), trace = TRUE)
#' small$history
#' @export
bilinear <- function(x, ar = NULL, ma = NULL, bl = NULL, intercept = TRUE, gamma = NULL,
                     start = NULL, b_grid = c(-0.4, -0.2, 0, 0.2, 0.4), c_grid = c(-0.5, 0, 0.5),
                     scale = 1, tol = 1e-8, maxit = 100, trace = FALSE) {
  cl <- match.call()
  xo <- x
  xv <- as.numeric(x)
  if (any(!is.finite(xv))) stop("x must not contain missing or infinite values")

  # ---- model specification ------------------------------------------------
  ar <- sort(unique(as.integer(ar)))
  ma <- sort(unique(as.integer(ma)))
  if (is.null(bl) || length(bl) == 0) {
    blm <- matrix(integer(0), ncol = 2)
  } else {
    if (!is.list(bl)) bl <- list(bl)
    blm <- do.call(rbind, lapply(bl, function(p) as.integer(p)))
    if (ncol(blm) != 2) stop("each element of bl must be a pair c(i, j)")
    blm <- blm[!duplicated(blm), , drop = FALSE]
  }
  if (any(c(ar, ma, blm) < 1)) stop("all lags must be positive integers")
  maxlag <- max(c(ar, ma, blm, 0L))
  if (maxlag == 0 && !intercept) stop("model has no parameters")
  if (is.null(gamma)) gamma <- maxlag
  if (gamma < maxlag) stop("gamma must be at least the largest lag in the model")
  N <- length(xv)
  n <- N - gamma
  pa <- length(ar); pc <- length(ma); pb <- nrow(blm)
  q <- pa + pc + pb + as.integer(intercept)
  if (n <= q) stop("not enough observations")
  nm <- c(if (pa) paste0("a", ar), if (pc) paste0("c", ma),
          if (pb) paste0("b(", blm[, 1], ",", blm[, 2], ")"), if (intercept) "alpha")
  xs <- xv / scale

  # ---- residual recursion and its derivatives (matrix J) -------------------
  # theta is ordered (a, c, b, alpha); returns Q, the residuals and J.
  core <- function(th, xx, arl, mal, blp, jac = TRUE) {
    na <- length(arl); nc <- length(mal); nb <- nrow(blp)
    a <- th[seq_len(na)]
    cc <- th[na + seq_len(nc)]
    b <- th[na + nc + seq_len(nb)]
    al <- if (intercept) th[length(th)] else 0
    k <- length(th)
    e <- numeric(N)
    D <- if (jac) matrix(0, N, k) else NULL
    for (t in (gamma + 1):N) {
      # e_t = X_t + sum a X - alpha - sum c e - sum b X e
      et <- xx[t] - al
      if (na) et <- et + sum(a * xx[t - arl])
      if (nc) et <- et - sum(cc * e[t - mal])
      if (nb) et <- et - sum(b * xx[t - blp[, 1]] * e[t - blp[, 2]])
      if (!is.finite(et) || abs(et) > 1e100) return(list(Q = Inf))
      e[t] <- et
      if (jac) {
        # direct term d_t, then the recursive part
        d <- c(if (na) xx[t - arl], if (nc) -e[t - mal],
               if (nb) -xx[t - blp[, 1]] * e[t - blp[, 2]], if (intercept) -1)
        if (nc) for (l in seq_len(nc)) d <- d - cc[l] * D[t - mal[l], ]
        if (nb) for (r in seq_len(nb)) d <- d - b[r] * xx[t - blp[r, 1]] * D[t - blp[r, 2], ]
        D[t, ] <- d
      }
    }
    idx <- (gamma + 1):N
    list(Q = sum(e[idx]^2), e = e, J = if (jac) D[idx, , drop = FALSE] else NULL)
  }

  # ---- Gauss-Newton with step-halving --------------------------------------
  gn <- function(th, arl, mal, blp, show = FALSE) {
    cur <- core(th, xs, arl, mal, blp)
    if (!is.finite(cur$Q)) return(list(Q = Inf))
    hist <- matrix(c(th, cur$Q), nrow = 1)
    conv <- FALSE
    it <- 0
    if (show) cat(sprintf("iter %3d  Q = %.8g\n", 0, cur$Q * scale^2))
    while (it < maxit) {
      # delta = (J'J)^{-1} J'e, computed as a least-squares regression of e on J
      qrJ <- qr(cur$J)
      dl <- qr.coef(qrJ, cur$e[(gamma + 1):N])
      dl[is.na(dl)] <- 0
      lam <- 1
      repeat {
        thn <- th - lam * dl
        nw <- core(thn, xs, arl, mal, blp, jac = FALSE)
        if (is.finite(nw$Q) && nw$Q <= cur$Q) break
        lam <- lam / 2
        if (lam < 1e-10) break
      }
      if (lam < 1e-10) { conv <- TRUE; break }
      it <- it + 1
      Qold <- cur$Q
      th <- thn
      cur <- core(th, xs, arl, mal, blp)
      hist <- rbind(hist, c(th, cur$Q))
      if (show) cat(sprintf("iter %3d  Q = %.8g\n", it, cur$Q * scale^2))
      if (abs(Qold - cur$Q) / max(Qold, .Machine$double.eps) < tol) { conv <- TRUE; break }
    }
    list(theta = th, Q = cur$Q, e = cur$e, J = cur$J, conv = conv, iter = it, hist = hist)
  }

  # ---- invertibility: growth of a perturbation of the starting residuals ---
  inv_check <- function(th, mal, blp) {
    na <- length(ar); nc <- length(mal); nb <- nrow(blp)
    cc <- th[na + seq_len(nc)]
    b <- th[na + nc + seq_len(nb)]
    dl <- numeric(N)
    dl[seq_len(gamma)] <- 1
    for (t in (gamma + 1):N) {
      v <- 0
      if (nc) v <- v - sum(cc * dl[t - mal])
      if (nb) v <- v - sum(b * xs[t - blp[, 1]] * dl[t - blp[, 2]])
      if (!is.finite(v) || abs(v) > 1e100) return(Inf)
      dl[t] <- v
    }
    log10(max(abs(dl[(N - max(c(mal, blp[, 2], 1L)) + 1):N])) + 1e-300)
  }

  # keep the fit started from zero unless another start is better AND invertible
  pick <- function(fits, base) {
    best <- fits[[base]]
    for (f in fits) {
      if (!is.finite(f$Q)) next
      if (!is.finite(best$Q) || (f$Q < best$Q && f$inv < 0)) best <- f
    }
    if (!is.finite(best$Q)) NULL else best
  }

  ols_start <- function() {
    idx <- (gamma + 1):N
    Z <- cbind(if (pa) sapply(ar, function(k) -xs[idx - k]), if (intercept) 1)
    if (is.null(Z) || ncol(Z) == 0) return(numeric(0))
    as.numeric(qr.coef(qr(Z), xs[idx]))
  }

  # ---- estimation ------------------------------------------------------------
  if (!is.null(start)) {
    if (length(start) != q) stop(sprintf("start must have length %d", q))
    th0 <- as.numeric(start)
    th0[pa + pc + seq_len(pb)] <- th0[pa + pc + seq_len(pb)] * scale
    if (intercept) th0[q] <- th0[q] / scale
    fit <- gn(th0, ar, ma, blm, show = trace)
    if (!is.finite(fit$Q)) stop("estimation failed: residual recursion diverged at the starting values")
  } else {
    s0 <- ols_start()
    a0 <- s0[seq_len(pa)]
    al0 <- if (intercept) s0[length(s0)] else numeric(0)
    cur_c <- numeric(0); cur_b <- numeric(0)
    mal <- integer(0); blp <- matrix(integer(0), ncol = 2)
    fit <- NULL
    if (pc == 0 && pb == 0) fit <- gn(c(a0, al0), ar, mal, blp, show = trace)
    # add MA terms one at a time
    for (l in seq_len(pc)) {
      mal <- ma[seq_len(l)]
      grid <- unique(c(0, c_grid))
      fits <- lapply(grid, function(c0) {
        f <- gn(c(a0, cur_c, c0, al0), ar, mal, blp, show = FALSE)
        if (is.finite(f$Q)) f$inv <- inv_check(f$theta, mal, blp)
        f
      })
      best <- pick(fits, 1)
      if (is.null(best)) stop("estimation failed while adding MA terms")
      fit <- best
      a0 <- fit$theta[seq_len(pa)]
      cur_c <- fit$theta[pa + seq_len(l)]
      al0 <- if (intercept) fit$theta[length(fit$theta)] else numeric(0)
    }
    # add bilinear terms one at a time, in the order given in `bl`
    sdx <- stats::sd(xs)
    for (r in seq_len(pb)) {
      blp <- blm[seq_len(r), , drop = FALSE]
      grid <- unique(c(0, b_grid / sdx))
      fits <- lapply(grid, function(b0) {
        f <- gn(c(a0, cur_c, cur_b, b0, al0), ar, ma, blp, show = FALSE)
        if (is.finite(f$Q)) f$inv <- inv_check(f$theta, ma, blp)
        f
      })
      best <- pick(fits, 1)
      if (is.null(best)) stop("estimation failed while adding bilinear terms")
      fit <- best
      a0 <- fit$theta[seq_len(pa)]
      cur_c <- fit$theta[pa + seq_len(pc)]
      cur_b <- fit$theta[pa + pc + seq_len(r)]
      al0 <- if (intercept) fit$theta[length(fit$theta)] else numeric(0)
    }
    if (trace) fit <- gn(fit$theta, ar, ma, blm, show = TRUE)
  }

  # ---- results on the original scale ------------------------------------------
  mult <- c(rep(1, pa + pc), rep(1 / scale, pb), if (intercept) scale)
  theta <- fit$theta * mult
  names(theta) <- nm
  sigma2 <- fit$Q * scale^2 / n
  XtX <- crossprod(fit$J)
  V <- tryCatch(solve(XtX), error = function(err) matrix(NA_real_, q, q))
  V <- (fit$Q / n) * V
  V <- V * outer(mult, mult)
  dimnames(V) <- list(nm, nm)
  se <- sqrt(pmax(diag(V), 0))
  zval <- theta / se
  pval <- 2 * stats::pnorm(-abs(zval))
  res <- rep(NA_real_, N)
  res[(gamma + 1):N] <- fit$e[(gamma + 1):N] * scale
  fitted <- xv - res
  if (stats::is.ts(xo)) {
    res <- stats::ts(res, start = stats::start(xo), frequency = stats::frequency(xo))
    fitted <- stats::ts(fitted, start = stats::start(xo), frequency = stats::frequency(xo))
  }
  invlog <- inv_check(fit$theta, ma, blm)
  aic <- n * log(sigma2) + 2 * q
  hist <- fit$hist
  hist[, seq_len(q)] <- sweep(hist[, seq_len(q), drop = FALSE], 2, mult, "*")
  hist[, q + 1] <- hist[, q + 1] * scale^2
  colnames(hist) <- c(nm, "Q")
  rownames(hist) <- paste0("s=", seq_len(nrow(hist)) - 1)
  J <- fit$J %*% diag(1 / mult, q) * scale
  colnames(J) <- nm
  structure(list(coef = theta, se = se, zvalue = zval, pvalue = pval, vcov = V,
                 sigma2 = sigma2, rss = fit$Q * scale^2,
                 loglik = -n / 2 * (1 + log(2 * pi) + log(sigma2)),
                 aic = aic, naic = aic / n, n = n, N = N, gamma = gamma, npar = q,
                 residuals = res, fitted.values = fitted,
                 converged = fit$conv, iter = fit$iter, history = hist,
                 invertibility = list(log10delta = invlog, invertible = invlog < 0),
                 jacobian = J,
                 model = list(ar = ar, ma = ma, bl = blm, intercept = intercept, scale = scale),
                 x = xo, call = cl),
            class = "bilinear")
}
