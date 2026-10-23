# =============================================================================
# Bilinear model for the monthly Nino 3 sea surface temperature (ERSST v5)
# Seasonal component (Fourier + trend) and subset bilinear model BL(2,2,1,1)
#
# Data   : NOAA CPC, ersst5.nino.mth.91-20.ascii (January 1950 - June 2026)
# Package: bilinear (github.com/zakiabiqu/Bilinear-Model)
# =============================================================================
library(bilinear)
library(ggplot2)
library(tseries)
library(parallel)

# ---- 1. Data and training/testing split (85:15) ------------------------------
nino <- read.table("ersst5.nino.mth.91-20.ascii", header = FALSE, skip = 1)
names(nino) <- c("tahun", "bulan", "nino12", "anom12", "nino3", "anom3",
                 "nino4", "anom4", "nino34", "anom34")
nino$tanggal <- as.Date(sprintf("%d-%02d-15", nino$tahun, nino$bulan))
all(diff((nino$tahun - 1950) * 12 + nino$bulan) == 1)    # no missing months
summary(nino$nino3)

x <- nino$nino3
N <- length(x)
n_train <- round(0.85 * N)
train <- x[1:n_train]
test <- x[(n_train + 1):N]
tgl_train <- nino$tanggal[1:n_train]
tgl_test <- nino$tanggal[(n_train + 1):N]
c(N = N, train = length(train), test = length(test))

p_data <- ggplot(nino, aes(tanggal, nino3)) +
  geom_line(linewidth = 0.3) +
  geom_vline(xintercept = tgl_test[1], linetype = "dashed", colour = "grey40") +
  labs(x = "Tahun", y = "SST Ni\u00f1o 3 (\u00b0C)") +
  theme_minimal(base_size = 10)
p_data

# ---- 2. Seasonal component: Fourier terms (period 12) and linear trend -------
t_all <- 0:(N - 1)                        # month number since January 1950
t_train <- t_all[1:n_train]
t_test <- t_all[(n_train + 1):N]
fourier_mat <- function(t, K, periode = 12) {
  Z <- do.call(cbind, lapply(1:K, function(k)
    cbind(sin(2 * pi * k * t / periode), cos(2 * pi * k * t / periode))))
  colnames(Z) <- paste0(c("sin", "cos"), rep(1:K, each = 2))
  Z
}
desain <- function(t, K, tren) {
  Z <- fourier_mat(t, K)
  if (tren) Z <- cbind(tren = t / 120, Z)   # trend in degrees per decade
  Z
}
krit <- expand.grid(K = 1:6, tren = c(FALSE, TRUE))
krit$BIC <- mapply(function(K, tren) BIC(lm(train ~ desain(t_train, K, tren))),
                   krit$K, krit$tren)
krit$dBIC <- krit$BIC - min(krit$BIC)
krit
K <- krit$K[which.min(krit$BIC)]
pakai_tren <- krit$tren[which.min(krit$BIC)]
fit_musim <- lm(train ~ desain(t_train, K, pakai_tren))
summary(fit_musim)
musim_train <- as.numeric(fitted(fit_musim))
musim_test <- as.numeric(cbind(1, desain(t_test, K, pakai_tren)) %*% coef(fit_musim))
y_train <- train - musim_train
y_test <- test - musim_test

p_musim <- ggplot(data.frame(tanggal = tgl_train, sst = train, musim = musim_train)) +
  geom_line(aes(tanggal, sst), colour = "grey55", linewidth = 0.25) +
  geom_line(aes(tanggal, musim), colour = "red", linewidth = 0.5) +
  labs(x = "Tahun", y = "SST Ni\u00f1o 3 (\u00b0C)") +
  theme_minimal(base_size = 10)
p_musim

# ---- 3. Identification on the deseasonalised training series Y ---------------
gamma <- 12                               # largest lag, fixed for every model
adf.test(y_train)
kpss.test(y_train, null = "Level")

acf_y <- acf(y_train, lag.max = 36, plot = FALSE)
pacf_y <- pacf(y_train, lag.max = 36, plot = FALSE)
batas <- qnorm(0.975) / sqrt(length(y_train))
df_acf <- rbind(data.frame(lag = acf_y$lag[-1], nilai = acf_y$acf[-1], fungsi = "ACF"),
                data.frame(lag = pacf_y$lag, nilai = pacf_y$acf, fungsi = "PACF"))
p_acf <- ggplot(df_acf, aes(lag, nilai)) +
  geom_hline(yintercept = 0, colour = "grey40") +
  geom_hline(yintercept = c(-batas, batas), linetype = "dashed", colour = "blue") +
  geom_segment(aes(xend = lag, yend = 0)) +
  facet_wrap(~fungsi, ncol = 2) +
  labs(x = "Lag (bulan)", y = "Korelasi") +
  theme_minimal(base_size = 10)
p_acf
which(abs(pacf_y$acf[1:12]) > batas)      # significant PACF lags up to 12

aic_ar <- sapply(1:12, function(p) bilinear(y_train, ar = 1:p, gamma = gamma)$aic)
round(aic_ar - min(aic_ar), 1)
p_ar_aic <- which.min(aic_ar)

# linearity check on the residuals of the AR model with the smallest AIC
fit_ar_id <- bilinear(y_train, ar = 1:p_ar_aic, gamma = gamma)
bl_diagnostics(fit_ar_id, lag = 12)
e_id <- as.numeric(fit_ar_id$residuals)
e_id[is.na(e_id)] <- 0
Box.test(e_id[-(1:gamma)]^2, lag = 12, type = "Ljung-Box")   # McLeod-Li test

# correlation between e_t and (Y_{t-i} - mean) e_{t-j}, i, j = 1, ..., 12
yc <- y_train - mean(y_train)
idx <- (gamma + 1):length(y_train)
kor <- expand.grid(i = 1:12, j = 1:12)
kor$r <- mapply(function(i, j) cor(e_id[idx], yc[idx - i] * e_id[idx - j]), kor$i, kor$j)
batas_kor <- qnorm(0.975) / sqrt(length(idx))
subset(kor, abs(r) > batas_kor)
p_heat <- ggplot(kor, aes(j, i, fill = r)) +
  geom_tile() +
  geom_point(data = subset(kor, abs(r) > batas_kor), shape = 4, size = 1.5) +
  scale_fill_gradient2(low = "blue", mid = "white", high = "red") +
  scale_x_continuous(breaks = 1:12) +
  scale_y_continuous(breaks = 1:12) +
  labs(x = "j (lag galat)", y = "i (lag Y)", fill = "korelasi") +
  coord_equal() +
  theme_minimal(base_size = 10)
p_heat

# Lagrange multiplier test (auxiliary regression, F version): residuals of a
# linear AR(p) model regressed on the AR regressors and on the products
# Y_{t-i} e_{t-1}, i = 1, ..., 5. H0: the coefficients of the products are zero.
uji_lm <- function(p) {
  f <- bilinear(y_train, ar = 1:p, gamma = gamma)
  e <- as.numeric(f$residuals)
  e[is.na(e)] <- 0
  t_ <- (gamma + 1):length(y_train)
  Ylag <- sapply(1:p, function(l) y_train[t_ - l])
  prod_ <- sapply(1:5, function(i) y_train[t_ - i] * e[t_ - 1])
  m0 <- lm(e[t_] ~ Ylag)
  m1 <- lm(e[t_] ~ Ylag + prod_)
  a <- anova(m0, m1)
  data.frame(model_nol = sprintf("AR(%d)", p), F = a$F[2], df1 = a$Df[2],
             df2 = a$Res.Df[2], p.value = a$`Pr(>F)`[2])
}
tabel_lm <- do.call(rbind, lapply(c(1, 2, 4, 5), uji_lm))   # candidate AR orders
tabel_lm

# ---- 4. Candidate models and rolling-origin evaluation ----------------------
# p: significant PACF lags (full AR 1..p), r = 0..12,
# (m, k): all b_ij with i <= m, j <= k, m <= 5, k <= 3, plus no bilinear part
mk_kand <- rbind(c(0, 0), as.matrix(expand.grid(m = 1:5, k = 1:3)))
kandidat <- merge(expand.grid(p = c(1, 2, 4, 5), r = 0:12),
                  data.frame(m = mk_kand[, 1], k = mk_kand[, 2]))
kandidat <- kandidat[order(kandidat$r, kandidat$p, kandidat$m, kandidat$k), ]
rownames(kandidat) <- NULL
nrow(kandidat)

pasangan_bl <- function(m, k) {
  if (m == 0 || k == 0) return(NULL)
  g <- expand.grid(i = 1:m, j = 1:k)
  lapply(seq_len(nrow(g)), function(s) c(g$i[s], g$j[s]))
}

# Five origins in the last 120 training months. At each origin the seasonal
# component and the model are re-estimated with the data up to the origin,
# and rolling one-step forecasts are made for the next 24 months.
h <- 24
K_ro <- 5
origin <- length(train) - (K_ro:1) * h

evaluasi <- function(s) {
  p <- kandidat$p[s]
  r <- kandidat$r[s]
  ma <- if (r > 0) 1:r else NULL
  bl <- pasangan_bl(kandidat$m[s], kandidat$k[s])
  out <- data.frame(kandidat[s, ], npar = NA, AIC = NA, sigma2 = NA, konvergen = NA,
                    iter_lanjut = FALSE, invertibel = NA, RMSE_RO = NA, waktu = NA)
  t0 <- Sys.time()
  fit <- tryCatch(bilinear(y_train, ar = 1:p, ma = ma, bl = bl, gamma = gamma,
                           b_grid = 0, c_grid = 0, maxit = 500),
                  error = function(err) NULL)
  if (!is.null(fit) && !fit$converged) {         # continue from the last estimate
    out$iter_lanjut <- TRUE
    fit <- tryCatch(bilinear(y_train, ar = 1:p, ma = ma, bl = bl, gamma = gamma,
                             start = coef(fit), maxit = 2000),
                    error = function(err) NULL)
  }
  if (!is.null(fit)) {
    out$npar <- fit$npar
    out$AIC <- fit$aic
    out$sigma2 <- fit$sigma2
    out$konvergen <- fit$converged
    out$invertibel <- fit$invertibility$invertible
    sq <- unlist(lapply(origin, function(o) {
      Z_o <- desain(t_train[1:o], K, pakai_tren)
      fm <- lm(train[1:o] ~ Z_o)
      y_o <- train[1:o] - fitted(fm)
      y_val <- train[(o + 1):(o + h)] -
        as.numeric(cbind(1, desain(t_train[(o + 1):(o + h)], K, pakai_tren)) %*% coef(fm))
      f <- tryCatch(bilinear(y_o, ar = 1:p, ma = ma, bl = bl, gamma = gamma,
                             start = coef(fit), maxit = 300),
                    error = function(err) NULL)
      if (is.null(f)) return(rep(NA, h))
      predict(f, newx = y_val)$onestep$table$error^2
    }))
    out$RMSE_RO <- sqrt(mean(sq))
  }
  out$waktu <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
  saveRDS(out, file.path(folder_kandidat, sprintf("kandidat_%03d.rds", s)))
  out
}

# every candidate is saved as soon as it is finished, so an interrupted
# search can be restarted (finished candidates are skipped)
folder_kandidat <- file.path(getwd(), "kandidat")
dir.create(folder_kandidat, showWarnings = FALSE)
file_kandidat <- file.path(folder_kandidat, sprintf("kandidat_%03d.rds", seq_len(nrow(kandidat))))
sisa <- which(!file.exists(file_kandidat))
length(sisa)
if (length(sisa) > 0) {
  cl <- makeCluster(max(1, detectCores() - 1))
  clusterEvalQ(cl, library(bilinear))
  clusterExport(cl, c("train", "y_train", "t_train", "K", "pakai_tren", "fourier_mat",
                      "desain", "gamma", "kandidat", "pasangan_bl", "origin", "h",
                      "folder_kandidat"))
  invisible(parLapplyLB(cl, sisa, evaluasi, chunk.size = 1))
  stopCluster(cl)
}
hasil <- do.call(rbind, lapply(file_kandidat, readRDS))
hasil$model <- with(hasil, sprintf("BL(%d,%d,%d,%d)", p, r, m, k))
hasil$status <- ifelse(!hasil$konvergen, "tidak konvergen",
                       ifelse(!hasil$invertibel, "tidak invertibel", "layak"))
table(hasil$status)
table(r = hasil$r, hasil$status)

# ---- 5. Model selection -------------------------------------------------------
# The ACF and PACF cannot separate the MA part from the bilinear part, so the
# model is chosen among the valid candidates by the smallest rolling-origin RMSE.
layak <- subset(hasil, status == "layak" & is.finite(RMSE_RO))
layak$peringkat <- rank(layak$RMSE_RO, ties.method = "min")
head(layak[order(layak$RMSE_RO), c("model", "npar", "AIC", "RMSE_RO")], 10)
head(layak[order(layak$AIC), c("model", "npar", "AIC", "RMSE_RO", "peringkat")], 10)
terpilih <- layak[which.min(layak$RMSE_RO), ]
terpilih

best_pr <- do.call(rbind, lapply(split(layak, list(layak$p, layak$r), drop = TRUE),
                                 function(d) d[which.min(d$RMSE_RO), ]))
p_pr <- ggplot(best_pr, aes(r, RMSE_RO, colour = factor(p))) +
  geom_hline(yintercept = min(layak$RMSE_RO), linetype = "dashed", colour = "grey50") +
  geom_line(linewidth = 0.7) +
  geom_point(size = 2) +
  scale_x_continuous(breaks = 0:12) +
  coord_cartesian(ylim = c(0.240, 0.258)) +
  labs(x = "Orde MA (r)", y = "RMSE titik asal bergulir terkecil", colour = "p") +
  theme_minimal(base_size = 10) +
  theme(legend.position = "top")
p_pr

# ---- 6. Final estimation of BL(2,2,1,1) on the training data ----------------
ar_t <- 1:terpilih$p
ma_t <- if (terpilih$r > 0) 1:terpilih$r else NULL
bl_t <- pasangan_bl(terpilih$m, terpilih$k)
duga <- function(y, start = NULL) {
  f <- bilinear(y, ar = ar_t, ma = ma_t, bl = bl_t, gamma = gamma, start = start,
                maxit = 500)
  if (!f$converged)
    f <- bilinear(y, ar = ar_t, ma = ma_t, bl = bl_t, gamma = gamma,
                  start = coef(f), maxit = 2000)
  f
}
fit_a <- duga(y_train)                                     # default multi-start
fit_b <- bilinear(y_train, ar = ar_t, ma = ma_t, bl = bl_t, gamma = gamma,
                  b_grid = 0, c_grid = 0, maxit = 500)      # as in the search
data.frame(versi = c("multi-start", "pencarian"), AIC = c(fit_a$aic, fit_b$aic),
           konvergen = c(fit_a$converged, fit_b$converged),
           invertibel = c(fit_a$invertibility$invertible, fit_b$invertibility$invertible))
ok <- function(f) f$converged && f$invertibility$invertible
fit_bl <- if (ok(fit_a) && (!ok(fit_b) || fit_a$aic <= fit_b$aic)) fit_a else fit_b
summary(fit_bl, lag = 12)
akar_ar <- polyroot(c(1, coef(fit_bl)[paste0("a", ar_t)]))
round(Mod(akar_ar), 3)                    # all > 1: the AR part is stationary

# ---- 7. Residual diagnostics ------------------------------------------------------
uji_res <- bl_diagnostics(fit_bl, lag = 12)
uji_res
e_bl <- as.numeric(na.omit(fit_bl$residuals))
mcleod_li <- Box.test(e_bl^2, lag = 12, type = "Ljung-Box")
mcleod_li
inovasi <- if (uji_res["Normality u1^2 + u2^2", "p.value"] >= 0.05) "normal" else "bootstrap"
inovasi

df_res <- data.frame(tanggal = tgl_train, residual = as.numeric(fit_bl$residuals))
p_res_t <- ggplot(df_res, aes(tanggal, residual)) +
  geom_hline(yintercept = 0, colour = "grey40") +
  geom_line(linewidth = 0.25, na.rm = TRUE) +
  labs(x = "Tahun", y = "Residual") +
  theme_minimal(base_size = 10)
p_res_h <- ggplot(data.frame(residual = e_bl), aes(residual)) +
  geom_histogram(aes(y = after_stat(density)), bins = 40, fill = "grey75", colour = "white") +
  stat_function(fun = dnorm, args = list(mean = 0, sd = sqrt(fit_bl$sigma2)), colour = "red") +
  labs(x = "Residual", y = "Densitas") +
  theme_minimal(base_size = 10)
acf_e <- acf(e_bl, lag.max = 24, plot = FALSE)
acf_e2 <- acf(e_bl^2, lag.max = 24, plot = FALSE)
batas_e <- qnorm(0.975) / sqrt(length(e_bl))
p_res_acf <- ggplot(rbind(data.frame(lag = 1:24, nilai = acf_e$acf[-1], fungsi = "ACF residual"),
                          data.frame(lag = 1:24, nilai = acf_e2$acf[-1],
                                     fungsi = "ACF kuadrat residual")), aes(lag, nilai)) +
  geom_hline(yintercept = 0, colour = "grey40") +
  geom_hline(yintercept = c(-batas_e, batas_e), linetype = "dashed", colour = "blue") +
  geom_segment(aes(xend = lag, yend = 0)) +
  facet_wrap(~fungsi, ncol = 2) +
  labs(x = "Lag (bulan)", y = "Korelasi") +
  theme_minimal(base_size = 10)

# ---- 8. Forecasting the testing period ----------------------------------------
# (a) multi-step forecast from the end of the training data, h = 1, ..., 138
fc_y <- predict(fit_bl, n.ahead = length(test), nsim = 5000, innov = inovasi,
                level = 95, seed = 123)
fc_x <- musim_test + as.numeric(fc_y$pred)
fc_lo <- musim_test + as.numeric(fc_y$lower[, "95%"])
fc_hi <- musim_test + as.numeric(fc_y$upper[, "95%"])
cakupan <- mean(test >= fc_lo & test <= fc_hi)
cakupan

# (b) rolling one-step-ahead forecasts: parameters fixed at the training
#     estimates, each forecast uses the observations up to the previous month
fit_arma <- bilinear(y_train, ar = ar_t, ma = ma_t, gamma = gamma)   # linear part
fit_ar5 <- bilinear(y_train, ar = 1:5, gamma = gamma)
galat_1 <- list(
  `BL(2,2,1,1)` = predict(fit_bl, newx = y_test)$onestep$table$error,
  `ARMA(2,2)` = predict(fit_arma, newx = y_test)$onestep$table$error,
  `AR(5)` = predict(fit_ar5, newx = y_test)$onestep$table$error,
  Persistensi = y_test - c(y_train[n_train], y_test[-length(y_test)]),
  Klimatologi = y_test)
akurasi <- function(g) c(RMSE = sqrt(mean(g^2)), MAE = mean(abs(g)))
tabel_1 <- t(sapply(galat_1, akurasi))
round(tabel_1, 4)

# accuracy in El Nino months (observed anomaly Y > 1 degree) and in other months
el <- y_test > 1
sum(el)
tabel_el <- t(sapply(galat_1[c("BL(2,2,1,1)", "ARMA(2,2)", "AR(5)")], function(g)
  c(El_Nino = sqrt(mean(g[el]^2)), lainnya = sqrt(mean(g[!el]^2)))))
round(tabel_el, 4)

ramal_1 <- test - galat_1[["BL(2,2,1,1)"]]
fitted_x <- musim_train + as.numeric(fit_bl$fitted.values)
df_plot <- rbind(
  data.frame(tanggal = nino$tanggal, nilai = x, seri = "Aktual"),
  data.frame(tanggal = tgl_train, nilai = fitted_x, seri = "Nilai dugaan (data latih)"),
  data.frame(tanggal = tgl_test, nilai = fc_x, seri = "Ramalan (data uji)"))
df_plot$seri <- factor(df_plot$seri, levels = c("Aktual", "Nilai dugaan (data latih)",
                                                "Ramalan (data uji)"))
p_fc <- ggplot() +
  geom_ribbon(data = data.frame(tanggal = tgl_test, lo = fc_lo, hi = fc_hi),
              aes(tanggal, ymin = lo, ymax = hi), fill = "blue", alpha = 0.12) +
  geom_line(data = df_plot, aes(tanggal, nilai, colour = seri), linewidth = 0.3, na.rm = TRUE) +
  geom_vline(xintercept = tgl_test[1], linetype = "dashed", colour = "grey40") +
  scale_colour_manual(values = c("Aktual" = "black", "Nilai dugaan (data latih)" = "red",
                                 "Ramalan (data uji)" = "blue")) +
  coord_cartesian(xlim = as.Date(c("2005-01-01", "2026-07-01"))) +
  labs(x = "Tahun", y = "SST Ni\u00f1o 3 (\u00b0C)", colour = NULL) +
  theme_minimal(base_size = 10) +
  theme(legend.position = "top")
p_fc

p_1 <- ggplot(rbind(data.frame(tanggal = tgl_test, nilai = test, seri = "Aktual"),
                    data.frame(tanggal = tgl_test, nilai = ramal_1,
                               seri = "Ramalan satu langkah")),
              aes(tanggal, nilai, colour = seri)) +
  geom_line(linewidth = 0.4) +
  scale_colour_manual(values = c("Aktual" = "black", "Ramalan satu langkah" = "blue")) +
  labs(x = "Tahun", y = "SST Ni\u00f1o 3 (\u00b0C)", colour = NULL) +
  theme_minimal(base_size = 10) +
  theme(legend.position = "top")
p_1

# ---- 9. Forecast for July - December 2026 using all data ----------------------
# same structure, parameters re-estimated on all 918 months
fit_musim_all <- lm(x ~ desain(t_all, K, pakai_tren))
musim_semua <- as.numeric(fitted(fit_musim_all))
y_semua <- x - musim_semua
fit_semua <- duga(y_semua, start = coef(fit_bl))
c(konvergen = fit_semua$converged, invertibel = fit_semua$invertibility$invertible)
uji_semua <- bl_diagnostics(fit_semua, lag = 12)
uji_semua
inovasi_semua <- if (uji_semua["Normality u1^2 + u2^2", "p.value"] >= 0.05) "normal" else
  "bootstrap"

h_depan <- 6
tgl_depan <- seq(as.Date("2026-07-15"), by = "month", length.out = h_depan)
musim_depan <- as.numeric(cbind(1, desain(N:(N + h_depan - 1), K, pakai_tren)) %*%
                            coef(fit_musim_all))
fc_depan <- predict(fit_semua, n.ahead = h_depan, nsim = 20000, innov = inovasi_semua,
                    level = c(80, 95), seed = 2026)

# naive chaining (forecasts plugged in as data, future errors set to zero)
cf <- coef(fit_semua)
y_ext <- c(y_semua, rep(NA, h_depan))
e_ext <- c(as.numeric(fit_semua$residuals), rep(0, h_depan))
e_ext[is.na(e_ext)] <- 0
for (s in N + 1:h_depan) {
  v <- cf["alpha"] - sum(cf[paste0("a", ar_t)] * y_ext[s - ar_t])
  if (length(ma_t)) v <- v + sum(cf[paste0("c", ma_t)] * e_ext[s - ma_t])
  for (u in seq_along(bl_t))
    v <- v + cf[grep("^b\\(", names(cf))][u] * y_ext[s - bl_t[[u]][1]] * e_ext[s - bl_t[[u]][2]]
  y_ext[s] <- v
}
fit_arma_semua <- bilinear(y_semua, ar = ar_t, ma = ma_t, gamma = gamma)
ar_depan <- as.numeric(predict(fit_arma_semua, n.ahead = h_depan, nsim = 20000,
                               innov = inovasi_semua, seed = 2026)$pred)

ramalan_2026 <- data.frame(
  bulan = format(tgl_depan, "%b %Y"),
  musiman = musim_depan,
  anomali = as.numeric(fc_depan$pred),
  SST = musim_depan + as.numeric(fc_depan$pred),
  bawah_95 = musim_depan + as.numeric(fc_depan$lower[, "95%"]),
  atas_95 = musim_depan + as.numeric(fc_depan$upper[, "95%"]),
  naif = musim_depan + y_ext[N + 1:h_depan],
  ARMA = musim_depan + ar_depan)
ramalan_2026[, -1] <- round(ramalan_2026[, -1], 2)
ramalan_2026

sejak <- nino$tanggal >= as.Date("2022-01-01")
df_fc6 <- data.frame(tanggal = c(nino$tanggal[N], tgl_depan),
                     nilai = c(x[N], ramalan_2026$SST),
                     lo95 = c(x[N], musim_depan + as.numeric(fc_depan$lower[, "95%"])),
                     hi95 = c(x[N], musim_depan + as.numeric(fc_depan$upper[, "95%"])),
                     lo80 = c(x[N], musim_depan + as.numeric(fc_depan$lower[, "80%"])),
                     hi80 = c(x[N], musim_depan + as.numeric(fc_depan$upper[, "80%"])))
p_2026 <- ggplot() +
  geom_ribbon(data = df_fc6, aes(tanggal, ymin = lo95, ymax = hi95), fill = "blue", alpha = 0.10) +
  geom_ribbon(data = df_fc6, aes(tanggal, ymin = lo80, ymax = hi80), fill = "blue", alpha = 0.15) +
  geom_line(data = data.frame(tanggal = c(nino$tanggal[sejak], tgl_depan),
                              nilai = c(musim_semua[sejak], musim_depan)),
            aes(tanggal, nilai), colour = "grey55", linetype = "dashed", linewidth = 0.4) +
  geom_line(data = data.frame(tanggal = nino$tanggal[sejak], nilai = x[sejak]),
            aes(tanggal, nilai), colour = "black", linewidth = 0.6) +
  geom_line(data = df_fc6, aes(tanggal, nilai), colour = "blue", linewidth = 0.8) +
  geom_point(data = df_fc6[-1, ], aes(tanggal, nilai), colour = "blue", size = 1.8) +
  labs(x = "Tahun", y = "SST Ni\u00f1o 3 (\u00b0C)") +
  theme_minimal(base_size = 10)
p_2026
