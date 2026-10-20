pkgload::load_all(".", quiet = TRUE)
library(ggplot2)

sun <- window(sunspot.year, 1700, 1955)
n_train <- floor(0.7 * length(sun))
train <- window(sun, end = time(sun)[n_train])
test  <- window(sun, start = time(sun)[n_train + 1])

fit <- bilinear(train, ar = c(1, 2, 9), gamma = 10,
                bl = list(c(2, 1), c(8, 1), c(1, 3), c(4, 3), c(1, 6), c(2, 4), c(3, 2)))
fc <- predict(fit, n.ahead = length(test), seed = 123)

df <- rbind(
  data.frame(year = as.numeric(time(sun)),   value = as.numeric(sun),               series = "Actual"),
  data.frame(year = as.numeric(time(train)), value = as.numeric(fit$fitted.values), series = "Fitted (training)"),
  data.frame(year = as.numeric(time(test)),  value = as.numeric(fc$pred),           series = "Forecast (testing)")
)
df$series <- factor(df$series, levels = c("Actual", "Fitted (training)", "Forecast (testing)"))

p <- ggplot(df, aes(year, value, colour = series)) +
  geom_line(linewidth = 0.6, na.rm = TRUE) +
  geom_vline(xintercept = max(time(train)) + 0.5, linetype = "dashed", colour = "grey60") +
  scale_colour_manual(values = c("Actual" = "black", "Fitted (training)" = "#2a78d6",
                                 "Forecast (testing)" = "#eb6834")) +
  labs(x = "Year", y = "Sunspot number", colour = NULL) +
  theme_minimal(base_size = 11) +
  theme(legend.position = "top")
ggsave("man/figures/README-forecast.png", p, width = 9, height = 4, dpi = 150, bg = "white")
