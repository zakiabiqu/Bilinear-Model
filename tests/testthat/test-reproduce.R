# Reproduce the models reported by Subba Rao and Gabr (1984, Chapter 6).

test_that("sunspot model reproduces sigma^2 = 124.33 and NAIC = 4.927", {
  sun <- window(sunspot.year, 1700, 1920)
  fit <- bilinear(sun, ar = c(1, 2, 9), gamma = 10,
                  bl = list(c(2, 1), c(8, 1), c(1, 3), c(4, 3), c(1, 6), c(2, 4), c(3, 2)))
  expect_equal(fit$n, 211)
  expect_equal(fit$sigma2, 124.33, tolerance = 1e-4)
  expect_equal(fit$naic, 4.927, tolerance = 1e-3)
  expect_equal(unname(coef(fit)[c("a1", "a2", "a9", "alpha")]),
               c(-1.5012, 0.7670, -0.1152, 6.886), tolerance = 1e-3)
  expect_equal(unname(coef(fit)["b(2,1)"]), -0.01458, tolerance = 1e-3)
  expect_true(fit$converged)
  expect_true(fit$invertibility$invertible)
})

test_that("one-step forecasts for 1921-1955 reproduce MSE = 123.77", {
  sun <- window(sunspot.year, 1700, 1920)
  fit <- bilinear(sun, ar = c(1, 2, 9), gamma = 10,
                  bl = list(c(2, 1), c(8, 1), c(1, 3), c(4, 3), c(1, 6), c(2, 4), c(3, 2)))
  fc <- predict(fit, newx = window(sunspot.year, 1921, 1955), nsim = 10)
  expect_equal(fc$onestep$mse, 123.77, tolerance = 1e-4)
  expect_equal(fc$onestep$table$forecast[1:3], c(24.55, 11.43, 10.08), tolerance = 1e-3)
})

test_that("log lynx model reproduces sigma^2 = 0.0223", {
  ly <- log10(window(lynx, 1821, 1920))
  fit <- bilinear(ly, ar = c(1, 2, 3, 4, 9, 12), gamma = 12,
                  bl = list(c(3, 9), c(9, 9), c(6, 2), c(1, 1), c(2, 7), c(4, 2)))
  expect_equal(fit$sigma2, 0.0223, tolerance = 2e-3)
  expect_equal(fit$naic, -3.508, tolerance = 1e-3)
})

test_that("scaling the data does not change the estimates", {
  ly <- window(lynx, 1821, 1920)
  f1 <- bilinear(ly, ar = c(1, 2), bl = list(c(1, 1)), scale = 1)
  f2 <- bilinear(ly, ar = c(1, 2), bl = list(c(1, 1)), scale = 100)
  expect_equal(coef(f1), coef(f2), tolerance = 1e-5)
  expect_equal(f1$sigma2, f2$sigma2, tolerance = 1e-6)
})
