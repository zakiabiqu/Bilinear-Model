test_that("Gauss-Newton path matches the hand-worked example", {
  x <- c(-0.6, -0.2, -1.2, -0.5, 0.0, 0.1, 1.6, 1.9, 0.8, -0.5, 0.4)
  fit <- bilinear(x, ar = 1, bl = list(c(1, 1)), intercept = FALSE, start = c(-0.5546, 0))
  h <- fit$history
  expect_equal(unname(h[1:3, "Q"]), c(6.1427, 5.0500, 4.7888), tolerance = 1e-4)
  expect_equal(unname(h[2, 1:2]), c(-0.3466, 0.3348), tolerance = 1e-3)
  expect_equal(unname(coef(fit)), c(-0.5238, 0.4136), tolerance = 1e-3)
})

test_that("recursive derivatives agree with finite differences", {
  x <- c(1.2, 0.8, 1.9, 0.4, 1.5, 1.1)
  fit <- bilinear(x, ar = 1, bl = list(c(1, 1)), start = c(-0.4, 0.3, 0.5), maxit = 0)
  expect_equal(unname(fit$jacobian[5, ]), c(1.406841, -1.919753, -0.580607), tolerance = 1e-5)
})
