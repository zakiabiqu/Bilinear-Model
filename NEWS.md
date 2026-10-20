# bilinear 0.1.0

* Initial release.
* `bilinear()` fits subset bilinear models by conditional least squares
  (Gauss–Newton with recursive derivatives and step-halving), with staged
  multi-start and an invertibility check.
* Methods: `print()`, `summary()`, `coef()`, `vcov()`, `predict()`
  (exact one-step forecast, Monte Carlo for longer horizons, rolling
  one-step forecasts on new data).
* `bl_diagnostics()`: W test on squared residuals, Lomnicki skewness and
  kurtosis tests, joint normality test, Ljung-Box test.
* Test suite reproducing results of Subba Rao & Gabr (1984).
