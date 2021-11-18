test_that("'identity' not supported for precision parameter of betareg",{
  library(betareg)
  data("GasolineYield", package = "betareg")
  fit_beta <- betareg(yield ~ batch + temp, data = GasolineYield)
  expect_error(ib(fit_beta))
})
