test_that("limited support of 'identity' link for precision parameter of betareg",{
  library(betareg)
  data("GasolineYield", package = "betareg")
  fit_beta <- betareg(yield ~ batch + temp | temp, link.phi = "identity", data = GasolineYield)
  expect_error(ib(fit_beta))
})
