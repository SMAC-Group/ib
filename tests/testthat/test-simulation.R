test_that("simulation does not work when parameters are outside their range",{
  # negative binomial regression
  library(MASS)
  fit_nb <- glm.nb(Days ~ Sex/(Age + Eth*Lrn), data = quine)
  fit_nb$theta <- -1 # theta>0!
  expect_error(simulation(fit_nb))
})
