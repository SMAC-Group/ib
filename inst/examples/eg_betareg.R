
## beta regression
library(betareg)
data("GasolineYield", package = "betareg")
fit_beta <- betareg(yield ~ batch + temp, data = GasolineYield)
fit_ib <- ib(fit_beta)

# precision parameter can also depend on covariates
fit_beta <- betareg(yield ~ batch + temp | temp, data = GasolineYield)
fit_ib <- ib(fit_beta)
