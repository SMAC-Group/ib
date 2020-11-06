
## linear regression
fit_lm <- lm(speed ~ dist, data = cars)
fit_ib <- ib(fit_lm)
summary(fit_ib)
## correct for variance of residuals
fit_ib <- ib(fit_lm,var=TRUE)
summary(fit_ib)
