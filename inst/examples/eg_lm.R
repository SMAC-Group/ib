
## linear regression
fit_lm <- lm(disp ~ cyl + hp + wt, data = mtcars)
fit_ib <- ib(fit_lm)
summary(fit_ib)
## correct for variance of residuals
fit_ib <- ib(fit_lm, extra_param = TRUE)
summary(fit_ib)
