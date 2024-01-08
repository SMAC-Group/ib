
## linear mixed-effects regression
\dontrun{
library(lme4)
fit_lmm <- lmer(Reaction ~ Days + (Days | Subject), data = sleepstudy, REML = FALSE)
fit_ib <- ib(fit_lmm)
summary(fit_ib)
## correct for variances and correlation
fit_ib <- ib(fit_lmm, extra_param = TRUE)
summary(fit_ib)
}
