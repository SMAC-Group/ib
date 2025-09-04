
## generalized linear mixed-effects regression
\dontrun{
  library(lme4)
  fit_glmm <- glmer(incidence / size ~ period + (1 | herd), weights = size,
                    family = binomial, data = cbpp)
  fit_ib <- ib(fit_glmm)
  summary(fit_ib)
  ## correct for variances and correlation
  fit_ib <- ib(fit_glmm, extra_param = TRUE)
  summary(fit_ib)
}
