
[![R-CMD-check](https://github.com/SMAC-Group/ib/workflows/R-CMD-check/badge.svg)](https://github.com/SMAC-Group/ib/actions)
[![Last-changedate](https://img.shields.io/badge/last%20change-2020--12--11-green.svg)](https://github.com/SMAC-Group/ib)
[![license](https://img.shields.io/badge/license-GPL--3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0.en.html)

# Bias correction via the iterative bootstrap

This is an under-development package that proposes the iterative
bootstrap algorithm of [Kuk
(1995)](https://doi.org/10.1111/j.2517-6161.1995.tb02035.x) and further
studied by [Guerrier et al
(2019)](https://doi.org/10.1080/01621459.2017.1380031) and [Guerrier et
al (2020)](https://arxiv.org/pdf/2002.08757.pdf).

In order to install the package

``` r
## if not installed
## install.packages("remotes")
remotes::install_github("SMAC-Group/ib")
```

The `ib` package is conceived as a wrapper: an `object` that needs a
bias correction is supplied to the `ib()` function. For example, for a
negative binomial regression:

``` r
library(ib)
library(MASS)
fit_nb <- glm.nb(Days ~ Sex/(Age + Eth*Lrn), data = quine)
fit_ib1 <- ib(fit_nb)
summary(fit_ib1)

## correct for overdispersion with H=100
fit_ib2 <- ib(fit_nb, control=list(H=100), extra_param = TRUE)
summary(fit_ib2)
```

Currently we support `lm`, `glm`, `glm.nb`, `lmer`, `nls` and `vglm`
classes, as shown in the example above with the overdispersion parameter
of the negative binomial regression. More details are in `help(ib)`.

On top of `simulate`, we also consider cases where the response variable
is generated using censoring, missing at random and outliers mechanisms
(see `help(ibControl)` for more details). For example

``` r
## suppose values above 30 are censored
quine2 <- transform(quine, Days=pmin(Days,30))
fit_nb <- glm.nb(Days ~ Sex/(Age + Eth*Lrn), data = quine2)
fit_ib1 <- ib(fit_nb, control = list(cens=TRUE, right=30))
summary(fit_ib1)

## correct for overdispersion with H=100
fit_ib2 <- ib(fit_nb, control=list(H=100, cens=TRUE, right=30), extra_param = TRUE)
summary(fit_ib2)
```
