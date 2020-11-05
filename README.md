
[![Travis-CI Build
Status](https://travis-ci.com/SMAC-Group/ib.svg?branch=master)](https://travis-ci.com/github/SMAC-Group/ib)
[![Last-changedate](https://img.shields.io/badge/last%20change-2020--11--05-green.svg)](https://github.com/SMAC-Group/ib)
[![license](https://img.shields.io/badge/license-GPL--3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0.en.html)

# Bias correction via the iterative bootstrap

This is an under-development package that proposes the iterative
bootstrap algorithm of [Kuk
(1995)](https://doi.org/10.1111/j.2517-6161.1995.tb02035.x) and further
studied by [Guerrier et al
(2019)](https://doi.org/10.1080/01621459.2017.1380031) and [Guerrier et
al (2020)](https://arxiv.org/pdf/2002.08757.pdf).

It is conceived as a wrapper: an `object` that needs a bias correction
is supplied to the `ib()` function. For example:

``` r
library(ib)

fit_lm <- fit(mtcars$mpg ~ mtcars$hp)
ib(fit_lm)   
```

Currently we support only `object` which possess respond to `getCall()`,
`model.matrix()` and `simulate()`. For example we support `lm()`,
`glm()` (with the exception of the `quasi-` families), `MASS::glm.nb()`.
More to come.
