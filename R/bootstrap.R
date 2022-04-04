# These functions are
# Copyright (C) 2020 S. Orso, University of Geneva
# All rights reserved.

## Parametric bootstrap

#' @title
#' Parametric bootstrap
#' @description
#' Method for generating parametric bootstrap estimates from
#' a fitted model.
#' @param object an \code{object} representing a fitted model (see 'Details').
#' @param B an \code{integer} for number of bootstrap replicates (default
#' 1,000).
#' @param extra_param if \code{TRUE}, bootstrap is also performed for extra parameters
#' (see 'Details').
#' @param ... additional optional arguments to pass to \code{ibControl}.
#' @return
#' A \code{matrix} p (size of parameter) times B of bootstrapped estimates.
#' @details
#' This method is a simple wrapper around the \code{ib} method
#' where number of iterations is set to 1.
#' @seealso \code{\link{ib}}, \code{\link{ibControl}}
#' @author Samuel Orso
#' @export
bootstrap <- function(object, B = 1e3, extra_param = FALSE, ...){
  fit <- ib(object, control = list(H = B, maxit = 1L, ...), extra_param = extra_param)
  extra <- getExtra(fit)
  extra$boot
}
