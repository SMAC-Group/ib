# These functions are
# Copyright (C) 2020 S. Orso, University of Geneva
# All rights reserved.

## Generic definition for the package

## Define method for accessing the object for class union "Ib"
#' @title Accessor to the object in class union "Ib"
#' @description
#' Method for obtaining a fitted model within
#' any object of class union \linkS4class{Ib}.
#' @param x an object of class union "Ib"
#' @seealso \linkS4class{Ib}
#' @export
setGeneric("getObject",
           function(x) standardGeneric("getObject"),
           signature = "x",
           package = "ib")

#' @rdname getObject
#' @export
setMethod("getObject",
          "Ib",
          definition = function(x) x@object)

## Define method for accessing an extra part for class union "Ib"
#' @title Accessor to an extra part in class union "Ib"
#' @param x an object of class union "Ib"
#' @description
#' Method for obtaining a extra values generated by
#' the iterative bootstrap procedure within
#' any object of class union \linkS4class{Ib}.
#' @seealso \linkS4class{Ib}
#' @export
setGeneric("getExtra",
           function(x) standardGeneric("getExtra"),
           signature = "x",
           package = "ib")

#' @rdname getExtra
#' @export
setMethod("getExtra",
          "Ib",
          definition = function(x) x@ib_extra)

## Define method for accessing the estimate for class union "Ib"
#' @title Accessor to the object in class union "Ib"
#' @description
#' Method for obtaining estimates from fitted model within
#' any object of class union \linkS4class{Ib}.
#' @param x an object of class union "Ib"
#' @seealso \linkS4class{Ib}
#' @export
setGeneric("getEst",
           function(x) standardGeneric("getEst"),
           signature = "x",
           package = "ib")

#' @rdname getEst
#' @export
setMethod("getEst",
          "Ib",
          definition = function(x) x@ib_extra$estimate)

## Define show method for class "Ib"
show.ib <- function(object){
  x <- getObject(object)
  print(x)
}

#' @title Method for printing object in class union "Ib"
#' @param object an object of class union "Ib"
#' @seealso \linkS4class{Ib}
#' @export
  #' @importFrom methods show
setMethod("show",
          "Ib",
          definition = show.ib)

## Generic for simulating from the object (internal use)
setGeneric("simulation",
           function(object, control=list(...), ...) standardGeneric("simulation"),
           signature = "object",
           package = "ib")

#' @importFrom stats simulate
simulation.default <- function(object, control=list(...), ...){
  control <- do.call("ibControl",control)
  sim <- simulate(object,nsim=control$H,seed=control$seed,...)
  if(control$cens) sim <- censoring(sim,control$right,control$left)
  if(control$mis) sim <- missing_at_random(sim, control$prop)
  if(control$out) sim <- outliers(sim, control$eps, control$G)
  data.matrix(sim)
}

## Define plot method for class "Ib"
plot.ib <- function(x, y = NULL, ...){
  plot(getObject(x), ...)
}

#' @title Method for plotting an object in class union "Ib"
#' @param x an object of class union "Ib"
#' @param y not used
#' @param ... further arguments to pass to \code{plot}
#' @seealso \linkS4class{Ib}, \link[stats]{plot.lm}
#' @export
setMethod("plot",
          "Ib",
          definition = plot.ib)

## Define residuals method for class "Ib"
residuals.ib <- function(object, ...){
  residuals(getObject(object), ...)
}

#' @title Method for extracting residuals from an object in class union "Ib"
#' @param object an object of class union "Ib"
#' @param ... further arguments to pass to \code{residuals}
#' @seealso \linkS4class{Ib}, \link[stats]{residuals}
#' @export
setMethod("residuals",
          "Ib",
          definition = residuals.ib)

## Define predict method for class "Ib"
predict.ib <- function(object, ...){
  predict(getObject(object), ...)
}

#' @title Method for making predictions from an object in class union "Ib"
#' @param object an object of class union "Ib"
#' @param ... further arguments to pass to \code{predict}
#' @seealso \linkS4class{Ib}, \link[stats]{predict}
#' @export
setMethod("predict",
          "Ib",
          definition = predict.ib)

## Define coef method for class "Ib"
coef.ib <- function(object, ...){
  coef(getObject(object), ...)
}

#' @title Method for extracting coefficients from an object in class union "Ib"
#' @param object an object of class union "Ib"
#' @param ... further arguments to pass to \code{coef}
#' @seealso \linkS4class{Ib}, \link[stats]{coef}
#' @export
setMethod("coef",
          "Ib",
          definition = coef.ib)

## Define fitted method for class "Ib"
fitted.ib <- function(object, ...){
  fitted(getObject(object), ...)
}

#' @title Method for extracting fitted values from an object in class union "Ib"
#' @param object an object of class union "Ib"
#' @param ... further arguments to pass to \code{fitted}
#' @seealso \linkS4class{Ib}, \link[stats]{fitted.values}
#' @export
setMethod("fitted",
          "Ib",
          definition = fitted.ib)

## Define effects method for class "Ib"
effects.ib <- function(object, ...){
  effects(getObject(object), ...)
}

#' @title Method for extracting effects from an object in class union "Ib"
#' @param object an object of class union "Ib"
#' @param ... further arguments to pass to \code{effects}
#' @seealso \linkS4class{Ib}, \link[stats]{effects}
#' @export
setMethod("effects",
          "Ib",
          definition = effects.ib)

## Define vcov method for class "Ib"
vcov.ib <- function(object, ...){
  vcov(getObject(object), ...)
}

#' @title Method for calculating covariance matrix from an object in class union "Ib"
#' @param object an object of class union "Ib"
#' @param ... further arguments to pass to \code{vcov}
#' @seealso \linkS4class{Ib}, \link[stats]{vcov}
#' @export
setMethod("vcov",
          "Ib",
          definition = vcov.ib)
