# These functions are
# Copyright (C) 2020 S. Orso, University of Geneva
# All rights reserved.

## Generic definition for the package

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
