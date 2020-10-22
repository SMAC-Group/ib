censoring <- function(y,right=NULL,left=NULL){
  if(all(!is.null(right),!is.null(left),right<left))
    stop("right-censoring must be greater than left-censoring!")
  if(!is.null(right)) y[y>right] <- right
  if(!is.null(left)) y[y<left] <- left
  y
}

# TODO: add missing, contamination

#' @title Simulation
simulation <- function(object, control=list(...), ...){
  UseMethod("simulation",object)
}

#' @importFrom stats simulate
simulation.default <- function(object, control=list(...), ...){
  control <- do.call("ibControl",control)
  sim <- simulate(object,nsim=control$H,seed=control$seed,...)
  if(control$cens) sim <- censoring(sim,control$right,control$left)
  data.matrix(sim)
}
