censoring <- function(y,right=NULL,left=NULL){
  if(all(!is.null(right),!is.null(left),right<left))
    stop("right-censoring must be greater than left-censoring!")
  if(!is.null(right)) y[y>right] <- right
  if(!is.null(left)) y[y<left] <- left
  y
}

#' @title Control for simulation
#' @export
sim_control <- function(cens=FALSE,right=NULL,left=NULL,seed=123L,H=1L){
  if(!is.logical(cens)) stop("cens must a boolean")
  list(cens=cens,right=right,left=left,seed=seed,H=H)
}

#' @title Simulation
simulation <- function(object, control=sim_control(), ...){
  UseMethod("simulation",object)
}

simulation.default <- function(object, control=sim_control(), ...){
  sim <- simulate(object,nsim=control$H,seed=control$seed,...)
  if(control$cens) sim <- censoring(sim,control$right,control$left)
  data.matrix(sim)
}
