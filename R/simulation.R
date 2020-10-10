censoring <- function(y,right=NULL,left=NULL){
  if(all(!is.null(right),!is.null(left),right<left))
    stop("right-censoring must be greater than left-censoring!")
  if(!is.null(right)) y[y>right] <- right
  if(!is.null(left)) y[y<left] <- left
  y
}

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

# adapted from stats::simulate.lm
simulation.lm <- function(object, control=list(...), std=NULL, ...){
  control <- do.call("ibControl",control)
  if (!exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE))
    runif(1)
  if (is.null(control$seed))
    RNGstate <- get(".Random.seed", envir = .GlobalEnv)
  else {
    R.seed <- get(".Random.seed", envir = .GlobalEnv)
    set.seed(control$seed)
    RNGstate <- structure(control$seed, kind = as.list(RNGkind()))
    on.exit(assign(".Random.seed", R.seed, envir = .GlobalEnv))
  }
  ftd <- fitted(object)
  n <- length(ftd)
  ntot <- n * control$H
  if(is.null(std)) std <- sigma(object)
  sim <- matrix(ftd + rnorm(ntot,sd=std), ncol=control$H)
  if(control$cens) sim <- censoring(sim,control$right,control$left)
  sim
}

# adapted from stats::simulate.lm
simulation.glm <- function(object, control=list(...), shape=NULL, ...){
  fam <- object$family$family
  if (is.null(object$family$simulate)) stop(gettextf("family '%s' not implemented",fam), domain = NA)

  control <- do.call("ibControl",control)
  if (!exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE))
    runif(1)
  if (is.null(control$seed))
    RNGstate <- get(".Random.seed", envir = .GlobalEnv)
  else {
    R.seed <- get(".Random.seed", envir = .GlobalEnv)
    set.seed(control$seed)
    RNGstate <- structure(control$seed, kind = as.list(RNGkind()))
    on.exit(assign(".Random.seed", R.seed, envir = .GlobalEnv))
  }
  sim <- switch(fam,
                Gamma = {
                  if(is.null(shape)){
                    matrix(object$family$simulate(object,control$H), ncol=control$H)
                    } else {
                    matrix(simulate_gamma(object,control$H,shape), ncol=control$H)}},
                matrix(object$family$simulate(object,control$H), ncol=control$H)
         )
  if(control$cens) sim <- censoring(sim,control$right,control$left)
  sim
}

simulate_gamma <- function (object, nsim, shape){
  wts <- object$prior.weights
  if (any(wts != 1))
    message("using weights as shape parameters")
  ftd <- fitted(object)
  shape <- shape * wts
  rgamma(nsim * length(ftd), shape = shape, rate = shape/ftd)
}
