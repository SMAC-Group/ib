# These functions are
# Copyright (C) 2020 S. Orso, University of Geneva
# All rights reserved.

ib.nls <- function(object, thetastart=NULL, control=list(...), extra_param = FALSE, ...){
  # controls
  control <- do.call("ibControl",control)

  # initial estimator:
  pi0 <- coef(object)

  if(extra_param) pi0 <- c(pi0, sigma(object))

  if(!is.null(thetastart)){
    if(is.numeric(thetastart) && length(thetastart) == length(pi0)){
      t0 <- thetastart
    } else {
      stop("`thetastart` must be a numeric vector of the same length as parameter
           of interest.", call.=FALSE)
    }
  } else {
    t0 <- pi0
  }


  # test diff between thetas
  p <- p0 <- length(t0)
  if(extra_param) p0 <- p - 1L
  test_theta <- control$tol + 1

  # iterator
  k <- 0L

  # create an environment for iterative bootstrap
  env_ib <- new.env(hash=F)

  # prepare data and formula for fit
  cl <- getCall(object)
  if(length(cl$formula)==1) cl$formula <- get(paste(cl$formula)) # get formula

  # retrieve data: 3 possibilities
  # 1. within object$model
  # 2. specified by an environment variable
  # 3. specified by a data.frame variable
  data <- object$model
  x_name <- attr(object$dataClasses, "names")
  y_name <- as.character(cl$formula[[2]])
  if(is.null(data)){
    data_name <- object$data
    data <- eval(data_name)
    if (!is.list(data) && !is.environment(data)) stop("'data' must be a list or an environment")
    if(is.environment(data)) data <- get(c(y_name,x_name), envir = data)
  }
  for(i in x_name) assign(i, data[[i]], env_ib)
  assign("y", data[[y_name]], env_ib)
  cl$formula[[2]] <- quote(y)
  cl$data <- quote(env_ib)
  # FIXME: add support for weights, subset, na.action, start

  # copy the object
  # FIXME: need a deep copy
  # current fix by new evaluation
  tmp_object <- eval(cl)

  if(!extra_param) std <- NULL

  # Iterative bootstrap algorithm:
  while(test_theta > control$tol && k < control$maxit){

    # update initial estimator
    tmp_object$m$setPars(t0[1:p0])
    if(extra_param) std <- t0[p]
    sim <- simulation(tmp_object,control,std)
    tmp_pi <- matrix(NA_real_,nrow=p,ncol=control$H)
    for(h in seq_len(control$H)){
      assign("y",sim[,h],env_ib)
      fit_tmp <- eval(cl)
      tmp_pi[1:p0,h] <- coef(fit_tmp)
      if(extra_param) tmp_pi[p,h] <- sigma(fit_tmp)
    }
    pi_star <- control$func(tmp_pi)

    # update value
    delta <- pi0 - pi_star
    if(extra_param) delta[p] <- exp(log(pi0[p])-log(pi_star[p]))
    t1 <- t0 + delta

    # test diff between thetas
    test_theta <- sqrt(drop(crossprod(t0-t1))/p)

    # initialize test
    if(!k) tt_old <- test_theta+1

    # Stop if no more progress
    if(tt_old <= test_theta) {break} else {tt_old <- test_theta}

    # update increment
    k <- k + 1L

    # Print info
    if(control$verbose){
      cat("Iteration:",k,"Norm between theta_k and theta_(k-1):",test_theta,"\n")
    }

    # update theta
    t0 <- t1
  }

  tmp_object$call <- object$call

  # additional metadata
  ib_warn <- NULL
  if(k>=control$maxit) ib_warn <- gettext("maximum number of iteration reached")
  if(tt_old<=test_theta) ib_warn <- gettext("objective function does not reduce")
  ib_extra <- list(
    iteration = k,
    of = sqrt(drop(crossprod(delta))),
    estimate = t0,
    test_theta = test_theta,
    ib_warn = ib_warn,
    boot = tmp_pi)

  new("IbNls",
      object = tmp_object,
      ib_extra = ib_extra)
}

#' @rdname ib
#' @details
#' For \link[stats]{nls}, if \code{extra_param=TRUE}: the variance of the residuals is
#' also corrected.
#' @example /inst/examples/eg_nls.R
#' @seealso \code{\link[stats]{nls}}
#' @export
setMethod("ib", signature = className("nls","stats"),
          definition = ib.nls)

# No native simulate method for class nls
# TODO: maybe compare with implementation in \pkg{nlraa}
simulation.nls <- function(object, control=list(...), std=NULL, ...){
  control <- do.call("ibControl",control)

  set.seed(control$seed)
  if(!exists(".Random.seed", envir = .GlobalEnv)) runif(1)

  ftd <- fitted(object)
  n <- length(ftd)
  ntot <- n * control$H
  if(is.null(std)) std <- sigma(object)
  sim <- matrix(ftd + rnorm(ntot,sd=std), ncol=control$H)
  if(control$cens) sim <- censoring(sim,control$right,control$left)
  if(control$mis) sim <- missing_at_random(sim, control$prop)
  if(control$out) sim <- outliers(sim, control$eps, control$G)
  sim
}

#' @title Simulation for nonlinear regression
#' @description simulation method for class \linkS4class{IbNls}
#' @param object an object of class \linkS4class{IbNls}
#' @param control a \code{list} of parameters for controlling the iterative procedure
#' (see \code{\link{ibControl}}).
#' @param std \code{NULL} by default; standard deviation to pass to simulation.
#' @param ... further arguments
#' @export
setMethod("simulation", signature = className("nls","stats"),
          definition = simulation.nls)
