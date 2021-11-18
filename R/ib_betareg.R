# These functions are
# Copyright (C) 2021 S. Orso, University of Geneva
# All rights reserved.

#' @importFrom betareg betareg
ib.betareg <- function(object, thetastart=NULL, control=list(...), ...){
  # controls
  control <- do.call("ibControl",control)

  if(!object$phi) stop("Only implemented with precision parameters")

  # initial estimator:
  # regression coefficients + phi (if phi is set to TRUE (default))
  pi0 <- coef(object)
  p <- length(pi0)

  # starting value
  if(!is.null(thetastart)){
    if(is.numeric(thetastart) && length(thetastart) == length(pi0)){
      t0 <- thetastart
    } else {
      stop("`thetastart` must be a numeric vector of the same length as
           parameter of interest.", call.=FALSE)
    }
  } else {
    t0 <- pi0
  }

  # test diff between thetas
  test_theta <- control$tol + 1

  # iterator
  k <- 0L

  # create an environment for iterative bootstrap
  env_ib <- new.env(hash=F)

  # prepare data and formula for fit
  cl <- getCall(object)
  if(length(cl$formula)==1) cl$formula <- get(paste(cl$formula)) # get formula
  intercept_only <- cl$formula[[3]] == 1 # check for intercept only models
  mf <- model.frame(object)
  mt <- terms(object)
  if(!intercept_only){
    x <- if(!is.empty.model(mt)) model.matrix(mt, mf, object$contrasts)
    # check if model has an intercept
    has_intercept <- attr(mt,"intercept")
    if(has_intercept){
      # remove intercept from design
      x <- x[,!grepl("Intercept",colnames(x))]
      cl$formula <- quote(y~x)
    } else {
      cl$formula <- quote(y~x-1)
    }
    assign("x",x,env_ib)
  } else{
    cl$formula <- quote(y~1)
  }
  o <- as.vector(model.offset(mf))
  if(!is.null(o)) assign("o",o,env_ib)
  cl$data <- NULL
  # add an offset
  if(!is.null(o)) cl$offset <- quote(o)
  # FIXME: add support for weights, subset, na.action, start,
  #        etastart, mustart, contrasts

  # copy the object
  tmp_object <- object

  # Iterative bootstrap algorithm:
  while(test_theta > control$tol && k < control$maxit){
    # update initial estimator
    tmp_object$coefficients <- t0
    sim <- simulation(tmp_object,control)
    tmp_pi <- matrix(NA_real_,nrow=p,ncol=control$H)
    for(h in seq_len(control$H)){
      assign("y",sim[,h],env_ib)
      fit_tmp <- tryCatch(error = function(cnd) NULL, {eval(cl,env_ib)})
      if(is.null(fit_tmp)) next
      tmp_pi[1:p0,h] <- coef(fit_tmp)
    }
    pi_star <- control$func(tmp_pi)

    # update value
    delta <- pi0 - pi_star
    t1 <- t0 + delta
    if(control$constraint) t1[p] <- exp(log(t0[p]) + log(pi0[p]) - log(pi_star[p]))

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

  # update betareg object
  eta <- predict.glm(tmp_object)
  mu <- object$family$linkinv(eta)
  dev <- sum(object$family$dev.resids(object$y,mu,object$prior.weights))

  tmp_object$linear.predictors <- eta
  tmp_object$fitted.values <- mu
  tmp_object$residuals <- (object$y - mu)/object$family$mu.eta(eta)
  tmp_object$call <- object$call
  tmp_object$deviance <- dev
  tmp_object$aic <- object$family$aic(object$y, length(object$prior.weights)-sum(object$prior.weights == 0),
                                      mu, object$prior.weights, dev) + 2 * object$rank

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

  new("IbBetareg",
      object = tmp_object,
      ib_extra = ib_extra)
}

#' @rdname ib
#' @details
#' For \link[betareg]{betareg}, \code{extra_param} is currently not available
#' as by default the precision parameter 'phi' is corrected.
#' @seealso \code{\link[betareg]{betareg}}
#' @example /inst/examples/eg_betareg.R
#' @export
setMethod("ib", className("betareg", "betareg"),
          definition = ib.betareg)

simulation.betareg <- function(object, control=list(...), extra=NULL, ...){
  control <- do.call("ibControl",control)

  set.seed(control$seed)
  if(!exists(".Random.seed", envir = .GlobalEnv)) runif(1)

  # user-defined simulation method
  if(!is.null(control$sim)){
    sim <- control$sim(object, control, extra, ...)
    return(sim)
  }



  if(control$cens) sim <- censoring(sim,control$right,control$left)
  if(control$mis) sim <- missing_at_random(sim, control$prop)
  if(control$out) sim <- outliers(sim, control$eps, control$G)
  sim
}

#' @title Simulation for a beta regression
#' @description simulation method for class \linkS4class{IbBetareg}
#' @param object an object of class \linkS4class{IbBetareg}
#' @param control a \code{list} of parameters for controlling the iterative procedure
#' (see \code{\link{ibControl}}).
#' @param extra \code{NULL} by default; extra parameters to pass to simulation.
#' @param ... further arguments
#' @export
setMethod("simulation", signature = className("betareg","betareg"),
          definition = simulation.betareg)
