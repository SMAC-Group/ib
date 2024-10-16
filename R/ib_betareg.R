# These functions are
# Copyright (C) 2021 S. Orso, University of Geneva
# All rights reserved.

#' @importFrom betareg betareg
#' @importFrom Formula Formula as.Formula
#' @importFrom stats delete.response
ib.betareg <- function(object, thetastart=NULL, control=list(...), ...){
  # Currently support is limited to model with precision parameters ...
  if(!object$phi) stop("Only implemented with precision parameters")

  # ... and without 'link.phi="identity"' (to avoid imposing positivity constraint)
  # if(object$link$precision$name == "identity") stop("'link.phi'='identity' not supported")

  # controls
  control <- do.call("ibControl",control)

  # initial estimator:
  # regression coefficients for mean and precision parameters
  pi0 <- coef(object)
  p <- length(pi0)
  p_mean <- length(object$coefficients$mean)
  p_prec <- length(object$coefficients$precision)
  if(p_mean + p_prec != p) stop("unexpected number of parameters")
  id_mean <- vector("logical",p)
  id_mean[seq_len(p_mean)] <- TRUE
  id_prec <- !id_mean
  phiIdentity <- p_prec == 1 && object$link$precision$name == "identity"
  if(object$link$precision$name == "identity" && p_prec > 1)
    stop("limited support of 'link.phi'='identity'")

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
  mf <- model.frame(object)
  names(mf)[1] <- "y"
  assign("data",mf,env_ib)
  cl <- getCall(object)
  cl$data <- quote(data)
  if(length(cl$formula)==1) cl$formula <- get(paste(cl$formula)) # get formula
  cl$formula[[2]] <- quote(y)
  # FIXME: add support for weights, subset, na.action, offset,
  #        contrasts
  # Extract x and z
  formula <- as.Formula(cl$formula)
  if(length(formula)[2L] < 2L) {
    formula <- as.Formula(formula(formula), ~ 1)
  } else {
    if(length(formula)[2L] > 2L) {
      formula <- Formula(formula(formula, rhs = 1:2))
    }
  }
  mtX <- terms(formula, rhs = 1L)
  mtZ <- delete.response(terms(formula, rhs = 2L))
  x <- model.matrix(mtX, mf)
  z <- model.matrix(mtZ, mf)

  # copy the object
  tmp_object <- object

  # copy the control
  control1 <- control
  control1$H <- 1L
  linkinv <- object$link$mean$linkinv

  # initial value
  diff <- rep(NA_real_, control$maxit)

  # Iterative bootstrap algorithm:
  while(test_theta > control$tol && k < control$maxit){
    # update object for simulation
    if(k!=0){
      eta <- as.vector(x %*% t0[id_mean])
      mu <- linkinv(eta)
      tmp_object$fitted.values <- mu
      tmp_object$coefficients$mean <- t0[id_mean]
      tmp_object$coefficients$precision <- t0[id_prec]
    }

    # approximate
    tmp_pi <- matrix(NA_real_,nrow=p,ncol=control$H)
    for(h in seq_len(control$H)){
      control1$seed <- control$seed + h
      sim <- simulation(tmp_object,control1)
      env_ib$data$y <- sim
      fit_tmp <- tryCatch(error = function(cnd) NULL, {eval(cl,env_ib)})
      iter <- 1L
      while(is.null(fit_tmp) && iter < 10L){
        control1$seed <- control$seed + control$H * h + iter
        sim <- simulation(tmp_object,control1)
        env_ib$data$y <- sim
        fit_tmp <- tryCatch(error = function(cnd) NULL, {eval(cl,env_ib)})
        iter <- iter + 1L
      }
      if(is.null(fit_tmp)) next
      tmp_pi[,h] <- coef(fit_tmp)
    }
    pi_star <- control$func(tmp_pi)

    # update value
    delta <- pi0 - pi_star
    t1 <- t0 + delta
    if(phiIdentity && control$constraint) t1[p] <- exp(log(t0[p]) + log(pi0[p]) - log(pi_star[p]))

    # test diff between thetas
    test_theta <- sum(delta^2)
    if(k>0) diff[k] <- test_theta

    # initialize test
    if(!k) tt_old <- test_theta+1

    # Alternative stopping criteria, early stop :
    if(control$early_stop){
      if(tt_old <= test_theta){
        warning("Algorithm stopped because the objective function does not reduce")
        break
      }
    }

    # Alternative stopping criteria, "statistically flat progress curve" :
    if(k > 10L){
      try1 <- diff[k:(k-10)]
      try2 <- k:(k-10)
      if(var(try1)<=1e-3) break
      mod <- lm(try1 ~ try2)
      if(summary(mod)$coefficients[2,4] > 0.2) break
    }

    # update increment
    k <- k + 1L

    # Print info
    if(control$verbose){
      cat("Iteration:",k,"Norm between theta_k and theta_(k-1):",test_theta,"\n")
    }

    # update theta
    t0 <- t1
  }
  # warning for reaching max number of iterations
  if(k>=control$maxit) warning("maximum number of iteration reached")

  # update betareg object
  # FIXME: update object$loglik, object$scoring, object$residuals,
  # object$fitted.values, object$pseudo.r.squared

  # eta <- predict(tmp_object)
  # mu <- object$link$mean$linkinv(eta)
  # dev <- sum(object$family$dev.resids(object$y,mu,object$prior.weights))
  #
  # tmp_object$linear.predictors <- eta
  # tmp_object$fitted.values <- mu
  # tmp_object$residuals <- (object$y - mu)/object$family$mu.eta(eta)
  # tmp_object$call <- object$call
  # tmp_object$deviance <- dev
  # tmp_object$aic <- object$family$aic(object$y, length(object$prior.weights)-sum(object$prior.weights == 0),
  #                                     mu, object$prior.weights, dev) + 2 * object$rank

  # additional metadata
 ib_extra <- list(
    iteration = k,
    of = sqrt(drop(crossprod(delta))),
    estimate = t0,
    test_theta = test_theta,
    boot = tmp_pi)

  new("IbBetareg",
      object = tmp_object,
      ib_extra = ib_extra)
}

#' @rdname ib
#' @details
#' For \link[betareg]{betareg}, \code{extra_param} is not available
#' as by default mean and precision parameters are corrected.
#' Currently the 'identity' link function is not supported for precision
#' parameters.
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

  sim <- matrix(simulate_betareg(object,control$H), ncol=control$H)

  if(control$cens) sim <- censoring(sim,control$right,control$left)
  if(control$mis) sim <- missing_at_random(sim, control$prop)
  if(control$out) sim <- outliers(sim, control$eps, control$G)
  sim
}

#' @importFrom stats rbeta
simulate_betareg <- function (object, nsim) {
  ftd <- fitted(object)
  phi <- predict(object,type="precision")
  if(length(ftd) != length(phi)) stop("dimension of precision and mean parameters mismatch")
  rbeta(n = nsim * length(ftd), shape1 = ftd * phi, shape2 = phi - ftd*phi)
}

#' @title Simulation for a beta regression
#' @description simulation method for class \code{IbBetareg}, see \linkS4class{Ib}
#' @param object an object of class \code{IbBetareg}, see \linkS4class{Ib}
#' @param control a \code{list} of parameters for controlling the iterative procedure
#' (see \code{\link{ibControl}}).
#' @param extra \code{NULL} by default; extra parameters to pass to simulation.
#' @param ... further arguments
#' @export
setMethod("simulation", signature = className("betareg","betareg"),
          definition = simulation.betareg)
