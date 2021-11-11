# These functions are
# Copyright (C) 2020 S. Orso, University of Geneva
# All rights reserved.

#' @importFrom stats glm predict.glm model.matrix model.frame model.offset is.empty.model terms
#' @importFrom MASS gamma.shape
#' @importFrom methods new
ib.glm <- function(object, thetastart=NULL, control=list(...), extra_param = FALSE, ...){
  # supports only glm.fit currently
  if(object$method != "glm.fit") stop("only implemented for `glm.fit`", call.=FALSE)

  # controls
  control <- do.call("ibControl",control)

  # initial estimator:
  # regression coefficients
  pi0 <- coef(object)
  p0 <- length(pi0)

  # extra parameters
  fam <- object$family$family
  # adjust for negbin
  isNegbin <- FALSE
  if(grepl("Negative Binomial",fam)){
    fam <- "negbin"
    isNegbin <- TRUE
  }

  if(extra_param){
    pi0 <- switch(fam,
                  gaussian = {c(pi0, sigma(object))},
                  Gamma = {c(pi0, gamma.shape(object)$alpha)},
                  negbin = {c(pi0, object$theta)}
    )
    if(is.null(pi0))
      stop(gettextf("extra_param for family '%s' is not implemented", fam), domain = NA)
  }
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

  extra <- NULL

  # Iterative bootstrap algorithm:
  while(test_theta > control$tol && k < control$maxit){
    # update initial estimator
    tmp_object$coefficients <- t0[1:p0]
    if(extra_param) switch (fam,
                            Gamma = {extra <- t0[p]},
                            gaussian = {extra <- t0[p]},
                            negbin = {tmp_object$theta <- t0[p]})
    sim <- simulation(tmp_object,control,extra)
    tmp_pi <- matrix(NA_real_,nrow=p,ncol=control$H)
    for(h in seq_len(control$H)){
      assign("y",sim[,h],env_ib)
      fit_tmp <- tryCatch(error = function(cnd) NULL, {eval(cl,env_ib)})
      if(is.null(fit_tmp)) next
      tmp_pi[1:p0,h] <- coef(fit_tmp)
      if(extra_param)
        tmp_pi[p,h] <- switch(fam,
                              Gamma = {gamma.shape(fit_tmp)$alpha},
                              gaussian = {sigma(fit_tmp)},
                              negbin = {fit_tmp$theta})
    }
    pi_star <- control$func(tmp_pi)

    # update value
    delta <- pi0 - pi_star
    t1 <- t0 + delta
    if(extra_param && control$constraint) t1[p] <- exp(log(t0[p]) + log(pi0[p]) - log(pi_star[p]))

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

  # update glm object
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

  if(isNegbin){
    return(
      new("IbNegbin",
          object = tmp_object,
          ib_extra = ib_extra)
    )
  }

  new("IbGlm",
      object = tmp_object,
      ib_extra = ib_extra)
}

#' @rdname ib
#' @details
#' For \link[stats]{glm}, if \code{extra_param=TRUE}: the shape parameter for the
#' \code{\link[stats:family]{Gamma}}, the variance of the residuals in \code{\link[stats]{lm}} or
#' the overdispersion parameter of the negative binomial regression in \code{\link[MASS]{glm.nb}},
#' are also corrected. Note that the \code{\link[stats:family]{quasi}} families
#' are not supported for the moment as they have no simulation method
#' (see \code{\link[stats]{simulate}}). Bias correction for extra parameters
#' of the \code{\link[stats:family]{inverse.gaussian}} is not yet implemented.
#' @seealso \code{\link[stats]{glm}}, \code{\link[MASS]{glm.nb}}
#' @example /inst/examples/eg_glm.R
#' @export
setMethod("ib", className("glm", "stats"),
          definition = ib.glm)

# inspired from stats::simulate.lm
simulation.glm <- function(object, control=list(...), extra=NULL, ...){
  control <- do.call("ibControl",control)

  fam <- object$family$family
  if(fam!="gaussian" && is.null(object$family$simulate))
    stop(gettextf("simulation not implemented for family '%s'",fam),
         call.=FALSE, domain=NA)

  if(grepl("Negative Binomial",fam)) fam <- "negbin"

  set.seed(control$seed)
  if(!exists(".Random.seed", envir = .GlobalEnv)) runif(1)

  # user-defined simulation method
  if(!is.null(control$sim)){
    sim <- control$sim(object, control, extra, ...)
    return(sim)
  }

  sim <- switch(fam,
                Gamma = {
                  if(is.null(extra)){
                    matrix(object$family$simulate(object,control$H), ncol=control$H)
                  } else {
                    matrix(simulate_gamma(object,control$H,extra), ncol=control$H)}},
                gaussian = if(is.null(extra)){
                  matrix(fitted(object) + rnorm(length(object$y) * control$H, sd=sigma(object)), ncol=control$H)
                } else {
                  matrix(fitted(object) + rnorm(length(object$y) * control$H, sd=extra), ncol=control$H)
                },
                negbin = {matrix(simulate_negbin(object,control$H), ncol=control$H)},
                matrix(object$family$simulate(object,control$H), ncol=control$H)
  )
  if(control$cens) sim <- censoring(sim,control$right,control$left)
  if(control$mis) sim <- missing_at_random(sim, control$prop)
  if(control$out) sim <- outliers(sim, control$eps, control$G)
  sim
}

#' @title Simulation for a Generalized Linear Model regression
#' @description simulation method for class \linkS4class{IbGlm}
#' @param object an object of class \linkS4class{IbGlm}
#' @param control a \code{list} of parameters for controlling the iterative procedure
#' (see \code{\link{ibControl}}).
#' @param extra \code{NULL} by default; extra parameters to pass to simulation.
#' @param ... further arguments
#' @export
setMethod("simulation", signature = className("glm","stats"),
          definition = simulation.glm)

# inspired from stats::family::Gamma::simulate
# which does not support "shape" as an argument
#' @importFrom stats rgamma
simulate_gamma <- function (object, nsim, shape){
  if(shape<0) stop("'shape' must be positive")
  wp <- object$prior.weights
  ftd <- fitted(object)
  shp <- shape * wp
  rgamma(n = nsim * length(ftd), shape = shp, rate = shp/ftd)
}

# ib.negbin (MASS)
ib.negbin <- ib.glm

#' \code{\link{ib}} method for \code{negbin} object
#' from \code{\link[MASS]{glm.nb}} function of \pkg{MASS}
#' package.
#' @inheritParams ib,glm-method
#' @export
setMethod("ib", signature = "negbin",
          definition = ib.negbin)

simulation.negbin <- simulation.glm

# inspired from MASS::simulate.negbin
#' @importFrom MASS rnegbin
simulate_negbin <- function (object, nsim) {
  if(object$theta<0) stop("'theta' must be positive")
  ftd <- fitted(object)
  rnegbin(n = nsim * length(ftd), mu = ftd, theta = object$theta)
}

#' @title Simulation for a negative binomial regression
#' @description simulation method for class \linkS4class{IbNegbin}
#' @param object an object of class \linkS4class{IbNegbin}
#' @param control a \code{list} of parameters for controlling the iterative procedure
#' (see \code{\link{ibControl}}).
#' @param extra \code{NULL} by default; extra parameters to pass to simulation.
#' @param ... further arguments
#' @export
setMethod("simulation", signature = "negbin",
          definition = simulation.negbin)
