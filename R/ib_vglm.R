# These functions are
# Copyright (C) 2020 S. Orso, University of Geneva
# All rights reserved.

#' @rdname ib
#' @details
#' For \link[VGAM]{vglm}, \code{extra_param} is currently not used.
#' Indeed, the philosophy of a vector generalized linear model is to
#' potentially model all parameters of a distribution with a linear predictor.
#' Hence, what would be considered as an extra parameter in \code{\link[stats]{glm}}
#' for instance, may already be captured by the default \code{coefficients}.
#' However, correcting the bias of a \code{coefficients} does not imply
#' that the bias of the parameter of the distribution is corrected
#' (by \href{https://en.wikipedia.org/wiki/Jensen's_inequality}{Jensen's inequality}),
#' so we may use this feature in a future version of the package.
#' Note that we currently only support distributions
#' with a \code{simslot} (see \code{\link[VGAM]{simulate.vlm}}).
#' @example /inst/examples/eg_vglm.R
#' @seealso \code{\link[VGAM]{vglm}}
#' @importFrom VGAM Coef has.intercept model.framevlm predictvglm vbacksub vchol vglm vforsub
#' @importFrom methods slot `slot<-` .hasSlot
#' @importFrom stats model.weights
#' @export
ib.vglm <- function(object, thetastart=NULL, control=list(...), extra_param = FALSE,...){
  # controls
  control <- do.call("ibControl",control)

  # initial estimator:
  pi0 <- coef(object)

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
  p <- p0 <- length(t0)
  test_theta <- control$tol + 1

  # iterator
  k <- 0L

  # create an environment for iterative bootstrap
  env_ib <- new.env(hash=F)

  # prepare data and formula for fit
  cl <- getCall(object)
  intercept_only <- cl$formula[[3]] == 1 # check for intercept only models
  # alternatively: intercept_only <- object@misc$intercept.only
  mf <- model.framevlm(object) # ? problem with mf <- model.frame(object)
  mt <- attr(mf, "terms")
  if(!intercept_only){
    x <- if(!is.empty.model(mt)) model.matrix(mt, mf, attr(mf,"contrasts"))
    # x <- model.matrixvlm(object)
    # remove intercept from design
    if(has.intercept(object)) x <- x[,!grepl("Intercept",colnames(x))]
    assign("x",x,env_ib)
    cl$formula <- quote(y~x)
  } else {
    cl$formula <- quote(y~1)
  }
  cl$data <- NULL
  o <- as.vector(model.offset(mf))
  if(!is.null(o)) assign("o",o,env_ib)
  # add an offset
  if(!is.null(o)) cl$offset <- quote(o)
  # FIXME: add support for subset, na.action, start,
  #        etastart, mustart, contrasts, constraints
  n <- nrow(mf)
  w <- model.weights(mf)
  if(!length(w)){
    w <- rep_len(1, n)
    } else {
    cl$weights <- quote(w)
    assign("w",w,env_ib)
  }
  if(is.null(cl$etastart)) etastart <- NULL

  # copy the object
  tmp_object <- object

  # Iterative bootstrap algorithm:
  while(test_theta > control$tol && k < control$maxit){
    # browser()
    # update initial estimator
    slot(tmp_object, "coefficients") <- t0[1:p0]
    sim <- simulation(tmp_object,control)
    tmp_pi <- matrix(NA_real_,nrow=p,ncol=control$H)
    for(h in seq_len(control$H)){
      assign("y",sim[,h],env_ib)
      # FIXME: deal with warnings from vglm.fitter
      fit_tmp <- eval(cl,env_ib)
      tmp_pi[1:p0,h] <- coef(fit_tmp)
    }
    pi_star <- control$func(tmp_pi)

    # update value
    delta <- pi0 - pi_star
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

  # update vglm object
  extra <- slot(object, "extra")
  fam <- slot(object, "family")
  y <- slot(object, "y")
  M <- NULL # to avoid undefined variable
  # w <- c(slot(object, "prior.weights"))
  # w <- drop(weights(object, "prior"))
  # if(!length(w)==0) w <- rep_len(1,n)
  eval(slot(fam,"initialize")) # initialize different parameters (among which M)
  eta <- predictvglm(tmp_object)
  mu <- slot(fam,"linkinv")(eta, extra)
  u <- eval(slot(fam,"deriv"))
  W <- eval(slot(fam,"weight"))

  U <- vchol(W, M, n, silent = TRUE)
  tvfor <- vforsub(U, as.matrix(u), M, n)
  res <- vbacksub(U, tvfor, M, n)

  if(.hasSlot(fam, "deviance") && !is.null(body(slot(fam,"deviance"))))
    tmp_object@criterion$deviance <- slot(fam,"deviance")(mu,y,w,residuals=FALSE,eta,extra)
  if(.hasSlot(fam, "loglikelihood") && !is.null(body(slot(fam,"loglikelihood"))))
    tmp_object@criterion$loglikelihood <- slot(fam,"loglikelihood")(mu,y,w,residuals=FALSE,eta,extra)

  slot(tmp_object, "predictors") <- eta
  slot(tmp_object, "fitted.values") <- as.matrix(mu)
  slot(tmp_object, "residuals") <- res
  slot(tmp_object, "call") <- slot(object,"call")

  tmp_object
}

# inspired from VGAM::simulate.vlm
#' @importFrom VGAM familyname
simulation.vglm <- function(object, control=list(...), extra_param = NULL, ...){
  control <- do.call("ibControl",control)

  fam <- slot(object, "family")
  if(is.null(body(slot(fam, "simslot"))))
    stop(paste0("simulation not implemented for family ", familyname(object)), call.=FALSE)

  set.seed(control$seed)
  if(!exists(".Random.seed", envir = .GlobalEnv)) runif(1)

  sim <- matrix(slot(fam, "simslot")(object,control$H), ncol = control$H)

  if(control$cens) sim <- censoring(sim,control$right,control$left)
  if(control$mis) sim <- missing_at_random(sim, control$prop)
  if(control$out) sim <- outliers(sim, control$eps, control$G)
  sim
}
