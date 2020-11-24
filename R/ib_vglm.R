# These functions are
# Copyright (C) 2020 S. Orso, University of Geneva
# All rights reserved.

#' @rdname ib
#' @seealso \code{\link[VGAM]{vglm}}
#' @importFrom VGAM vglm Coef model.framevlm has.intercept vchol vforsub vbacksub
#' @importFrom methods slot `slot<-` .hasSlot
#' @export
ib.vglm <- function(object, thetastart=NULL, control=list(...), extra_param = FALSE,...){
  # initial estimator:
  pi0 <- Coef(object)

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

  control <- do.call("ibControl",control)

  # test diff between thetas
  p <- p0 <- length(t0)
  test_theta <- control$tol + 1

  # test at iteration k-1
  tt_old <- test_theta

  # iterator
  k <- 0L

  # create an environment for iterative bootstrap
  env_ib <- new.env(hash=F)

  # prepare data and formula for fit
  mf <- model.framevlm(object) # ? problem with mf <- model.frame(object)
  mt <- attr(mf, "terms")
  x <- if(!is.empty.model(mt)) model.matrix(mt, mf, attr(mf,"contrasts"))
  # x <- model.matrixvlm(object)
  # remove intercept from design
  if(!has.intercept(object)) x <- x[,!grepl("Intercept",colnames(x))]
  assign("x",x,env_ib)
  o <- as.vector(model.offset(mf))
  if(!is.null(o)) assign("o",o,env_ib)
  cl <- getCall(object)
  cl$formula <- quote(y~x)
  cl$data <- NULL
  # add an offset
  if(!is.null(o)) cl$offset <- quote(o)
  # FIXME: add support for weights, subset, na.action, start,
  #        etastart, mustart, contrasts, constraints

  # copy the object
  tmp_object <- object

  # Iterative bootstrap algorithm:
  while(test_theta > control$tol && k < control$maxit){

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

    # update increment
    k <- k + 1L

    # test diff between thetas
    test_theta <- sqrt(drop(crossprod(t0-t1)))/p

    # Stop if no more progress
    if(tt_old <= test_theta) {break} else {tt_old <- test_theta}

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
  w <- drop(slot(object, "prior.weights"))
  y <- slot(object, "y")
  M <- NULL # to avoid undefined variable
  eval(slot(fam,"initialize")) # initialize M
  eta <- predict(tmp_object)
  mu <- slot(fam,"linkinv")(eta, extra)
  u <- eval(slot(fam,"deriv"))
  W <- eval(slot(fam,"weight"))
  n <- nrow(x)

  U <- vchol(W, M, n, silent = TRUE)
  tvfor <- vforsub(U, as.matrix(u), M, n)
  res <- vbacksub(U, tvfor, M, n)

  if(.hasSlot(fam, "deviance"))
    tmp_object@criterion$deviance <- slot(fam,"deviance")(mu,y,w,residuals=FALSE,eta,extra)
  if(.hasSlot(fam, "loglikelihood"))
    tmp_object@criterion$loglikelihood <- slot(fam,"loglikelihood")(mu,y,w,residuals=FALSE,eta,extra)

  slot(tmp_object, "predictors") <- eta
  slot(tmp_object, "fitted.values") <- mu
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

