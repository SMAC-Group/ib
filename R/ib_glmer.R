# These functions are
# Copyright (C) 2024 S. Orso, University of Geneva
# All rights reserved.

ib.glmerMod <- function(object, thetastart=NULL, control=list(...), extra_param = FALSE, ...){
  # controls
  control <- do.call("ibControl",control)

  # for lme4::glmer:
  # parameters are beta (coefficients),
  # theta (lower tri of Cholesky),
  # and sigma (check ?lme4::getME for more informations)
  par_tmp <- getParam(object,extra_param)
  nbeta <- length(par_tmp$beta)
  ntheta <- length(par_tmp$theta)

  # extra parameters to extract variances and correlations
  # see lme4:::VarrCor.merMod
  cnms <- getME(object, "cnms")
  nc <- lengths(cnms)
  fl <- getME(object, "flist")
  nms <- names(fl)[attr(fl, "assign")]

  # initial estimator:
  init_est <- Param_to_Est(par_tmp,extra_param,cnms,nc,nms,all=TRUE)
  pi0 <- init_est$est
  nvar <- init_est$nvar # number of variance components
  ncor <- init_est$ncor # number of correlation components

  # define identifier for different components of (stack) estimator
  p <- length(pi0)
  if(nvar>0) id_var <- c(seq(nbeta+1,nbeta+nvar), p)
  if(ncor>0) id_cor <- seq(nbeta+nvar+1,nbeta+nvar+ncor)


  # if(extra_param && p != nbeta + nvar + ncor + 1)
    # stop("sigma is assumed to be a scalar")

  if(!is.null(thetastart)){
    if(is.numeric(thetastart) && length(thetastart) == p){
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
  mf <- model.frame(object)
  names(mf)[1] <- "y"
  isWEIGHTS <- names(mf) %in% "(weights)"
  cl <- getCall(object)
  cl$data <- quote(data)
  if(any(isWEIGHTS)) {
    weights_name <- as.character(cl$weights)
    names(mf)[isWEIGHTS] <- weights_name
  }
  assign("data",mf,env_ib)
  if(length(cl$formula)==1) cl$formula <- get(paste(cl$formula)) # get formula
  cl$formula[[2]] <- quote(y)
  # FIXME: add support for subset, na.action, start, offset,
  #        contrasts

  # FIXME: We need a deep copy (see ?lme4::modular and ? methods::ReferenceClasses):
  # With the following line, we keep modifying the original object
  # tmp_object <- object
  # temporary solution by new evaluation:
  tmp_object <- eval(cl,env_ib)

  # initial value
  diff <- rep(NA_real_, control$maxit)

  # Iterative bootstrap algorithm:
  while(test_theta > control$tol && k < control$maxit){

    # update initial estimator
    if(k>0) par_tmp <- Est_to_Param(t0,extra_param,nbeta,nvar,ncor,nc)
    tmp_object <- setParam(tmp_object,par_tmp)
    sim <- simulation(tmp_object,control)
    tmp_pi <- matrix(NA_real_,nrow=p,ncol=control$H)
    for(h in seq_len(control$H)){
      env_ib$data$y <- sim[,h]
      # FIXME: deal with warnings from checkConv,
      # see https://stats.stackexchange.com/questions/110004/how-scared-should-we-be-about-convergence-warnings-in-lme4
      tmp <- getParam(eval(cl,env_ib),extra_param)
      tmp_pi[,h] <- Param_to_Est(tmp,extra_param,cnms,nc,nms)
    }
    pi_star <- control$func(tmp_pi)

    # update value
    delta <- pi0 - pi_star
    t1 <- t0 + delta
    if(extra_param && control$constraint) t1[id_var] <- exp(log(t0[id_var]) + log(pi0[id_var])-log(pi_star[id_var]))
    if(ncor>0 && control$constraint) t1[id_cor] <- tanh(atanh(t0[id_cor]) + atanh(pi0[id_cor] - atanh(pi_star[id_cor])))

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

  # update glmerMod object
  updateLmer(tmp_object, extra_param)

  # additional metadata
  ib_extra <- list(
    iteration = k,
    of = sqrt(drop(crossprod(delta))),
    estimate = t0,
    test_theta = test_theta,
    boot = tmp_pi)

  new("IbGlmer",
      object = tmp_object,
      ib_extra = ib_extra)
}

#' @rdname ib
#' @details
#' For \code{\link[lme4]{glmer}}, by default, only the fixed effects are corrected.
#' If \code{extra_param=TRUE}: all the random effects
#' (variances and correlations) are also corrected.
#' @example /inst/examples/eg_glmer.R
#' @seealso \code{\link[lme4]{glmer}}
#' @importFrom lme4 glmer
#' @export
setMethod("ib", signature = className("glmerMod", "lme4"),
          definition = ib.glmerMod)

getParam2 <- function(object, Sigma=FALSE){
  list(beta = getME(object,"beta"),
       theta = if(Sigma) unname(getME(object,"theta")))
}

simulation.glmerMod <- simulation.default

#' @title Simulation for linear mixed model regression
#' @description simulation method for class \code{IbGlmer}, see \linkS4class{Ib}
#' @param object an object of class \code{IbGlmer}, see \linkS4class{Ib}
#' @param control a \code{list} of parameters for controlling the iterative procedure
#' (see \code{\link{ibControl}}).
#' @param ... further arguments.
#' @export
setMethod("simulation", signature = className("glmerMod","lme4"),
          definition = simulation.glmerMod)

# Useful resources for `lme4`:
# * https://arxiv.org/pdf/1406.5823.pdf
# * https://stats.stackexchange.com/a/155500
