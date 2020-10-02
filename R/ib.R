#' @title Bias correction via iterative bootstrap
#' @description Correct the bias of a fitted model object with the iterative
#' bootstrap procedure.
#' @param object an object representing a fitted model.
#' @param ... additional optional arguments.
#' @return The \code{object} with bias corrected estimators.
#' @export
ib <- function(object, thetastart=NULL, control=list(...), ...){
  UseMethod("ib",object)
}

#' @export
ib.default <- function(object, thetastart=NULL, control=list(...), ...){
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

  control <- do.call("ibControl",control)

  # test diff between thetas
  p <- length(t0)
  test_theta <- control$tol + 1

  # test at iteration k-1
  tt_old <- test_theta

  # iterator
  k <- 0L

  # create an environment for iterative bootstrap
  env_ib <- new.env(hash=F)
  assign("x",unname(model.matrix(object))[,-1],env_ib)
  cl <- call(as.character(getCall(object)[[1]]),quote(y~x))
  tmp_object <- object

  # Iterative bootstrap algorithm:
  while(test_theta > control$tol & k < control$maxit){
    # update initial estimator
    tmp_object$coefficients <- t0
    sim <- simulation(tmp_object,control)
    tmp_pi <- matrix(NA_real_,nrow=p,ncol=control$H)
    for(h in seq_len(control$H)){
      assign("y",sim[,h],env_ib)
      tmp_pi[,h] <- coef(eval(cl,env_ib))
    }
    pi_star <- rowMeans(tmp_pi)

    # update value
    t1 <- t0 + pi0 - pi_star

    # update increment
    k <- k + 1L

    # test diff between thetas
    test_theta <- sqrt(crossprod(t0-t1))/p

    # Stop if no more progress
    if(tt_old <= test_theta) {break} else {tt_old <- test_theta}

    # Print info
    if(control$verbose){
      cat("Iteration:",k,"Norm between theta_k and theta_(k-1):",test_theta,"\n")
    }

    # update theta
    t0 <- t1
  }

  tmp_object$fitted.values <- predict.lm(tmp_object)
  tmp_object$residuals <- unname(model.frame(object))[,1] - tmp_object$fitted.values
  tmp_object$call <- object$call
  tmp_object
}

#' @title Auxiliary for controlling IB
#' @export
ibControl <- function(tol = 1e-5, maxit = 25, verbose = FALSE,
                      cens=FALSE,right=NULL,left=NULL,seed=123L,H=1L){
  if(!is.logical(cens)) stop("cens must a boolean")
  list(tol=tol,maxit=maxit,verbose=verbose,
       cens=cens,right=right,left=left,seed=seed,H=H)
}
