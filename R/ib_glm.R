#' @importFrom stats glm predict.glm model.matrix
#' @export
ib.glm <- function(object, thetastart=NULL, control=list(...), shape=FALSE, ...){
  # supports only glm.fit currently
  if(object$method != "glm.fit") stop("only implemented for `glm.fit`", call.=FALSE)

  # shape for gamma regression
  if(!is.null(shape) && object$family$family != "Gamma") stop("`shape` is for gamma regression", call.=FALSE)

  # initial estimator:
  pi0 <- coef(object)

  if(shape) pi0 <- c(pi0, MASS::gamma.shape(object)$alpha)

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
  if(shape) p0 <- p - 1L
  test_theta <- control$tol + 1

  # test at iteration k-1
  tt_old <- test_theta

  # iterator
  k <- 0L

  # create an environment for iterative bootstrap
  env_ib <- new.env(hash=F)
  assign("x",unname(model.matrix(object))[,-1],env_ib)
  cl <- getCall(object)
  cl[[2]] <- quote(y~x)
  tmp_object <- object

  if(!shape) shp <- NULL

  # Iterative bootstrap algorithm:
  while(test_theta > control$tol & k < control$maxit){

    # update initial estimator
    tmp_object$coefficients <- t0[1:p0]
    if(shape) shp <- t0[p]
    sim <- simulation(tmp_object,control,shp)
    tmp_pi <- matrix(NA_real_,nrow=p,ncol=control$H)
    for(h in seq_len(control$H)){
      assign("y",sim[,h],env_ib)
      fit_tmp <- eval(cl,env_ib)
      tmp_pi[1:p0,h] <- coef(fit_tmp)
      if(shape) tmp_pi[p,h] <- MASS::gamma.shape(fit_tmp)$alpha
    }
    pi_star <- rowMeans(tmp_pi)

    # update value
    delta <- pi0 - pi_star
    if(shape) delta[p] <- exp(log(pi0[p])-log(pi_star[p]))
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
  tmp_object
}
