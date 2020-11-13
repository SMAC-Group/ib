#' @rdname ib
#' @param var if \code{TRUE}, the variance of the residuals in \code{\link[stats]{lm}} is
#'        also corrected
#' @example /inst/examples/eg_lm.R
#' @seealso \code{\link[stats]{lm}}
#' @importFrom stats lm predict.lm model.matrix
#' @export
ib.lm <- function(object, thetastart=NULL, control=list(...), var = FALSE, ...){
  # initial estimator:
  pi0 <- coef(object)

  if(var) pi0 <- c(pi0, sigma(object))

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

  control <- do.call("ibControl",control)

  # test diff between thetas
  p <- p0 <- length(t0)
  if(var) p0 <- p - 1L
  test_theta <- control$tol + 1

  # test at iteration k-1
  tt_old <- test_theta

  # iterator
  k <- 0L

  # create an environment for iterative bootstrap
  env_ib <- new.env(hash=F)

  # prepare data and formula for fit
  mf <- model.frame(object)
  x <- if(!is.empty.model(object$terms)) model.matrix(object$terms, mf, object$contrasts)
  assign("x",x,env_ib)
  o <- as.vector(model.offset(mf))
  if(!is.null(o)) assign("o",o,env_ib)
  cl <- getCall(object)
  cl$formula <- quote(y~0+x)
  cl$data <- NULL
  # add an offset
  if(!is.null(o)) cl$offset <- quote(o)

  # copy the object
  tmp_object <- object

  if(!var) std <- NULL

  # Iterative bootstrap algorithm:
  while(test_theta > control$tol && k < control$maxit){

    # update initial estimator
    tmp_object$coefficients <- t0[1:p0]
    if(var) std <- t0[p]
    sim <- simulation(tmp_object,control,std)
    tmp_pi <- matrix(NA_real_,nrow=p,ncol=control$H)
    for(h in seq_len(control$H)){
      assign("y",sim[,h],env_ib)
      fit_tmp <- eval(cl,env_ib)
      tmp_pi[1:p0,h] <- coef(fit_tmp)
      if(var) tmp_pi[p,h] <- sigma(fit_tmp)
    }
    pi_star <- control$func(tmp_pi)

    # update value
    delta <- pi0 - pi_star
    if(var) delta[p] <- exp(log(pi0[p])-log(pi_star[p]))
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

  tmp_object$fitted.values <- predict.lm(tmp_object)
  tmp_object$residuals <- unname(model.frame(object))[,1] - tmp_object$fitted.values
  tmp_object$call <- object$call
  tmp_object
}

# inspired from stats::simulate.lm
#' @importFrom stats fitted sigma rnorm runif
simulation.lm <- function(object, control=list(...), std=NULL, ...){
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
