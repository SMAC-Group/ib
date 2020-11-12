#' @rdname ib
#' @param shape if \code{TRUE}, the shape parameter for the \code{\link[stats:family]{Gamma}}
#'        family is also corrected
#' @param overdispersion if \code{TRUE}, the overdispersion parameter of the
#'        negative binomial regression in \code{\link[MASS]{glm.nb}} is also corrected
#' @param var if \code{TRUE}, the variance of the residuals in \code{\link[stats]{lm}} is
#'        also corrected
#' @seealso \code{\link[stats]{glm}}, \code{\link[MASS]{glm.nb}}
#' @example /inst/examples/eg_glm.R
#' @importFrom stats glm predict.glm model.matrix
#' @importFrom MASS gamma.shape
#' @export
ib.glm <- function(object, thetastart=NULL, control=list(...), shape=FALSE, overdispersion=FALSE, var=FALSE, ...){
  # supports only glm.fit currently
  if(object$method != "glm.fit") stop("only implemented for `glm.fit`", call.=FALSE)

  # shape for gamma regression
  if(shape && !grepl("Gamma",object$family$family)) stop("`shape` is for gamma regression", call.=FALSE)

  # overdispersion for negative binomial regression
  if(overdispersion && !grepl("Negative Binomial",object$family$family)) stop("`overdispersion` is for negative binomial regression", call.=FALSE)

  extra <- FALSE
  if(any(shape,overdispersion)) extra <- TRUE

  # initial estimator:
  pi0 <- coef(object)

  if(shape) pi0 <- c(pi0, MASS::gamma.shape(object)$alpha)
  if(overdispersion) pi0 <- c(pi0, object$theta)

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
  if(extra) p0 <- p - 1L
  test_theta <- control$tol + 1

  # test at iteration k-1
  tt_old <- test_theta

  # iterator
  k <- 0L

  # create an environment for iterative bootstrap
  env_ib <- new.env(hash=F)
  assign("x",unname(model.matrix(object)),env_ib)
  cl <- getCall(object)
  cl$formula <- quote(y~0+x)
  cl$data <- NULL
  tmp_object <- object

  if(!shape) shp <- NULL
  if(!var) std <- NULL

  # Iterative bootstrap algorithm:
  while(test_theta > control$tol && k < control$maxit){

    # update initial estimator
    tmp_object$coefficients <- t0[1:p0]
    if(shape) shp <- t0[p]
    if(overdispersion) tmp_object$theta <- t0[p]
    if(var) std <- t0[p]
    sim <- simulation(tmp_object,control,shape=shp,std=std)
    tmp_pi <- matrix(NA_real_,nrow=p,ncol=control$H)
    for(h in seq_len(control$H)){
      assign("y",sim[,h],env_ib)
      fit_tmp <- eval(cl,env_ib)
      tmp_pi[1:p0,h] <- coef(fit_tmp)
      if(shape) tmp_pi[p,h] <- MASS::gamma.shape(fit_tmp)$alpha
      if(overdispersion) tmp_pi[p,h] <- fit_tmp$theta
      if(var) tmp_pi[p,h] <- sigma(fit_tmp)
    }
    pi_star <- control$func(tmp_pi)

    # update value
    delta <- pi0 - pi_star
    if(extra) delta[p] <- exp(log(pi0[p])-log(pi_star[p]))
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

# inspired from stats::simulate.lm
simulation.glm <- function(object, control=list(...), shape=NULL, std=NULL, ...){
  control <- do.call("ibControl",control)

  fam <- object$family$family
  if(fam!="gaussian" && is.null(object$family$simulate)) stop(paste0("simulation not implemented for family ",fam), call.=FALSE)

  set.seed(control$seed)
  if(!exists(".Random.seed", envir = .GlobalEnv)) runif(1)

  if(is.null(std)) std <- sigma(object)

  sim <- switch(fam,
                Gamma = {
                  if(is.null(shape)){
                    matrix(object$family$simulate(object,control$H), ncol=control$H)
                  } else {
                    matrix(simulate_gamma(object,control$H,shape), ncol=control$H)}},
                gaussian = {matrix(fitted(object) + rnorm(length(object$y) * control$H,sd=std), ncol=control$H)},
                matrix(object$family$simulate(object,control$H), ncol=control$H)
  )
  if(control$cens) sim <- censoring(sim,control$right,control$left)
  if(control$mis) sim <- missing_at_random(sim, control$prop)
  if(control$out) sim <- outliers(sim, control$eps, control$G)
  sim
}

# inspired from stats::family::Gamma::simulate
# which does not support "shape" as an argument
#' @importFrom stats rgamma
simulate_gamma <- function (object, nsim, shape){
  wp <- object$prior.weights
  ftd <- fitted(object)
  shp <- shape * wp
  rgamma(nsim * length(ftd), shape = shp, rate = shp/ftd)
}

