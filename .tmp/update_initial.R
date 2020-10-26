update_initial <- function(object,theta,control){
  object2 <- change_est(object,theta)
  # sim <- simulation(object2,nsim=control$H,seed=control$seed)
  # drop(rowMeans(
  #   apply(data.matrix(sim),2,function(y,obj)update_est(obj,y),obj=object2)
  #   ))
  update_est(object2,control)
}

change_est <- function(object, theta, ...){
  UseMethod("change_est",object)
}

change_est.default <- function(object, theta, ...){
  object$coefficients <- theta
  object
}

update_est <- function(object, control){
  UseMethod("update_est",object)
}

# adapted from stats:::update.defaults
update_est.default <- function(object, control){
  # verify existence of a call object
  if (is.null(call <- getCall(object)))
    stop("need an object with call component")

  # generate new observations
  sim <- simulation(object,nsim=control$H,seed=control$seed)
  # names(sim) <- "y_new"
  y_new <- data.matrix(sim)

  # modify the call
  call$formula <- update.formula(formula(object), y_new~.)

  # create a new environment for evaluation
  # update_env <- new.env(hash=FALSE, parent=env)
  # update_env$y_new <- y_new

  eval(call,envir=current_fn())
  # get_est(eval(call))
  # get_est(eval(call))
  # get_est(update(object,y_new~.))
}
