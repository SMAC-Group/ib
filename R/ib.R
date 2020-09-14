#' @export
ib_control <- function(tol = 1e-5, maxit = 25, verbose = TRUE){
  list(tol=tol,maxit=maxit,verbose=verbose)
}

#' @title Bias  correction via the iterative bootstrap algorithm
#' @author Samuel Orso
#' @
#' @export
ib <- function(object, ibcontrol=ib_control(), method="IB",
               theta_start=NULL, simcontrol=sim_control(), ...){
  # initial estimator:
  pi0 <- coef(object)
  if(!is.null(theta_start)){
    t0 <- theta_start
  } else {
    t0 <- pi0
  }

  # test diff between thetas
  p <- length(t0)
  test_theta <- ibcontrol$tol + 1

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
  while(test_theta > ibcontrol$tol & k < ibcontrol$maxit){
    # update initial estimator
    tmp_object$coefficients <- t0
    sim <- simulation(tmp_object,simcontrol)
    tmp_pi <- matrix(nr=p,nc=simcontrol$H)
    for(h in seq_len(ncol(sim))){
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
    if(ibcontrol$verbose){
      cat("Iteration:",k,"Norm between theta_k and theta_(k-1):",test_theta,"\n")
    }

    # update theta
    t0 <- t1
  }

  tmp_object
}
