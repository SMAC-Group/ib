ib <- function(object, ...){
  UseMethod("ib")
}

ib.default <- function(object, control, ...){

}

  #starting_value,theta_tilde,X,H,seed,robust=FALSE,pseudo=FALSE,maxIt=100,tol=1e-4,verbose=FALSE){
  # Initial values:
  # current value
  t0 <- starting_value

  # test diff between thetas
  p <- length(t0)
  test_theta <- tol + 1

  # test at iteration k-1
  tt_old <- test_theta

  # iterations
  it <- 0L

  # Iterative bootstrap algorithm:
  while(test_theta > tol & it < maxIt){
    # new pi_star
    new_pi_star <- pi_star(t0,X,H,seed,robust,pseudo)

    # update value
    t1 <- t0 + theta_tilde - new_pi_star

    # update increment
    it <- it + 1L

    # test diff between thetas
    test_theta <- sqrt(crossprod(t0-t1))/p

    # Stop if no more progress
    if(tt_old < test_theta) {break} else {tt_old <- test_theta}

    # update theta
    t0 <- t1
    if(verbose) print(cbind(it,test_theta))
  }

  # Save results
  res <- NULL
  res$est <- t0
  res$test_theta <- test_theta
  res$iter <- it
  return(res)
}
