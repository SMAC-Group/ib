#' @importFrom lme4 getME mkVarCorr lmer
#' @export
ib.lmerMod <- function(object, thetastart=NULL, control=list(...), Sigma=FALSE, ...){
  # for lme4::lmer:
  # parameters are beta (coefficients),
  # theta (lower tri of Cholesky),
  # and sigma (check ?getME for more informations)
  par_tmp <- getParam(object,Sigma)
  nbeta <- length(par_tmp$beta)
  ntheta <- length(par_tmp$theta)

  # extra parameters to extract variances and correlations
  # see lme4:::VarrCor.merMod
  cnms <- getME(object, "cnms")
  nc <- lengths(cnms)
  fl <- getME(object, "flist")
  nms <- names(fl)[attr(fl, "assign")]

  # initial estimator:
  init_est <- Param_to_Est(par_tmp,Sigma,cnms,nc,nms,all=TRUE)
  pi0 <- init_est$est
  nvar <- init_est$nvar # number of variance components
  ncor <- init_est$ncor # number of correlation components

  # define identifier for different components of (stack) estimator
  p <- length(pi0)
  if(nvar>0) id_var <- c(seq(nbeta+1,nbeta+nvar), p)
  if(ncor>0) id_cor <- seq(nbeta+nvar+1,nbeta+nvar+ncor)


  if(Sigma && p != nbeta + nvar + ncor + 1)
    stop("sigma is assumed to be a scalar")

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

  control <- do.call("ibControl",control)

  # test diff between thetas
  test_theta <- control$tol + 1

  # test at iteration k-1
  tt_old <- test_theta

  # iterator
  k <- 0L

  # create an environment for iterative bootstrap
  env_ib <- new.env(hash=F)
  mf <- model.frame(object)
  names(mf)[1] <- "y"
  assign("data",mf,env_ib)
  cl <- getCall(object)
  cl$data <- quote(data)
  cl$formula[[2]] <- quote(y)
  tmp_object <- object

  # Iterative bootstrap algorithm:
  while(test_theta > control$tol && k < control$maxit){

    # update initial estimator
    if(k>0) par_tmp <- Est_to_Param(t0,Sigma,nbeta,nvar,ncor,nc)
    tmp_object <- setParam(tmp_object,par_tmp)
    sim <- simulation(tmp_object,control)
    tmp_pi <- matrix(NA_real_,nrow=p,ncol=control$H)
    for(h in seq_len(control$H)){
      env_ib$data$y <- sim[,h]
      tmp <- getParam(eval(cl,env_ib),Sigma)
      tmp_pi[,h] <- Param_to_Est(tmp,Sigma,cnms,nc,nms)
    }
    pi_star <- control$func(tmp_pi)

    # update value
    delta <- pi0 - pi_star
    if(Sigma) delta[id_var] <- exp(log(pi0[id_var])-log(pi_star[id_var]))
    t1 <- t0 + delta
    if(ncor>0) t1[id_cor] <- tanh(atanh(t0[id_cor]) + atanh(pi0[id_cor] - atanh(pi_star[id_cor])))

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

  # update lmerMod object
  # TODO: update AIC, BIC, logLik, deviance, residuals
  tmp_object
}

getParam <- function(object, Sigma=FALSE){
  list(beta = getME(object,"beta"),
       theta = if(Sigma) unname(getME(object,"theta")),
       sigma = if(Sigma) sigma(object))
}

# not used at the moment
mkParam <- function(params, nbeta, ntheta){
  if(ntheta > 0 && (length(params) != nbeta + ntheta + 1))
    stop("sigma is assumed to be a scalar")
  list(beta = params[seq_len(nbeta)],
       theta = if(ntheta>0) params[seq(nbeta+1,nbeta+ntheta)],
       sigma = if(ntheta>0) params[length(params)])
}

# inspired from lme4:::setParam
#' @importFrom methods slot<-
setParam <- function(object, params) {
  if(!is.null(params$beta)) {
    slot(object, "beta") <- params$beta
    object@pp$setBeta0(params$beta)
  }
  if(!is.null(params$theta)) {
    slot(object, "theta") <- params$theta
    object@pp$setTheta(params$theta)
  }
  if (!is.null(params$sigma)) {
    snm <- if (object@devcomp$dims[["REML"]])
      "sigmaREML"
    else "sigmaML"
    object@devcomp[["cmp"]][snm] <- params$sigma
  }
  object
}

# Find variances and correlations from theta
Param_to_Est <- function(params, Sigma, cnms, nc, nms, all=FALSE) {
  if(!Sigma && all){
    return(
      list(est = params$beta,
           nvar = 0,
           ncor = 0)
    )
  }
  if(!Sigma && !all){
    return(unlist(params, use.names = FALSE))
  }

  # make a VarCorr object (see ? lme4::VarCorr)
  vc <- mkVarCorr(params$sigma, cnms, nc, params$theta, nms)
  # extract variances components
  vars <- unname(sapply(vc,function(x)diag(x)))
  # extract correlations components
  cors <- unlist(lapply(vc,function(x){
                  tmp <- attr(x,"correlation")
                  if(ncol(tmp)>1) tmp[lower.tri(tmp)]
                }))
  if(all){
    return(
      list(est=c(params$beta, vars, cors, params$sigma),
           nvar = length(vars),
           ncor = length(cors))
    )
  }
  c(params$beta, vars, cors, params$sigma)
}

# Find theta from variances and correlations
# ! Requires Cholesky decomposition
Est_to_Param <- function(est, Sigma, nbeta, nvar, ncor, nc) {
  if(!Sigma){
    return(
      list(beta = est,
           theta = NULL,
           sigma = NULL)
    )
  }

  vars <- est[seq(nbeta+1,nbeta+nvar)]
  if(ncor>0) cors <- est[seq(nbeta+nvar+1,nbeta+nvar+ncor)]
  sc <- est[length(est)]

  # inspired from lme4::mkVarCorr
  ncseq <- seq_along(nc)
  vhl <- split(vars, rep.int(ncseq, nc))
  if(ncor>0) chl <- split(cors, rep.int(ncseq, (nc * (nc - 1))/2))
  ans <- sapply(ncseq, function(i) {
    M <- tcrossprod(sqrt(vhl[[i]]) / sc)
    corr <- diag(nrow = nc[i])
    if(ncor>0) corr[upper.tri(corr)] <- corr[lower.tri(corr)] <- chl[[i]]
    R <- chol(corr * M)
    R[upper.tri(R,diag=TRUE)]
  })

  list(beta = est[seq_len(nbeta)],
       theta = c(ans),
       sigma = sc)
}

# Useful resource for `lme4`
# https://stats.stackexchange.com/a/155500
# https://arxiv.org/pdf/1406.5823.pdf
