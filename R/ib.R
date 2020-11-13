#' @title
#' Bias correction via iterative bootstrap
#' @description
#' \code{ib} is used to correct the bias of a fitted model \code{object}
#' with the iterative bootstrap procedure.
#' @param object an \code{object} representing a fitted model (see 'Details').
#' @param thetastart an optional starting value for the iterative procedure.
#' If \code{NULL} (default), the procedure starts at the estimates in \code{object}.
#' @param control a \code{list} of parameters for controlling the iterative procedure
#' (see \code{\link{ibControl}}).
#' @param ... additional optional arguments (see 'Details').
#' @return
#' A fitted model \code{object} where estimates are bias corrected with the \code{ib}.
#' @details
#' The iterative bootstrap procedure is described in
#' \insertCite{kuk1995;textual}{ib} and further
#' studied by \insertCite{guerrier2019;textual}{ib} and
#' \insertCite{guerrier2020;textual}{ib}. The \emph{k}th iteration of this
#' algorithm is
#' \deqn{\hat{\theta}^{k} = \hat{\theta}^{k-1} + \hat{\pi} -
#' \frac{1}{H}\sum_{h=1}^H\hat{\pi}_h(\hat{\theta}^{k-1})}{%
#' \theta^(k)=\theta^(k-1) + \pi - \sum \pi_h(\theta^(k-1)) / H
#' }
#' for \eqn{k=1,2,\ldots} and where the sum is over \eqn{h=1,\ldots,H}.
#' The estimate \eqn{\pi} is provided by the \code{object}.
#' The value \eqn{\pi_h(\theta)} is a parametric bootstrap
#' estimate where the bootstrap sample is generated from \eqn{\theta}
#' and a fixed \code{seed} (see \code{\link{ibControl}}).
#' The greater the parameter value \eqn{H} generally the better bias correction
#' but the more computation it requires (see \code{\link{ibControl}}).
#' If \code{thetastart=NULL}, the initial value of the procedure is \eqn{\theta^(0)=\pi}.
#' The number of iterations are controlled by \code{maxit} parameter of \code{\link{ibControl}}.
#'
#' By default, the method correct \code{\link[stats:coef]{coefficients}} only. For
#' extra parameters, it depends on the model. Currently, \code{ib}
#' supports the following \code{object}:
#' \tabular{ll}{
#'    \code{\link[stats]{glm}}\tab
#'        with \code{shape=TRUE}, the shape parameter for the \code{\link[stats:family]{Gamma}}
#'        family is also corrected. Note that the \code{\link[stats:family]{quasi}} families
#'        are not supported for the moment as they have no simulation method
#'        (see \code{\link[stats]{simulate}}).
#'    \cr
#'    \code{\link[MASS]{glm.nb}}\tab
#'        with \code{overdispersion=TRUE}, the overdispersion parameter of the
#'        negative binomial regression is also corrected.
#'    \cr
#'    \code{\link[stats]{lm}}\tab
#'        with \code{vars=TRUE}, the variance of the residuals is
#'        also corrected. Note that using the \code{ib} is not useful as coefficients
#'        are already unbiased, unless one considers different
#'        data generating mechanism such as censoring, missing values
#'        and outliers (see \code{\link{ibControl}}).
#'    \cr
#'    \code{\link[lme4]{lmer}}\tab
#'        by default, only the fixed effects are corrected. With \code{Sigma=TRUE},
#'        all the random effects (variances and correlations) and the variance
#'        of the residuals are also corrected. Note that using the \code{ib} is
#'        certainly not useful with the argument \code{REML=TRUE} in
#'        \code{\link[lme4]{lmer}} as the bias of variance components is
#'        already addressed, unless one considers different
#'        data generating mechanism such as censoring, missing values
#'        and outliers (see \code{\link{ibControl}}).
#'    \cr
#' }
#' @references
#' \insertAllCited{}
#' @importFrom Rdpack reprompt
#' @export
ib <- function(object, thetastart=NULL, control=list(...), ...){
  UseMethod("ib",object)
}

#' @importFrom stats coef model.matrix getCall predict model.frame is.empty.model model.offset
#' @export
ib.default <- function(object, thetastart=NULL, control=list(...), ...){
  # check control
  control <- do.call("ibControl",control)

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

  # test diff between thetas
  p <- length(t0)
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

  # Iterative bootstrap algorithm:
  while(test_theta > control$tol && k < control$maxit){
    # update initial estimator
    tmp_object$coefficients <- t0
    sim <- simulation(tmp_object,control)
    tmp_pi <- matrix(NA_real_,nrow=p,ncol=control$H)
    for(h in seq_len(control$H)){
      assign("y",sim[,h],env_ib)
      tmp_pi[,h] <- coef(eval(cl,env_ib))
    }
    pi_star <- control$func(tmp_pi)

    # update value
    t1 <- t0 + pi0 - pi_star

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

  tmp_object$fitted.values <- predict(tmp_object)
  tmp_object$residuals <- unname(model.frame(object))[,1] - tmp_object$fitted.values
  tmp_object$call <- object$call
  class(tmp_object) <- c("ib",class(object))
  tmp_object
}

#' @title Auxiliary for controlling IB
#' @description
#' Auxiliary function for \code{\link{ib}} bias correction.
#' @param tol positive convergence tolerance \eqn{\epsilon}.
#' The \code{\link{ib}} procedure converges when
#' \eqn{||\theta^{k+1}-\theta^k||_2/p<\epsilon},
#' where \eqn{p} is the dimension of \eqn{\theta}.
#' @param maxit \code{integer} representing the maximal number of iterations.
#' @param verbose if \code{TRUE}, it prints some output in the console
#' at each iteration.
#' @param seed \code{integer} to set the seed (see \code{\link[base]{Random}}).
#' @param H \code{integer} representing the number of bootstrap estimates
#' (see \code{\link{ib}}).
#' @param cens if \code{TRUE} the simulated responses are censored according to
#' \code{left} and \code{right} values.
#' @param right \code{double} for right-censoring (only used if \code{cens=TRUE}).
#' @param left \code{double} for left-censoring (only used if \code{cens=TRUE}).
#' @param mis if \code{TRUE} the simulated responses have missing data at random.
#' @param prop \code{double} between 0 and 1 representing the proportion of
#' missing data (only used if \code{mis=TRUE}).
#' @param out if \code{TRUE} the simulated responses are also generated with an
#' contamination mechanism
#' @param eps \code{double} between 0 and 1 representing the proportion of
#' outliers in the data (only used if \code{out=TRUE}).
#' @param G a \code{function} to generate outliers. It takes only
#' a sample size as argument.
#' @param func a \code{function} to reduce the \code{H} bootstrap estimates (rowwise).
#' By default, the average is computed. The user can supply a function.
#' One could imagine using other function such as the median or a trimmed mean.
#' @return a list with components named as the arguments.
#' @seealso \code{\link{ib}}, the iterative procedure for bias correction.
#' @export
ibControl <- function(tol = 1e-5, maxit = 25, verbose = FALSE,
                      seed=123L,H=1L,
                      cens=FALSE,right=NULL,left=NULL,
                      mis=FALSE,prop=NULL,
                      out=FALSE,eps=NULL,G=NULL,
                      func=function(x)rowMeans(x,na.rm=T)){
  if(!is.numeric(tol)) stop("`tol` must be numeric")
  if(!is.numeric(maxit)) stop("`maxit` must be numeric")
  if(!is.logical(verbose)) stop("`verbose` must be a boolean")
  if(!is.numeric(seed)) stop("`seed` must be numeric")
  if(!is.numeric(H)) stop("`H` must be numeric")
  if(!is.logical(cens)) stop("`cens` must be a boolean")
  if(!is.logical(mis)) stop("`mis` must be a boolean")
  if(!is.logical(out)) stop("`out` must be a boolean")
  if(!is.function(func)) stop("`func` must be a function")
  list(tol=tol,maxit=maxit,verbose=verbose,
       cens=cens,right=right,left=left,seed=seed,
       H=H,func=func,mis=mis,prop=prop,out=out,eps=eps,G=G)
}
