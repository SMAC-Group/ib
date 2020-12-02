# These functions are
# Copyright (C) 2020 S. Orso, University of Geneva
# All rights reserved.

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
#' @param extra_param if \code{TRUE}, the bias of estimation of extra parameters
#' is performed (see 'Details').
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
#' extra parameters, it depends on the model.
#' @references
#' \insertAllCited{}
#' @importFrom Rdpack reprompt
#' @importFrom stats coef model.matrix getCall predict model.frame is.empty.model model.offset
#' @author Samuel Orso
#' @export
setGeneric("ib",
           function(object, thetastart = NULL, control=list(...), extra_param = FALSE, ...) standardGeneric("ib"),
           signature = "object",
           package = "ib")


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
#' @param out if \code{TRUE} the simulated responses are also generated with a
#' contamination mechanism.
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

setGeneric("simulation",
           function(object, control=list(...), ...) standardGeneric("simulation"),
           signature = "object",
           package = "ib")

#' @importFrom stats simulate
simulation.default <- function(object, control=list(...), ...){
  control <- do.call("ibControl",control)
  sim <- simulate(object,nsim=control$H,seed=control$seed,...)
  if(control$cens) sim <- censoring(sim,control$right,control$left)
  if(control$mis) sim <- missing_at_random(sim, control$prop)
  if(control$out) sim <- outliers(sim, control$eps, control$G)
  data.matrix(sim)
}
