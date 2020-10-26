setParams <-
function (object, params, inplace = FALSE, subset = FALSE)
{
  pNames <- c("beta", "theta")
  if (object@devcomp$dims["useSc"])
    pNames <- c(pNames, "sigma")
  if (!is.list(params) || length(setdiff(names(params), pNames)) >
      0)
    stop("params should be specifed as a list with elements from ",
         "{", paste(shQuote(pNames), collapse = ", "), "}")
  if (!subset && length(setdiff(pNames, names(params))) > 0) {
    warning("some parameters not specified in setParams()")
  }
  nbeta <- length(object@pp$beta(1))
  ntheta <- length(object@pp$theta)
  if (!is.null(beta <- params$beta) && length(beta) != nbeta)
    stop("length mismatch in beta (", length(beta), "!=",
         nbeta, ")")
  if (!is.null(theta <- params$theta) && length(theta) != ntheta)
    stop("length mismatch in theta (", length(theta), "!=",
         ntheta, ")")
  matchNames <- function(x, tn, vecname = "theta") {
    if (!is.null(pn <- names(x))) {
      if (!setequal(pn, tn)) {
        stop("mismatch between ", shQuote(vecname), " parameter vector names and internal names (",
             paste(tn, collapse = ","), ")")
      }
      x <- x[tn]
    }
    else {
      message(vecname, " parameter vector not named: assuming same order as internal vector")
    }
    x
  }
  theta <- matchNames(theta, tnames(object), "theta")
  beta <- matchNames(beta, colnames(getME(object, "X")), "beta")
  sigma <- params$sigma
  if (inplace) {
    stop("modification in place (copy=FALSE) not yet implemented")
  }
  else {
    newObj <- object
    newObj@pp <- newObj@pp$copy()
    newObj@resp <- newObj@resp$copy()
    if (!is.null(beta)) {
      newObj@pp$setBeta0(beta)
      newObj@beta <- beta
    }
    if (!is.null(theta)) {
      newObj@theta <- theta
      newObj@pp$setTheta(theta)
    }
    if (!is.null(sigma)) {
      snm <- if (object@devcomp$dims[["REML"]])
        "sigmaREML"
      else "sigmaML"
      newObj@devcomp[["cmp"]][snm] <- sigma
    }
    return(newObj)
  }
}
