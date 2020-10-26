function (object, nsim)
{
  wts <- object$prior.weights
  if (any(wts != 1))
    message("using weights as shape parameters")
  ftd <- fitted(object)
  shape <- MASS::gamma.shape(object)$alpha * wts
  rgamma(nsim * length(ftd), shape = shape, rate = shape/ftd)
}
