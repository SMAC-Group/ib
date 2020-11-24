#' @export
summary.ib <- function(object, ...){
  ## call summary from superclass
  summ <- NextMethod()
  ## add metadata from ib
  summ$iteration <- object$ib$iteration
  summ$of <- object$ib$of
  summ$ib_warn <- object$ib$ib_warn
  class(summ) <- c("summary.ib", class(summ))
  summ
}

#' @export
print.summary.ib <- function(x, ...){
  NextMethod()
  cat("\nIterative bootstrap procedure:")
  cat("\n\n     * number of iterations:", x$iteration)
  cat("\n     * objective function:", x$of)
  if(!is.null(x$ib_warn))
    cat("\n\nWarning while correcting the bias:", x$ib_warn)
}
