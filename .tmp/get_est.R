get_est <- function(object, ...){
  UseMethod("get_est",object)
}

get_est.default <- function(object, ...){
  coef(object)
}
