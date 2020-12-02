# right and left are censoring values
censoring <- function(y,right=NULL,left=NULL){
  if(all(!is.null(right),!is.null(left),right<left))
    stop("right-censoring must be greater than left-censoring!")
  if(!is.null(right)) y[y>right] <- right
  if(!is.null(left)) y[y<left] <- left
  y
}

# prop is the proportion of missing data
missing_at_random <- function(y, prop=NULL){
  if(!is.null(prop)){
    if(any(!is.numeric(prop),prop<=0,prop>=1))
      stop("`prop` must be a number between 0 and 1!")
    y[sample.int(length(y),floor(prop*length(y)))] <- NA_real_
  }
  y
}

# eps is the proportion of outliers in the data
# G is generating mechanism for outliers
# G must take only one argument, sample size G(n)
# extra parameters should be included within G
outliers <- function(y, eps=NULL, G=NULL){
  if(all(!is.null(eps),!is.null(G))){
    if(any(!is.numeric(eps),eps<=0,eps>=1))
      stop("`eps` must be a number between 0 and 1!")
    if(!is.function(G)) stop("`G` must be a function!")
    if(length(formals(G))!=1) stop("`G` must have only one argument!")
    n <- length(y)
    m <- floor(eps * n)
    y[sample.int(n,m)] <- G(m)
  }
  y
}
