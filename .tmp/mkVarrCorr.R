mkVarCorr <-
function (sc, cnms, nc, theta, nms)
{
  ncseq <- seq_along(nc)
  thl <- split(theta, rep.int(ncseq, (nc * (nc + 1))/2))
  if (!all(nms == names(cnms)))
    warning("nms != names(cnms)  -- whereas lme4-authors thought they were --\n",
            "Please report!", immediate. = TRUE)
  ans <- lapply(ncseq, function(i) {
    Li <- diag(nrow = nc[i])
    Li[lower.tri(Li, diag = TRUE)] <- thl[[i]]
    rownames(Li) <- cnms[[i]]
    val <- tcrossprod(sc * Li)
    stddev <- sqrt(diag(val))
    corr <- t(val/stddev)/stddev
    diag(corr) <- 1
    structure(val, stddev = stddev, correlation = corr)
  })
  if (is.character(nms)) {
    if (anyDuplicated(nms))
      nms <- make.names(nms, unique = TRUE)
    names(ans) <- nms
  }
  structure(ans, sc = sc)
}

VarCorr.merMod <-
function (x, sigma = 1, ...)
{
  if (is.null(cnms <- x@cnms))
    stop("VarCorr methods require reTrms, not just reModule")
  if (missing(sigma))
    sigma <- sigma(x)
  nc <- lengths(cnms)
  structure(mkVarCorr(sigma, cnms = cnms, nc = nc, theta = x@theta,
                      nms = {
                        fl <- x@flist
                        names(fl)[attr(fl, "assign")]
                      }), useSc = as.logical(x@devcomp$dims[["useSc"]]), class = "VarCorr.merMod")
}
