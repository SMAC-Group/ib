bootMer <-
function (x, FUN, nsim = 1, seed = NULL, use.u = FALSE, re.form = NA,
          type = c("parametric", "semiparametric"), verbose = FALSE,
          .progress = "none", PBargs = list(), parallel = c("no", "multicore",
                                                            "snow"), ncpus = getOption("boot.ncpus", 1L), cl = NULL)
{
  stopifnot((nsim <- as.integer(nsim[1])) > 0)
  if (.progress != "none") {
    pbfun <- get(paste0(.progress, "ProgressBar"))
    setpbfun <- get(paste0("set", .simpleCap(.progress),
                           "ProgressBar"))
    pb <- do.call(pbfun, PBargs)
  }
  do_parallel <- have_mc <- have_snow <- NULL
  eval(initialize.parallel)
  if (do_parallel && .progress != "none")
    message("progress bar disabled for parallel operations")
  FUN <- match.fun(FUN)
  type <- match.arg(type)
  if (!is.null(seed))
    set.seed(seed)
  else if (!exists(".Random.seed", envir = .GlobalEnv))
    runif(1)
  mc <- match.call()
  t0 <- FUN(x)
  if (!is.numeric(t0))
    stop("bootMer currently only handles functions that return numeric vectors")
  mle <- list(beta = getME(x, "beta"), theta = getME(x, "theta"))
  if (isLMM(x))
    mle <- c(mle, list(sigma = sigma(x)))
  if (type == "parametric") {
    argList <- list(x, nsim = nsim, na.action = na.exclude)
    if (!missing(re.form)) {
      argList <- c(argList, list(re.form = re.form))
    }
    else {
      argList <- c(argList, list(use.u = use.u))
    }
    ss <- do.call(simulate, argList)
  }
  else {
    if (!missing(re.form))
      stop(paste(sQuote("re.form")), "cannot be used with semiparametric bootstrapping")
    if (use.u) {
      if (isGLMM(x))
        warning("semiparametric bootstrapping is questionable for GLMMs")
      ss <- replicate(nsim, fitted(x) + sample(residuals(x,
                                                         "response"), replace = TRUE), simplify = FALSE)
    }
    else {
      stop("semiparametric bootstrapping with use.u=FALSE not yet implemented")
    }
  }
  ffun <- local({
    FUN
    refit
    x
    ss
    verbose
    do_parallel
    length.t0 <- length(t0)
    f1 <- factory(function(i) FUN(refit(x, ss[[i]])), errval = rep(NA,
                                                                   length.t0))
    function(i) {
      ret <- f1(i)
      if (verbose) {
        cat(sprintf("%5d :", i))
        str(ret)
      }
      if (!do_parallel && .progress != "none") {
        setpbfun(pb, i/nsim)
      }
      ret
    }
  })
  simvec <- seq_len(nsim)
  res <- if (do_parallel) {
    if (have_mc) {
      parallel::mclapply(simvec, ffun, mc.cores = ncpus)
    }
    else if (have_snow) {
      if (is.null(cl)) {
        cl <- parallel::makePSOCKcluster(rep("localhost",
                                             ncpus))
        parallel::clusterExport(cl, varlist = getNamespaceExports("lme4"))
        if (RNGkind()[1L] == "L'Ecuyer-CMRG")
          parallel::clusterSetRNGStream(cl)
        res <- parallel::parLapply(cl, simvec, ffun)
        parallel::stopCluster(cl)
        res
      }
      else parallel::parLapply(cl, simvec, ffun)
    }
  }
  else lapply(simvec, ffun)
  t.star <- do.call(cbind, res)
  rownames(t.star) <- names(t0)
  msgs <- list()
  for (mtype in paste0("factory-", c("message", "warning",
                                     "error"))) {
    msgs[[mtype]] <- trimws(unlist(lapply(res, attr, mtype)))
    msgs[[mtype]] <- table(msgs[[mtype]])
  }
  if ((numFail <- sum(msgs[["factory-error"]])) > 0) {
    warning("some bootstrap runs failed (", numFail, "/",
            nsim, ")")
  }
  fail.msgs <- if (numFail == 0)
    NULL
  else msgs[["factory-error"]]
  s <- structure(list(t0 = t0, t = t(t.star), R = nsim, data = model.frame(x),
                      seed = .Random.seed, statistic = FUN, sim = "parametric",
                      call = mc, ran.gen = "simulate(<lmerMod>, 1, *)", mle = mle),
                 class = c("bootMer", "boot"))
  attr(s, "bootFail") <- numFail
  attr(s, "boot.fail.msgs") <- fail.msgs
  attr(s, "boot.all.msgs") <- msgs
  attr(s, "boot_type") <- "boot"
  s
}
