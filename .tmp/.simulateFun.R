.simulateFun <-
function (object, nsim = 1, seed = NULL, use.u = FALSE, re.form = NA,
          ReForm, REForm, REform, newdata = NULL, newparams = NULL,
          formula = NULL, family = NULL, weights = NULL, offset = NULL,
          allow.new.levels = FALSE, na.action = na.pass, cond.sim = TRUE,
          ...)
{
  nullWts <- FALSE
  if (is.null(weights)) {
    if (is.null(newdata))
      weights <- weights(object)
    else {
      nullWts <- TRUE
      weights <- rep(1, nrow(newdata))
    }
  }
  if (missing(object)) {
    if (is.null(formula) || is.null(newdata) || is.null(newparams)) {
      stop("if ", sQuote("object"), " is missing, must specify all of ",
           sQuote("formula"), ", ", sQuote("newdata"), ", and ",
           sQuote("newparams"))
    }
    if (is.character(family))
      family <- get(family, mode = "function", envir = parent.frame())
    if (is.function(family))
      family <- family()
    if (is.null(family) || (family$family == "gaussian" &&
                            family$link == "identity")) {
      lmod <- lFormula(formula, newdata, weights = weights,
                       offset = offset, control = lmerControl(check.formula.LHS = "ignore"))
      devfun <- do.call(mkLmerDevfun, lmod)
      object <- mkMerMod(environment(devfun), opt = list(par = NA,
                                                         fval = NA, conv = NA), lmod$reTrms, fr = lmod$fr)
    }
    else {
      glmod <- glFormula(formula, newdata, family = family,
                         weights = weights, offset = offset, control = glmerControl(check.formula.LHS = "ignore"))
      devfun <- do.call(mkGlmerDevfun, glmod)
      object <- mkMerMod(environment(devfun), opt = list(par = NA,
                                                         fval = NA, conv = NA), glmod$reTrms, fr = glmod$fr)
    }
  }
  stopifnot((nsim <- as.integer(nsim[1])) > 0, is(object, "merMod"))
  if (!is.null(newparams)) {
    object <- setParams(object, newparams)
  }
  re.form.miss <- missing(re.form)
  re.form <- reFormHack(re.form, ReForm, REForm, REform)
  if (!missing(use.u)) {
    if (!re.form.miss) {
      stop("should specify only one of ", sQuote("use.u"),
           " and ", sQuote("re.form"))
    }
    re.form <- if (use.u)
      NULL
    else ~0
  }
  if (is.null(re.form)) {
    re.form <- noLHSform(formula(object))
  }
  if (!is.null(seed))
    set.seed(seed)
  if (!exists(".Random.seed", envir = .GlobalEnv))
    runif(1)
  RNGstate <- .Random.seed
  sigma <- sigma(object)
  etapred <- predict(object, newdata = newdata, re.form = re.form,
                     type = "link", na.action = na.omit, allow.new.levels = allow.new.levels)
  n <- length(etapred)
  makeOp <- function(x, y, op = NULL) {
    if (is.null(op)) {
      substitute(OP(X), list(X = x, OP = y))
    }
    else substitute(OP(X, Y), list(X = x, OP = op, Y = y))
  }
  compReForm <- reOnly(formula(object))
  if (isRE(re.form)) {
    rr <- reOnly(re.form)[[2]]
    ftemplate <- substitute(. ~ . - XX, list(XX = rr))
    compReForm <- update.formula(compReForm, ftemplate)[-2]
  }
  sim.reff <- if (!is.null(findbars(compReForm))) {
    newRE <- mkNewReTrms(object, newdata, compReForm, na.action = na.action,
                         allow.new.levels = allow.new.levels)
    U <- t(newRE$Lambdat %*% newRE$Zt)
    u <- rnorm(ncol(U) * nsim)
    as(U %*% matrix(u, ncol = nsim), "matrix")
  }
  else 0
  val <- if (isLMM(object)) {
    etapred + sigma * (sim.reff + if (cond.sim)
      matrix(rnorm(n * nsim), ncol = nsim)
      else 0)
  }
  else if (isGLMM(object)) {
    etasim <- etapred + sim.reff
    family <- normalizeFamilyName(object@resp$family)
    musim <- family$linkinv(etasim)
    if (family$family == "binomial" && is.matrix(r <- model.response(object@frame))) {
      if (nullWts)
        weights <- rowSums(r)
    }
    if (is.null(sfun <- simfunList[[family$family]])) {
      stop("simulation not implemented for family ", sQuote(family$family))
    }
    if (cond.sim) {
      val <- sfun(object, nsim = 1, ftd = rep_len(musim,
                                                  n * nsim), wts = weights)
    }
    else {
      val <- rep_len(musim, n * nsim)
    }
    if (family$family == "binomial" && is.matrix(r <- model.response(object@frame))) {
      lapply(split(val[[1]], gl(nsim, n, 2 * nsim * n)),
             matrix, ncol = 2, dimnames = list(NULL, colnames(r)))
    }
    else if (family$family == "binomial" && is.factor(val[[1]])) {
      split(val[[1]], gl(nsim, n))
    }
    else split(val, gl(nsim, n))
  }
  else stop("simulate method for NLMMs not yet implemented")
  if (!is.list(val)) {
    dim(val) <- c(n, nsim)
    val <- as.data.frame(val)
  }
  else class(val) <- "data.frame"
  names(val) <- paste("sim", seq_len(nsim), sep = "_")
  f <- fitted(object)
  nm <- names(f)[!is.na(f)]
  if (length(nm) == 0) {
    nm <- as.character(seq(n))
  }
  else if (!is.null(newdata)) {
    nm <- rownames(newdata)
  }
  row.names(val) <- nm
  fit.na.action <- attr(model.frame(object), "na.action")
  if (!missing(na.action) && !is.null(fit.na.action)) {
    class.na.action <- class(attr(na.action(NA), "na.action"))
    if (!identical(class.na.action, class(fit.na.action))) {
      class(fit.na.action) <- class.na.action
    }
  }
  nafun <- function(x) {
    x[] <- apply(x, 2L, napredict, omit = fit.na.action)
    x
  }
  val <- if (is.matrix(val[[1]])) {
    structure(lapply(val, nafun), class = "data.frame")
  }
  else {
    as.data.frame(lapply(val, napredict, omit = fit.na.action))
  }
  nm2 <- if (is.null(newdata))
    names(napredict(na.omit(f), omit = fit.na.action))
  else rownames(napredict(newdata, omit = fit.na.action))
  if (length(nm2) > 0)
    row.names(val) <- nm2
  structure(val, na.action = fit.na.action, seed = RNGstate)
}
