# These functions are
# Copyright (C) 2020 S. Orso, University of Geneva
# All rights reserved.

## Define summary method for different classes of the union "SummaryIb"
show.summary.ib <- function(object){
  digits <- max(3L, getOption("digits")) - 3L
  print(object@summ)
  cat("\nIterative bootstrap procedure:")
  cat("\n\n     * number of iterations:", object@ib_extra$iteration)
  cat("\n     * objective function:", format(object@ib_extra$of, digits = digits))
  if(!is.null(object@ib_extra$ib_warn))
    cat("\n\nWarning while correcting the bias:", object@ib_extra$ib_warn)
  invisible(object)
}

#' @title Summarizing a fitted model corrected by the ib procedure
#' @description Method for printing a \code{summary} of
#' class union \linkS4class{SummaryIb}.
#' @param object a summary object of member of \linkS4class{SummaryIb}
#' @seealso \linkS4class{SummaryIb}
#' @export
setMethod("show",
          "SummaryIb",
          definition = show.summary.ib)

## IbBetareg
summaryIbBetareg <- function(object, ...){
  summary.betareg <- getFromNamespace("summary.betareg", ns = "betareg")
  x <- getObject(object)
  y <- getExtra(object)
  summ <- summary.betareg(x, ...)
  new("SummaryIbBetareg",
      summ = summ,
      ib_extra = y)
}

#' @title Summarizing a beta regression fit corrected by
#' the iterative bootstrap
#' @description summary method for class \linkS4class{IbBetareg}
#' @param object an object of class \linkS4class{IbBetareg}
#' @param ... further arguments passed to \code{summary.betareg}
#' @seealso \link[betareg]{summary.betareg}
#' @export
setMethod("summary",
          "IbBetareg",
          definition = summaryIbBetareg)

## IbGlm
#' @importFrom stats summary.glm
summaryIbGlm <- function(object, ...){
  x <- getObject(object)
  y <- getExtra(object)
  summ <- summary.glm(x, ...)
  new("SummaryIbGlm",
      summ = summ,
      ib_extra = y)
}

#' @title Summarizing a Generalized Linear Model regression fit corrected by
#' the iterative bootstrap
#' @description summary method for class \linkS4class{IbGlm}
#' @param object an object of class \linkS4class{IbGlm}
#' @param ... further arguments passed to \code{summary.glm}
#' @seealso \link[stats]{summary.glm}
#' @export
setMethod("summary",
          "IbGlm",
          definition = summaryIbGlm)

## IbGlmer
#' @importFrom utils getFromNamespace
summaryIbGlmer <- function(object, ...){
  summary.glmer <- getFromNamespace("summary.merMod", ns = "lme4")
  x <- getObject(object)
  y <- getExtra(object)
  summ <- summary.glmer(x, ...)
  new("SummaryIbGlmer",
      summ = summ,
      ib_extra = y)
}

#' @title Summarizing a generalized linear mixed model regression fit corrected by
#' the iterative bootstrap
#' @description summary method for class \linkS4class{IbGlmer}
#' @param object an object of class \linkS4class{IbGlmer}
#' @param ... further arguments passed to \code{summary.merMod} of \pkg{lme4}
#' @export
setMethod("summary",
          "IbGlmer",
          definition = summaryIbGlmer)


## IbLm
#' @importFrom stats summary.lm
summaryIbLm <- function(object, ...){
  x <- getObject(object)
  y <- getExtra(object)
  summ <- summary.lm(x, ...)
  new("SummaryIbLm",
      summ = summ,
      ib_extra = y)
}

#' @title Summarizing a linear regression fit corrected by
#' the iterative bootstrap
#' @description summary method for class \linkS4class{IbLm}
#' @param object an object of class \linkS4class{IbLm}
#' @param ... further arguments passed to \code{summary.lm}
#' @seealso \link[stats]{summary.lm}
#' @export
setMethod("summary",
          "IbLm",
          definition = summaryIbLm)

## IbLmer
#' @importFrom utils getFromNamespace
summaryIbLmer <- function(object, ...){
  summary.lmer <- getFromNamespace("summary.merMod", ns = "lme4")
  x <- getObject(object)
  y <- getExtra(object)
  summ <- summary.lmer(x, ...)
  new("SummaryIbLmer",
      summ = summ,
      ib_extra = y)
}

#' @title Summarizing a linear mixed model regression fit corrected by
#' the iterative bootstrap
#' @description summary method for class \linkS4class{IbLmer}
#' @param object an object of class \linkS4class{IbLmer}
#' @param ... further arguments passed to \code{summary.merMod} of \pkg{lme4}
#' @export
setMethod("summary",
          "IbLmer",
          definition = summaryIbLmer)

## IbNegbin
summaryIbNegbin <- function(object, ...){
  summary.negbin <- getFromNamespace("summary.negbin", ns = "MASS")
  x <- getObject(object)
  y <- getExtra(object)
  summ <- summary.negbin(x, ...)
  new("SummaryIbNegbin",
      summ = summ,
      ib_extra = y)
}

#' @title Summarizing a negative binomial regression fits corrected by
#' the iterative bootstrap
#' @description summary method for class \linkS4class{IbNegbin}
#' @param object an object of class \linkS4class{IbNegbin}
#' @param ... further arguments passed to \code{summary.negbin}
#' @seealso \link[MASS]{summary.negbin}
#' @export
setMethod("summary",
          "IbNegbin",
          definition = summaryIbNegbin)

## IbNls
summaryIbNls <- function(object, ...){
  summary.nls <- getFromNamespace("summary.nls", ns = "stats")
  x <- getObject(object)
  y <- getExtra(object)
  summ <- summary.nls(x, ...)
  new("SummaryIbNls",
      summ = summ,
      ib_extra = y)
}

#' @title Summarizing a nonlinear regression fit corrected by
#' the iterative bootstrap
#' @description summary method for class \linkS4class{IbNls}
#' @param object an object of class \linkS4class{IbNls}
#' @param ... further arguments passed to \code{summary.nls} of \pkg{stats}
#' @export
setMethod("summary",
          "IbNls",
          definition = summaryIbNls)

## IbVglm
#' @importFrom VGAM summaryvglm
summaryIbVglm <- function(object, ...){
  x <- getObject(object)
  y <- getExtra(object)
  summ <- summaryvglm(x, ...)
  new("SummaryIbVglm",
      summ = summ,
      ib_extra = y)
}

#' @title Summarizing a vector generalized linear model regression
#' fit corrected by the iterative bootstrap
#' @description summary method for class \linkS4class{IbVglm}
#' @param object an object of class \linkS4class{IbVglm}
#' @param ... further arguments passed to \code{summary.merMod} of \pkg{VGAM}
#' @export
setMethod("summary",
          "IbVglm",
          definition = summaryIbVglm)
