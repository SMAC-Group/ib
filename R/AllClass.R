# These functions are
# Copyright (C) 2020 S. Orso, University of Geneva
# All rights reserved.

#' @import methods
#' @importFrom methods setOldClass setClass setClassUnion className
NULL

## ---- Base old-class registrations (S3/S4 interop) -------------------------
setOldClass("betareg")
setOldClass("summary.betareg")
setOldClass("summary.glm")
setOldClass("summary.lm")
setOldClass("negbin")
setOldClass("summary.negbin")
setOldClass("nls")
setOldClass("summary.nls")

## lme4 old classes (helps R CMD check be explicit)
setOldClass("summary.merMod")

## ---- Documentation: Ib* classes (one page, many aliases) ------------------

#' Fitted-model wrappers for the iterative bootstrap (\code{ib})
#'
#' S4 classes that wrap fitted model objects and attach extra metadata from the
#' iterative bootstrap procedure.
#'
#' Members include wrappers around models from \pkg{betareg}, \pkg{stats},
#' \pkg{lme4}, \pkg{MASS}, and \pkg{VGAM}.
#'
#' @slot object  The original fitted model object, corrected by \code{ib}.
#' @slot ib_extra A \code{list} of extra metadata from \code{ib}.
#'
#' @return Objects of the respective S4 classes.
#'
#' @details
#' These classes are light wrappers: they store the corrected fitted object in
#' \code{@object} (see \code{\link{getObject}}) and additional information in
#' \code{@ib_extra} (see \code{\link{getExtra}}).
#'
#' @seealso \code{\link{getObject}}, \code{\link{getExtra}},
#' \linkS4class{Ib}, \linkS4class{SummaryIb}
#'
#' @name Ib-class
#' @aliases
#' Ib-class
#' IbBetareg-class IbGlm-class IbLm-class IbLmer-class IbGlmer-class
#' IbNegbin-class IbNls-class IbVglm-class
NULL

#' @rdname Ib-class
#' @export
setClass("IbBetareg",
         slots = list(object = "betareg", ib_extra = "list"))

#' @rdname Ib-class
#' @export
setClass("IbGlm",
         slots = list(object = "glm", ib_extra = "list"))

#' @rdname Ib-class
#' @export
setClass("IbLm",
         slots = list(object = "lm", ib_extra = "list"))

#' @rdname Ib-class
#' @export
setClass("IbLmer",
         slots = list(object = className("lmerMod", "lme4"),
                      ib_extra = "list"))

#' @rdname Ib-class
#' @export
setClass("IbGlmer",
         slots = list(object = className("glmerMod", "lme4"),
                      ib_extra = "list"))

#' @rdname Ib-class
#' @export
setClass("IbNegbin",
         slots = list(object = "negbin", ib_extra = "list"))

#' @rdname Ib-class
#' @export
setClass("IbNls",
         slots = list(object = "nls", ib_extra = "list"))

#' @rdname Ib-class
#' @export
setClass("IbVglm",
         slots = list(object = className("vglm", "VGAM"),
                      ib_extra = "list"))

## ---- Documentation: SummaryIb* classes (one page, many aliases) -----------

#' Summaries for \code{ib} wrappers
#'
#' S4 classes that wrap summary objects corresponding to the \code{Ib*} classes.
#'
#' @slot summ    The \code{summary.*} object.
#' @slot ib_extra A \code{list} of extra metadata from \code{ib}.
#'
#' @seealso \linkS4class{SummaryIb}, \linkS4class{Ib}
#'
#' @name SummaryIb-class
#' @aliases
#' SummaryIb-class
#' SummaryIbBetareg-class SummaryIbGlm-class SummaryIbLm-class
#' SummaryIbLmer-class SummaryIbGlmer-class SummaryIbNegbin-class
#' SummaryIbNls-class SummaryIbVglm-class
NULL

#' @rdname SummaryIb-class
#' @export
setClass("SummaryIbBetareg",
         slots = list(summ = "summary.betareg", ib_extra = "list"))

#' @rdname SummaryIb-class
#' @export
setClass("SummaryIbGlm",
         slots = list(summ = "summary.glm", ib_extra = "list"))

#' @rdname SummaryIb-class
#' @export
setClass("SummaryIbLm",
         slots = list(summ = "summary.lm", ib_extra = "list"))

#' @rdname SummaryIb-class
#' @export
setClass("SummaryIbLmer",
         slots = list(summ = "summary.merMod", ib_extra = "list"))

#' @rdname SummaryIb-class
#' @export
setClass("SummaryIbGlmer",
         slots = list(summ = "summary.merMod", ib_extra = "list"))

#' @rdname SummaryIb-class
#' @export
setClass("SummaryIbNegbin",
         slots = list(summ = "summary.negbin", ib_extra = "list"))

#' @rdname SummaryIb-class
#' @export
setClass("SummaryIbNls",
         slots = list(summ = "summary.nls", ib_extra = "list"))

#' @rdname SummaryIb-class
#' @export
setClass("SummaryIbVglm",
         slots = list(summ = className("summary.vglm", "VGAM"),
                      ib_extra = "list"))

## ---- Class unions (keep their own topics) ---------------------------------

#' An S4 class union for \code{Ib}
#'
#' Members: \linkS4class{IbBetareg}, \linkS4class{IbGlm},
#' \linkS4class{IbLm}, \linkS4class{IbLmer}, \linkS4class{IbGlmer},
#' \linkS4class{IbNegbin}, \linkS4class{IbNls}, \linkS4class{IbVglm}.
#'
#' @name Ib-class
#' @aliases Ib
#' @export
setClassUnion(name = "Ib",
              members = c("IbBetareg","IbGlm","IbLm","IbLmer","IbGlmer",
                          "IbNegbin","IbNls","IbVglm"))

#' An S4 class union for \code{SummaryIb}
#'
#' Members: \linkS4class{SummaryIbBetareg}, \linkS4class{SummaryIbGlm},
#' \linkS4class{SummaryIbLm}, \linkS4class{SummaryIbLmer},
#' \linkS4class{SummaryIbGlmer}, \linkS4class{SummaryIbNegbin},
#' \linkS4class{SummaryIbNls}, \linkS4class{SummaryIbVglm}.
#'
#' @name SummaryIb-class
#' @aliases SummaryIb
#' @export
setClassUnion(name = "SummaryIb",
              members = c("SummaryIbBetareg","SummaryIbGlm","SummaryIbLm",
                          "SummaryIbLmer","SummaryIbGlmer",
                          "SummaryIbNegbin","SummaryIbNls","SummaryIbVglm"))
