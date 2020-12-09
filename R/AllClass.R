# These functions are
# Copyright (C) 2020 S. Orso, University of Geneva
# All rights reserved.

## Class definition for the package

## NOTE: due to compatibility issues between S3/S4 classes,
##       we consider supplying the fitted model in a new slot.

#' @describeIn Ib fitted model by \code{glm} from \pkg{stats}
#' @export
setClass("IbGlm",
         slots = list(
           object = "glm",
           ib_extra = "list"))

setOldClass("summary.glm")
#' @describeIn SummaryIb summary of class \code{summary.glm} from \pkg{stats}
#' @export
setClass("SummaryIbGlm",
         slots = list(
           summ = "summary.glm",
           ib_extra = "list"))

#' @describeIn Ib fitted model by \code{lm} from \pkg{stats}
#' @export
setClass("IbLm",
         slots = list(
           object = "lm",
           ib_extra = "list"))

setOldClass("summary.lm")
#' @describeIn SummaryIb summary of class \code{summary.lm} from \pkg{stats}
#' @export
setClass("SummaryIbLm",
         slots = list(
           summ = "summary.lm",
           ib_extra = "list"))

#' @describeIn Ib fitted model by \code{lmer} from \pkg{lme4}
#' @export
setClass("IbLmer",
         slots = list(
           object = className("lmerMod","lme4"),
           ib_extra = "list"))

setOldClass("summary.merMod")
#' @describeIn SummaryIb summary of class \code{summary.merMod} from \pkg{lme4}
#' @export
setClass("SummaryIbLmer",
         slots = list(
           summ = "summary.merMod",
           ib_extra = "list"))

setOldClass("negbin")
#' @describeIn Ib fitted model by \code{glm.nb} from \pkg{MASS}
#' @export
setClass("IbNegbin",
         slots = list(
           object = "negbin",
           ib_extra = "list"))

setOldClass("summary.negbin")
#' @describeIn SummaryIb summary of class \code{summary.negbin} from \pkg{MASS}
#' @export
setClass("SummaryIbNegbin",
         slots = list(
           summ = "summary.negbin",
           ib_extra = "list"))

#' @describeIn Ib fitted model by \code{vglm} from \pkg{VGAM}
#' @export
setClass("IbVglm",
         slots = list(
           object = className("vglm","VGAM"),
           ib_extra = "list"))

#' @describeIn SummaryIb summary of class \code{summary.vglm} from \pkg{VGAM}
#' @export
setClass("SummaryIbVglm",
         slots = list(
           summ = className("summary.vglm","VGAM"),
           ib_extra = "list"))


#' @title
#' An S4 class union for \code{ib}
#' @description
#' Members of the union are \linkS4class{IbGlm}, \linkS4class{IbLm},
#' \linkS4class{IbLmer}, \linkS4class{IbNegbin}, \linkS4class{IbVglm}
#' @details
#' The `Functions` section describes members of the class union.
#' @author Samuel Orso
#' @export
setClassUnion(name = "Ib",
              members = c("IbGlm","IbLm","IbLmer","IbNegbin","IbVglm"))

#' @title An S4 class union for \code{summary}
#' @description
#' Members of the union are \linkS4class{SummaryIbGlm}, \linkS4class{SummaryIbLm},
#' \linkS4class{SummaryIbLmer}, \linkS4class{SummaryIbNegbin}, \linkS4class{SummaryIbVglm}
#' iterative bootstrap procedure
#' @details
#' The `Functions` section describes members of the class union.
#' @author Samuel Orso
#' @export
setClassUnion(name = "SummaryIb",
              members = c("SummaryIbGlm","SummaryIbLm","SummaryIbLmer",
                          "SummaryIbNegbin","SummaryIbVglm"))
