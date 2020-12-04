# These functions are
# Copyright (C) 2020 S. Orso, University of Geneva
# All rights reserved.

## Class definition for the package

#' @title
#' An S4 class to represent a glm object to which ib is applied
#'
#' @slot object modified \link[stats]{glm} with bias corrected quantities
#' @slot ib_extra a \code{list} with additional information from the
#' iterative bootstrap procedure
#' @author Samuel Orso
#' @export
setClass("IbGlm",
         slots = list(
           object = "glm",
           ib_extra = "list"
         ),
         contains = "glm")

#' @title
#' An S4 class to represent a lm object to which ib is applied
#'
#' @slot object modified from \link[stats]{lm} with bias corrected quantities
#' @slot ib_extra a \code{list} with additional information from the
#' iterative bootstrap procedure
#' @author Samuel Orso
#' @export
setClass("IbLm",
         slots = list(
           object = "lm",
           ib_extra = "list"
         ),
         contains = "lm")

#' @title
#' An S4 class to represent a lmerMod object to which ib is applied
#'
#' @slot object modified from \link[lme4]{lmer} with bias corrected quantities
#' @slot ib_extra a \code{list} with additional information from the
#' iterative bootstrap procedure
#' @author Samuel Orso
#' @export
setClass("IbLmer",
         slots = list(
           object = "lmerMod",
           ib_extra = "list"
         ),
         contains = "lmerMod")

setOldClass("negbin")
#' @title
#' An S4 class to represent a negbin object to which ib is applied
#'
#' @slot object modified from \link[MASS]{glm.nb} with bias corrected quantities
#' @slot ib_extra a \code{list} with additional information from the
#' iterative bootstrap procedure
#' @author Samuel Orso
#' @export
setClass("IbNegbin",
         slots = list(
           object = "negbin",
           ib_extra = "list"
         ),
         contains = "negbin")

#' @title
#' An S4 class to represent a vglm object to which ib is applied
#'
#' @slot object a \link[VGAM]{vglm} object with bias corrected quantities
#' @slot ib_extra a \code{list} with additional information from the
#' iterative bootstrap procedure
#' @author Samuel Orso
#' @export
setClass("IbVglm",
         slots = list(
           object = "vglm",
           ib_extra = "list"
         ),
         contains = "vglm")


#' @title
#' An S4 class union
#'
#' @description
#' Members of the union are \linkS4class{IbGlm}, \linkS4class{IbLm},
#' \linkS4class{IbLmer}, \linkS4class{IbNegbin}, \linkS4class{IbVglm}
#'
#' @export
setClassUnion(name = "Ib",
              members = c("IbGlm","IbLm","IbLmer","IbNegbin","IbVglm"))
