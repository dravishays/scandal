
### =========================================================================
### Genercis for class ScandalDataSet
### -------------------------------------------------------------------------
###

#' @export
setGeneric("logtpm", function(object, ...) standardGeneric("logtpm"))

#' @export
setGeneric("logtpm<-", function(object, ..., value) standardGeneric("logtpm<-"))

#' @export
setGeneric("childNodes", function(object, ...) standardGeneric("childNodes"))

#' @export
setGeneric("childNodes<-", function(object, ..., value) standardGeneric("childNodes<-"))

#' @export
setGeneric("parentNode", function(object, ...) standardGeneric("parentNode"))

#' @export
setGeneric("parentNode<-", function(object, ..., value) standardGeneric("parentNode<-"))

#' @export
setGeneric("identifier", function(object, ...) standardGeneric("identifier"))

##' @export
#setGeneric("identifier<-", function(object, ..., value) standardGeneric("identifier<-"))

#' @export
setGeneric("complexity", function(object, ..., return_sorted = FALSE) standardGeneric("complexity"))

#' @export
setGeneric("sampleNames", function(object, ..., return_sorted = FALSE) standardGeneric("sampleNames"))

#' @export
setGeneric("meanExpression", function(object, ..., assay_name = NULL, by = "row", log_transform_res = TRUE, genes_subset = NULL) standardGeneric("meanExpression"))

#' @export
setGeneric("unprocessedData", function(object, ...) standardGeneric("unprocessedData"))

#' @export
setGeneric("configParam", function(object, ...) standardGeneric("configParam"))

#' @export
setGeneric("corrMatrix", function(object, ...) standardGeneric("corrMatrix"))

#' @export
setGeneric("corrMatrix<-", function(object, ..., value) standardGeneric("corrMatrix<-"))

### =========================================================================
### Genercis for class ConfigParam
### -------------------------------------------------------------------------
###

#' @export
setGeneric("complexityCutoff", function(x, ...) standardGeneric("complexityCutoff"))

#' @export
setGeneric("complexityCutoff<-", function(x, ..., value) standardGeneric("complexityCutoff<-"))

#' @export
setGeneric("expressionCutoff", function(x, ...) standardGeneric("expressionCutoff"))

#' @export
setGeneric("expressionCutoff<-", function(x, ..., value) standardGeneric("expressionCutoff<-"))

#' @export
setGeneric("housekeepingCutoff", function(x, ...) standardGeneric("housekeepingCutoff"))

#' @export
setGeneric("housekeepingCutoff<-", function(x, ..., value) standardGeneric("housekeepingCutoff<-"))

#' @export
setGeneric("scalingFactor", function(x, ...) standardGeneric("scalingFactor"))

#' @export
setGeneric("scalingFactor<-", function(x, ..., value) standardGeneric("scalingFactor<-"))

#' @export
setGeneric("logBase", function(x, ...) standardGeneric("logBase"))

#' @export
setGeneric("logBase<-", function(x, ..., value) standardGeneric("logBase<-"))
