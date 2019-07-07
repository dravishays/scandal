
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
setGeneric("nodeID", function(object, ...) standardGeneric("nodeID"))

#' @export
setGeneric("projectID", function(object, ...) standardGeneric("projectID"))

#' @export
setGeneric("sampleNames", function(object, ..., return_sorted = FALSE) standardGeneric("sampleNames"))

#' @export
setGeneric("unprocessedData", function(object, ...) standardGeneric("unprocessedData"))

#' @export
setGeneric("preprocConfig", function(object, ...) standardGeneric("preprocConfig"))

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
