#' @include AllGenerics.R
NULL


### =========================================================================
### ScandalDataSet objects (start)
### -------------------------------------------------------------------------
###

setClassUnion("ScandalDataSetOrNULL", c("NULL"))
setClassUnion("MatrixOrNULL", c("Matrix", "matrix", "NULL"))

#'
#' @title ScandalDataSet class
#' @description An S4 class for storing single-cell seqeuncing data, analysis
#'
#' @slot childNodes
#' @slot parentNode
#' @slot unprocessedData
#' @slot preprocConfig
#' @slot nodeID
#' @slot projectID
#'
#' @method
#'
#' @author Avishay Spitzer
#'
#' @export
#' @exportClass ScandalDataSet
setClass("ScandalDataSet",
         slots = c(childNodes = "SimpleList",
                   parentNode = "ScandalDataSetOrNULL",
                   unprocessedData = "MatrixOrNULL",
                   preprocConfig = "PreprocConfig",
                   nodeID = "character",
                   projectID = "character"),
         contains = "SingleCellExperiment"
)

setIs("ScandalDataSet", "ScandalDataSetOrNULL")

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Constructor
###

#'
#' @title
#'
#' @description
#'
#' @param ...
#' @param childNodes
#' @param parentNode
#' @param preprocConfig
#' @param nodeID
#' @param projectID
#'
#' @details
#'
#' @return A ScandalDataSet object is returned from the constructor
#'
#' @author Avishay Spitzer
#'
#' @export
ScandalDataSet <- function(..., childNodes = S4Vectors::SimpleList(), parentNode = NULL, preprocConfig = DefaultPreprocConfig(), nodeID = NODE_ID(), projectID = PROJ_ID()) {

  sce <- SingleCellExperiment::SingleCellExperiment(...)

  if(!is(sce, "SingleCellExperiment")) {
    sce <- as(sce, "SingleCellExperiment")
  }

  if (is.null(parentNode))
    unprocessedData <- assay(sce)
  else
    unprocessedData <- NULL

  object <- new("ScandalDataSet", sce,
                childNodes = childNodes,
                parentNode = parentNode,
                unprocessedData = unprocessedData,
                preprocConfig = preprocConfig,
                nodeID = nodeID,
                projectID = projectID)

  int_colData(object)$Scandal <- S4Vectors::DataFrame(row.names = colnames(object))
  int_elementMetadata(object)$Scandal <- S4Vectors::DataFrame(row.names = rownames(object))
  int_metadata(object)$Scandal <- list()

  return (object)
}

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Validity
###

setValidity("ScandalDataSet", function(object) {

  is_child_sds <- sapply(childNodes(object), function(c) is(c, "ScandalDataSet"))

  if (!(base::all(is_child_sds) == TRUE))
    return (sprintf("Every child node must be a ScandalDataSet object"))

  return (TRUE)
})

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Getters and setters.
###

#' @include AllGenerics.R
#' @export
setMethod("logtpm", "ScandalDataSet",   function(object, ...) {
  return (assay(object, i = "logtpm", ...))
})

#' @include AllGenerics.R
#' @export
setReplaceMethod("logtpm", c("ScandalDataSet", "ANY"),   function(object, ..., value) {
  assay(object, i = "logtpm", ...) <- value
  return (object)
})

#' @export
setMethod("childNodes", "ScandalDataSet", function(object) {
  return(object@childNodes)
})

#' @export
setReplaceMethod("childNodes", "ScandalDataSet", function(object, value) {
  object@childNodes <- value
  return(object)
})

#' @export
setMethod("parentNode", "ScandalDataSet", function(object) {
  return(object@parentNode)
})

#' @export
setReplaceMethod("parentNode", "ScandalDataSet", function(object, value) {
  object@parentNode <- value
  return(object)
})

#' @export
setMethod("nodeID", "ScandalDataSet", function(object) {
  return(object@nodeID)
})

#' @export
setMethod("projectID", "ScandalDataSet", function(object) {
  return(object@projectID)
})

#' @export
setMethod("sampleNames", "ScandalDataSet", function(object) {
  return(unique(gsub("*-.*", "", colnames(object))))
})

#' @export
setMethod("unprocessedData", "ScandalDataSet", function(object) {

  o <- object
  while(!is.null(parentNode(o)))
    o <- parentNode(o)

  sname <- base::unique(.cell2tumor(colnames(object)))

  # return only the cells that belong to this specific tumor
  return (o@unprocessedData[, .subset_cells(colnames(o@unprocessedData), sname), drop = FALSE])
})

#' @export
setMethod("preprocConfig", "ScandalDataSet", function(object) {
  return(object@preprocConfig)
})

##' @export
#setReplaceMethod("preprocConfig", "ScandalDataSet", function(x, value) {
#  x@preprocConfig <- value
#  return(x)
#})

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Display
###

scat <- function(fmt, vals=character(), exdent=2, ...) {
  vals <- ifelse(nzchar(vals), vals, "''")
  lbls <- paste(S4Vectors:::selectSome(vals), collapse=" ")
  txt <- sprintf(fmt, length(vals), lbls)
  cat(strwrap(txt, exdent=exdent, ...), sep="\n")
}

#' @export
setMethod("show", "ScandalDataSet", function(object) {
  callNextMethod()
  scat("childNodes(%d): %s\n", names(childNodes(object)))
  cat("parentNode:", ifelse(is.null(parentNode(object)), "None", nodeID(parentNode(object))), "\n")
  cat("nodeID:", nodeID(object), "\n")
  cat("projectID:", projectID(object), "\n")
  cat("unprocessedData:", class(unprocessedData(object)), "with", NROW(unprocessedData(object)), "rows and", NCOL(unprocessedData(object)), "columns\n")
})

is_scandal_object <- function(object) { return (is.null(object) | !is(object, "ScandalDataSet")) }

is_valid_assay <- function(x) { return (!(is.null(x)) & (is(x, "Matrix") | is.matrix(x))) }

# Generates a random experiment ID by sampling a random integer
NODE_ID <- function() { paste0("NODE", base::sample(1:1e9, 1, replace = FALSE)) }
PROJ_ID <- function() { paste0("PROJ", base::sample(1:1e9, 1, replace = FALSE)) }

### -------------------------------------------------------------------------
### ScandalDataSet objects (end)
### =========================================================================
###

### =========================================================================
### PreprocConfig objects (start)
### -------------------------------------------------------------------------
###

#' @export
setClass("PreprocConfig",
         slots = c(complexity_cutoff = "vector",
                   expression_cutoff = "numeric",
                   housekeeping_cutoff = "numeric",
                   scaling_factor = "numeric",
                   log_base = "numeric"))

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Constructors
###

#' @export
PreprocConfig <- function(complexity_cutoff, expression_cutoff, housekeeping_cutoff, log_base, scaling_factor) {
  cp <- new("PreprocConfig")

  cp@complexity_cutoff <- complexity_cutoff
  cp@expression_cutoff <- expression_cutoff
  cp@housekeeping_cutoff <- housekeeping_cutoff
  cp@log_base <- log_base
  cp@scaling_factor <- scaling_factor

  return (cp)
}

#' @export
NucseqPreprocConfig <- function() {
  return (PreprocConfig(complexity_cutoff = c(2000, 6000), expression_cutoff = 5, housekeeping_cutoff = 7, log_base = 2, scaling_factor = 10))
}

#' @export
SS2PreprocConfig <- function() {
  return (PreprocConfig(complexity_cutoff = c(3000, 8000), expression_cutoff = 4, housekeeping_cutoff = 7, log_base = 2, scaling_factor = 10))
}

#' @export
DefaultPreprocConfig <- function() {
  return (NucseqPreprocConfig())
}

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Validity
###

setValidity("PreprocConfig", function(object) {
  return(TRUE)
})

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Getters and setters.
###

#' @export
setMethod("complexityCutoff", "PreprocConfig", function(x) {
  return(x@complexity_cutoff)
})

#' @export
setReplaceMethod("complexityCutoff", "PreprocConfig", function(x, value) {
  x@complexity_cutoff <- value
  return(x)
})

#' @export
setMethod("expressionCutoff", "PreprocConfig", function(x) {
  return(x@expression_cutoff)
})

#' @export
setReplaceMethod("expressionCutoff", "PreprocConfig", function(x, value) {
  x@expression_cutoff <- value
  return(x)
})

#' @export
setMethod("housekeepingCutoff", "PreprocConfig", function(x) {
  return(x@housekeeping_cutoff)
})

#' @export
setReplaceMethod("housekeepingCutoff", "PreprocConfig", function(x, value) {
  x@housekeeping_cutoff <- value
  return(x)
})

#' @export
setMethod("scalingFactor", "PreprocConfig", function(x) {
  return(x@scaling_factor)
})

#' @export
setReplaceMethod("scalingFactor", "PreprocConfig", function(x, value) {
  x@scaling_factor <- value
  return(x)
})

#' @export
setMethod("logBase", "PreprocConfig", function(x) {
  return(x@log_base)
})

#' @export
setReplaceMethod("logBase", "PreprocConfig", function(x, value) {
  x@log_base <- value
  return(x)
})

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Display
###

#' @export
setMethod("show", "PreprocConfig", function(object) {
  cat("Complexity cutoff:", complexityCutoff(object), "\n")
  cat("Expression cutoff:", expressionCutoff(object), "\n")
  cat("Housekeeping cutoff:", housekeepingCutoff(object), "\n")
  cat("Scaling factor:", scalingFactor(object), "\n")
  cat("Log base:", logBase(object), "\n")
})

is_config_object <- function(object) { return (!is.null(object) & is(object, "PreprocConfig")) }

### -------------------------------------------------------------------------
### PreprocConfig objects (end)
### =========================================================================
###
