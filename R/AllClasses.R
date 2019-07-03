
### =========================================================================
### ScandalDataSet objects (start)
### -------------------------------------------------------------------------
###

setClassUnion("ScandalDataSetOrNULL", c("ScandalDataSet", "NULL"))
setClassUnion("MatrixOrNULL", c("Matrix", "matrix", "NULL"))

#'
#' @description An S4 class for storing single-cell seqeuncing data, analysis
#'
#' @slot childNodes
#' @slot parentNode
#'
#' @author Avishay Spitzer
#'
#' @export
setClass("ScandalDataSet",
         slots = c(childNodes = "SimpleList",
                   parentNode = "ScandalDataSetOrNULL",
                   unprocessedData = "MatrixOrNULL"),
         contains = "SingleCellExperiment"
)

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
#' @param identifier
#' @param preprocConfig
#'
#' @details
#'
#' @return A ScandalDataSet object is returned from the constructor
#'
#' @author Avishay Spitzer
#'
#' @export
ScandalDataSet <- function(..., childNodes = SimpleList(), parentNode = NULL, identifier = exp_id(), preprocConfig = DefaultPreprocConfig()) {
  sce <- SingleCellExperiment(...)
  if(!is(sce, "SingleCellExperiment")) {
    sce <- as(sce, "SingleCellExperiment")
  }

  object <- new("ScandalDataSet", sce)
  object@childNodes <- childNodes
  object@parentNode <- parentNode

  int_colData(object)$Scandal <- S4Vectors::DataFrame(row.names = colnames(object))
  int_elementMetadata(object)$Scandal <- S4Vectors::DataFrame(row.names = rownames(object))
  int_metaData(object)$Scandal <- list()

  int_metaData(object)$Scandal[["identifier"]] <- identifier
  int_metaData(object)$Scandal[["preprocConfig"]] <- preprocConfig

  if (is.null(object@parentNode))
    object@unprocessedData <- assay(sce)
  else
    object@unprocessedData <- NULL

  return (object)
}

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Validity
###

setValidity2("ScandalDataSet", function(object) {

  if (is.null(int_colData(object)$Scandal))
    return (sprintf("int_colData$Scandal cannot be set to NULL"))
  if (!is(int_colData(object)$Scandal, "DataFrame"))
    return (sprintf("int_colData$Scandal must be a DataFrame object and not ", class(int_colData$Scandal)))

  if (is.null(int_elementMetadata(object)$Scandal))
    return (sprintf("int_elementMetadata$Scandal cannot be set to NULL"))
  if (!is(int_elementMetadata(object)$Scandal, "DataFrame"))
    return (sprintf("int_elementMetadata$Scandal must be a DataFrame object and not ", class(int_elementMetadata$Scandal)))

  if (is.null(int_metaData(object)$Scandal))
    return (sprintf("int_metaData$Scandal cannot be set to NULL"))
  if (!is(int_metaData(object)$Scandal, "list"))
    return (sprintf("int_metaData$Scandal must be a list object and not ", class(int_metaData$Scandal)))

  for(c in object@childNodes)
    if (!is(c, "ScandalDataSet"))
      return (sprintf("Every child node must be a ScandalDataSet object and not ", class(c)))

  if (!is.character(identifier(object)))
    return (sprintf("A ScandalDataSet object must have a character identifier"))

  if (!is(preprocConfig(object), "PreprocConfig"))
    return (sprintf("A ScandalDataSet object must have a preprocessing configuration object of class PreprocConfig and not ", class(preprocConfig(object))))

  return (TRUE)
})

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Getters and setters.
###

#' @export
setMethod("logtpm", "ScandalDataSet",   function(object, ...) {
  return (assay(object, i = "logtpm", ...))
})

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
setMethod("identifier", "ScandalDataSet", function(object) {
  return(int_metaData(object)$Scandal[["identifier"]])
})

#setReplaceMethod("identifier", "ExtendedSingleCellExperiment", function(x, value) {
#  x@identifier <- value
#  return(x)
#})

#' @export
setMethod("sampleNames", "ScandalDataSet", function(object) {
  return(unique(gsub("*-.*", "", colnames(object))))
})

#' @export
setMethod("unprocessedData", "ScandalDataSet", function(object) {

  o <- object
  while(!is.null(parentNode(o)))
    o <- parentNode(o)

  return (o@unprocessedData)
})

#' @export
setMethod("preprocConfig", "ScandalDataSet", function(object) {
  return(int_metaData(object)$Scandal[["preprocConfig"]])
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
  scat("indSampleData(%d): %s\n", names(indSampleData(object)))
  cat("identifier:", identifier(object), "\n")
  cat("unprocessedData:", class(unprocessedData(object)), "with", NROW(unprocessedData(object)), "rows and", NCOL(unprocessedData(object)), "columns\n")
})

is_scandal_object <- function(object) { return (is.null(object) | !is(object, "ScandalDataSet")) }

is_valid_assay <- function(x) { return (!(is.null(x)) & (is(x, "Matrix") | is.matrix(x)) & is.numeric(x)) }

# Generates a random experiment ID by sampling a random integer
exp_id <- function() { paste0("EXP", base::sample(1:1e9, 1, replace = FALSE)) }

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

setValidity2("PreprocConfig", function(object) {
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

is_config_object <- function(object) { return (is.null(object) | !is(object, "PreprocConfig")) }

### -------------------------------------------------------------------------
### PreprocConfig objects (end)
### =========================================================================
###
