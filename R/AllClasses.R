
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
#' @slot identifier
#' @slot unprocessedData
#' @slot configParam
#' @slot corrMatrix
#'
#' @author Avishay Spitzer
#'
#' @export
ScandalDataSet <- setClass("ScandalDataSet",
                           slots = c(childNodes = "SimpleList",
                                     parentNode = "ScandalDataSetOrNULL",
                                     identifier = "character",
                                     unprocessedData = "MatrixOrNULL",
                                     configParam = "ConfigParam",
                                     corrMatrix = "MatrixOrNULL"),
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
#' @param
#'
#' @details
#'
#' @return A ScandalDataSet object is returned from the constructor
#'
#' @author Avishay Spitzer
#'
#' @export
ScandalDataSet <- function(..., childNodes = SimpleList(), parentNode = NULL, identifier = exp_id(), configParam = DefaultConfigParam()) {
  sce <- SingleCellExperiment(...)
  if(!is(sce, "SingleCellExperiment")) {
    sce <- as(sce, "SingleCellExperiment")
  }

  object <- new("ScandalDataSet", sce)
  childNodes(object) <- childNodes
  object@parentNode <- parentNode
  object@identifier <- identifier
  object@unprocessedData <- assay(sce)
  object@configParam <- configParam
  object@corrMatrix <- NULL

  SingleCellExperiment::int_colData(object)$Scandal <- S4Vectors::DataFrame(row.names = colnames(object))

  return (object)
}

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Validity
###

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Getters and setters.
###

#' @export
setMethod("logtpm", "ExtendedSingleCellExperiment",   function(object, ...) {
  return (assay(object, i = "logtpm", ...))
})

#' @export
setReplaceMethod("logtpm", c("ExtendedSingleCellExperiment", "ANY"),   function(object, ..., value) {
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
  return(object@identifier)
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
  return (object@unprocessedData)
})

#' @export
setMethod("configParam", "ScandalDataSet", function(object) {
  return(object@configParam)
})

#' @export
setMethod("corrMatrix", "ScandalDataSet", function(object) {
  if (is.null(object@corrMatrix))
    return (NULL)

  return(object@corrMatrix)
})

#' @export
setReplaceMethod("corrMatrix", "ScandalDataSet", function(object, value) {

  if(is.null(value)) {
    object@corrMatrix <- NULL
    return (object)
  }

  stopifnot(NCOL(value) == NROW(value))
  stopifnot(NCOL(value) == NCOL(object))

  object@corrMatrix <- value

  return(object)
})

##' @export
#setReplaceMethod("configParam", "ScandalDataSet", function(x, value) {
#  x@configParam <- value
#  return(x)
#})

#' @export
setMethod("complexity", "ScandalDataSet", function(object, return_sorted = FALSE) {

  # Complexity is calculated using the unprocessed (un-centered) data
  c <- base::apply(unprocessedData(object), 2, function(x) base::sum(object != 0))

  if(isTRUE(return_sorted)) {
    c <- base::sort(c)
  }

  return(c)
})

#' @export
setMethod("meanExpression", "ScandalDataSet", function(object, assay_name = NULL, by = "row", log_transform_res = TRUE, genes_subset = NULL) {

  if (is.null(assay_name))
    m <- assay(object)
  else if (assay_name == "unprocessed")
    m <- unprocessedData(object)
  else if (assay_name %in% assayNames(object))
    m <- assay(object, i = assay_name)
  else
    stop(paste0("assay ", assay_name, " does not exist"))

  if (!(by %in% c("row", "col")))
    stop(paste0("by can be either row col, got ", by))

  if (!is.null(genes_subset))
    m <- m[rownames(m) %in% genes_subset, ]

  m <- apply(m, ifelse(by == "row", 1, 2), mean)

  if (isTRUE(log_transform_res))
    m <- log2(m + 1)

  return (m)
})

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

  if (!is.null(corrMatrix(object)))
    cat("corrMatrix:", class(corrMatrix(object)), "with", NROW(corrMatrix(object)), "rows and", NCOL(corrMatrix(object)), "columns\n")
  else
    cat("corrMatrix:\n")
})

is_scandal_object <- function(object) { return (is.null(object) | !is(object, "ScandalDataSet")) }

# Generates a random experiment ID by sampling a random integer
exp_id <- function() { paste0("EXP", base::sample(1:1e9, 1, replace = FALSE)) }

### -------------------------------------------------------------------------
### ScandalDataSet objects (end)
### =========================================================================
###

#' @export
ConfigParam <- setClass("ConfigParam",
                        slots = c(complexity_cutoff = "vector",
                                  expression_cutoff = "numeric",
                                  housekeeping_cutoff = "numeric",
                                  scaling_factor = "numeric",
                                  log_base = "numeric"))

#' @export
ConfigParam <- function(complexity_cutoff, expression_cutoff, housekeeping_cutoff, log_base, scaling_factor) {
  cp <- new("ConfigParam")

  cp@complexity_cutoff <- complexity_cutoff
  cp@expression_cutoff <- expression_cutoff
  cp@housekeeping_cutoff <- housekeeping_cutoff
  cp@log_base <- log_base
  cp@scaling_factor <- scaling_factor

  return (cp)
}

#' @export
NucseqConfigParam <- function() {
  return (ConfigParam(complexity_cutoff = c(2000, 6000), expression_cutoff = 5, housekeeping_cutoff = 7, log_base = 2, scaling_factor = 10))
}

#' @export
SS2ConfigParam <- function() {
  return (ConfigParam(complexity_cutoff = c(3000, 8000), expression_cutoff = 4, housekeeping_cutoff = 7, log_base = 2, scaling_factor = 10))
}

#' @export
DefaultConfigParam <- function() {
  return (NucseqConfigParam())
}

#' @export
setMethod("complexityCutoff", "ConfigParam", function(x) {
  return(x@complexity_cutoff)
})

#' @export
setReplaceMethod("complexityCutoff", "ConfigParam", function(x, value) {
  x@complexity_cutoff <- value
  return(x)
})

#' @export
setMethod("expressionCutoff", "ConfigParam", function(x) {
  return(x@expression_cutoff)
})

#' @export
setReplaceMethod("expressionCutoff", "ConfigParam", function(x, value) {
  x@expression_cutoff <- value
  return(x)
})

#' @export
setMethod("housekeepingCutoff", "ConfigParam", function(x) {
  return(x@housekeeping_cutoff)
})

#' @export
setReplaceMethod("housekeepingCutoff", "ConfigParam", function(x, value) {
  x@housekeeping_cutoff <- value
  return(x)
})

#' @export
setMethod("scalingFactor", "ConfigParam", function(x) {
  return(x@scaling_factor)
})

#' @export
setReplaceMethod("scalingFactor", "ConfigParam", function(x, value) {
  x@scaling_factor <- value
  return(x)
})

#' @export
setMethod("logBase", "ConfigParam", function(x) {
  return(x@log_base)
})

#' @export
setReplaceMethod("logBase", "ConfigParam", function(x, value) {
  x@log_base <- value
  return(x)
})

#' @export
setMethod("show", "ConfigParam", function(object) {
  cat("Complexity cutoff:", complexityCutoff(object), "\n")
  cat("Expression cutoff:", expressionCutoff(object), "\n")
  cat("Housekeeping cutoff:", housekeepingCutoff(object), "\n")
  cat("Scaling factor:", scalingFactor(object), "\n")
  cat("Log base:", logBase(object), "\n")
})

is_config_object <- function(object) { return (is.null(object) | !is(object, "ConfigParam")) }
