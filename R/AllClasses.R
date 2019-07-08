
#' @title S4 classes definition
#' @details Includes the generics definition file
#' @include AllGenerics.R
NULL # Do not remove me!!!

### =========================================================================
### PreprocConfig objects (start)
### -------------------------------------------------------------------------
###

#'
#' @title PreprocConfig class
#'
#' @description An S4 class for storing preprocessing configuration parameters.
#'
#' @slot complexityCutoff A numeric vector of length 2 representing the lower and
#' upper bounds of complexity (i.e. the number of detected genes per cell).
#' @slot expressionCutoff A numeric representing the minimal log2 mean expression
#' per gene below which a gene is considered lowly expressed.
#' @slot housekeepingCutoff A numeric representing the log2 mean expression of
#' house-keeping genes (i.e. genes that are highly expressed in all cells) per
#' cell below which a cell is considered low quality.
#' @slot logBase A numeric representing the logarithm base for performing log
#' transformation on the data.
#' @slot scalingFactor A numeric representing a scaling factor by which to divide
#' each data point before log transformation.
#' @slot pseudoCount A numeric representing the pseudo count added when performing
#' log transformation to avoid taking the log of zero.
#' @slot typeMatrix A logical indicating if the dataset should be represented using
#' the S4 Matrix class (instead of base R matrix) to reduce memory overhead using
#' sparse matrix representation.
#'
#' @section Methods:
#' \describe{
#'   \item{\code{complexityCutoff}}{Getter/setter for the complexity cutoff}
#'   \item{\code{expressionCutoff}}{Getter/setter for the expression cutoff}
#'   \item{\code{housekeepingCutoff}}{Getter/setter for the housekeeping cutoff}
#'   \item{\code{logBase}}{Getter/setter for the log base}
#'   \item{\code{scalingFactor}}{Getter/setter for the scaling factor}
#'   \item{\code{pseudoCount}}{Getter/setter for the pseudo count}
#'   \item{\code{typeMatrix}}{Getter/setter for the Matrix type}
#' }
#'
#' @examples
#' pc <- PreprocConfig(complexityCutoff = c(0, 10000), expressionCutoff = 5, housekeepingCutoff = 7, logBase = 2, scalingFactor = 10, pseudoCount = 1, typeMatrix = TRUE)
#'
#' @aliases PreprocConfig
#'
#' @author Avishay Spitzer
#'
#' @export
setClass("PreprocConfig",
         slots = c(complexityCutoff = "vector",
                   expressionCutoff = "numeric",
                   housekeepingCutoff = "numeric",
                   logBase = "numeric",
                   scalingFactor = "numeric",
                   pseudoCount = "numeric",
                   typeMatrix = "logical"))

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Constructors
###

#'
#' @describeIn PreprocConfig-class Constructs a new \code{PreprocConfig} object.
#'
#' @export
PreprocConfig <- function(complexityCutoff, expressionCutoff, housekeepingCutoff, logBase, scalingFactor, pseudoCount, typeMatrix) {
  cp <- new("PreprocConfig",
            complexityCutoff = complexityCutoff,
            expressionCutoff = expressionCutoff,
            housekeepingCutoff = housekeepingCutoff,
            logBase = logBase,
            scalingFactor = scalingFactor,
            pseudoCount = pseudoCount,
            typeMatrix = typeMatrix)

  return (cp)
}

#'
#' @describeIn PreprocConfig-class constructs default preprocessing configuration for data
#' generated from frozen samples using SmartSeq2 protocol (single-nuclei sequencing).
#'
#' @export
NucseqPreprocConfig <- function() {
  return (PreprocConfig(complexityCutoff = c(2000, 6000), expressionCutoff = 5, housekeepingCutoff = 7, logBase = 2, scalingFactor = 10, pseudoCount = 1, typeMatrix = TRUE))
}

#'
#' @describeIn PreprocConfig-class constructs default preprocessing configuration for data
#' generated from fresh samples using SmartSeq2 protocol.
#'
#' @export
SS2PreprocConfig <- function() {
  return (PreprocConfig(complexityCutoff = c(3000, 8000), expressionCutoff = 4, housekeepingCutoff = 7, logBase = 2, scalingFactor = 10, pseudoCount = 1, typeMatrix = TRUE))
}

#'
#' @describeIn PreprocConfig-class constructs default preprocessing configuration
#' (currently trhe selected default configuration is Nucseq).
#'
#' @export
DefaultPreprocConfig <- function() {
  return (NucseqPreprocConfig())
}

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Validity
###

setValidity("PreprocConfig", function(object) {

  if (!is.numeric(object@complexityCutoff))
    return (sprintf("complexity cutoff must be of numeric type"))
  if (length(object@complexityCutoff) != 2)
    return (sprintf("complexity cutoff must be a numeric vector of length equals to 2"))
  if (object@complexityCutoff[1] < 0)
    return (sprintf("lower bound of complexity cutoff must be greater than or equal to 0"))
  if (object@complexityCutoff[2] <= object@complexityCutoff[1])
    return (sprintf("upper bound of complexity cutoff must be greater than lower bound"))

  if (object@expressionCutoff <= 0)
    return (sprintf("expression cutoff must be greater than or equal to 0"))

  if (object@housekeepingCutoff <= 0)
    return (sprintf("housekeeping cutoff must be greater than 0"))

  if (object@logBase <= 0)
    return (sprintf("log base must be greater than 0"))

  if (object@scalingFactor <= 0)
    return (sprintf("scaling factor must be greater than 0"))

  if (object@pseudoCount <= 0)
    return (sprintf("Pseudo count must be greater than 0"))

  return(TRUE)
})

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Getters and setters.
###

#'
#' @rdname PreprocConfig-class
#'
#' @export
setMethod("complexityCutoff", "PreprocConfig", function(x) {
  return(x@complexityCutoff)
})

#'
#' @rdname PreprocConfig-class
#'
#' @export
setReplaceMethod("complexityCutoff", "PreprocConfig", function(x, value) {
  x@complexityCutoff <- value
  return(x)
})

#'
#' @rdname PreprocConfig-class
#'
#' @export
setMethod("expressionCutoff", "PreprocConfig", function(x) {
  return(x@expressionCutoff)
})

#'
#' @rdname PreprocConfig-class
#'
#' @export
setReplaceMethod("expressionCutoff", "PreprocConfig", function(x, value) {
  x@expressionCutoff <- value
  return(x)
})

#'
#' @rdname PreprocConfig-class
#'
#' @export
setMethod("housekeepingCutoff", "PreprocConfig", function(x) {
  return(x@housekeepingCutoff)
})

#'
#' @rdname PreprocConfig-class
#'
#' @export
setReplaceMethod("housekeepingCutoff", "PreprocConfig", function(x, value) {
  x@housekeepingCutoff <- value
  return(x)
})

#'
#' @rdname PreprocConfig-class
#'
#' @export
setMethod("logBase", "PreprocConfig", function(x) {
  return(x@logBase)
})

#'
#' @rdname PreprocConfig-class
#'
#' @export
setReplaceMethod("logBase", "PreprocConfig", function(x, value) {
  x@logBase <- value
  return(x)
})

#'
#' @rdname PreprocConfig-class
#'
#' @export
setMethod("scalingFactor", "PreprocConfig", function(x) {
  return(x@scalingFactor)
})

#'
#' @rdname PreprocConfig-class
#'
#' @export
setReplaceMethod("scalingFactor", "PreprocConfig", function(x, value) {
  x@scalingFactor <- value
  return(x)
})

#'
#' @rdname PreprocConfig-class
#'
#' @export
setMethod("pseudoCount", "PreprocConfig", function(x) {
  return(x@pseudoCount)
})

#'
#' @rdname PreprocConfig-class
#'
#' @export
setReplaceMethod("pseudoCount", "PreprocConfig", function(x, value) {
  x@pseudoCount <- value
  return(x)
})

#'
#' @rdname PreprocConfig-class
#'
#' @export
setMethod("typeMatrix", "PreprocConfig", function(x) {
  return(x@typeMatrix)
})

#'
#' @rdname PreprocConfig-class
#'
#' @export
setReplaceMethod("typeMatrix", "PreprocConfig", function(x, value) {
  x@typeMatrix <- value
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
  cat("Log base:", logBase(object), "\n")
  cat("Scaling factor:", scalingFactor(object), "\n")
  cat("Pseudo count:", pseudoCount(object), "\n")
  cat("Matrix type:", typeMatrix(object), "\n")
})

is_config_object <- function(object) { return (!is.null(object) & is(object, "PreprocConfig")) }

### -------------------------------------------------------------------------
### PreprocConfig objects (end)
### =========================================================================
###

### =========================================================================
### ScandalDataSet objects (start)
### -------------------------------------------------------------------------
###

setClassUnion("ScandalDataSetOrNULL", c("NULL"))
setClassUnion("MatrixOrNULL", c("Matrix", "matrix", "NULL"))

#'
#' @title ScandalDataSet class
#'
#' @description An S4 class for storing single-cell seqeuncing data, analysis
#'
#' @details The S4 class \code{ScandalDataSet}
#'
#' @slot childNodes a \code{SimpleList} object containing the child nodes of the
#' constructed \code{ScandalDataSet} object.
#' @slot parentNode the parent node of the constructed \code{ScandalDataSet}
#' object.
#' @slot unprocessedData a read-only matrix that contains the unprocessed data that
#' allows re-accessing this data without the need to read it from file. Sparse matrix
#' representation as well as maintaining a single copy for the entire objects tree
#' decreases the memory overhead of this approach.
#' @slot preprocConfig a configuration object of class \code{PreprocConfig}.
#' @slot nodeID a unique character identifier of the constructed \code{ScandalDataSet}
#' object that should represent the specific sample.
#' @slot projectID a character identifier common to all the nodes in the constructed
#' \code{ScandalDataSet} object.
#'
#' @section Constructor:
#' Constructs a \code{ScandalDataSet} object.
#'
#' @section Methods:
#' \describe{
#'   \item{\code{logtpm}}{Getter/setter for the logtpm assay}
#'   \item{\code{childNodes}}{Getter/setter for childNodes}
#'   \item{\code{parentNode}}{Getter/setter for the parentNode}
#'   \item{\code{unprocessedData}}{Getter for the unprocessedData (read-only)}
#'   \item{\code{preprocConfig}}{Getter for the preprocConfig (read-only)}
#'   \item{\code{nodeID}}{Getter/setter for the nodeID}
#'   \item{\code{projectID}}{Getter/setter for the projectID}
#' }
#'
#' @examples
#'
#'
#'
#' @rdname ScandalDataSet
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

#  Will be added by Roxygen to the class documentation

#' @usage
#' ## Constructor
#' ScandalDataSet(..., childNodes = S4Vectors::SimpleList(),
#'   parentNode = NULL, preprocConfig = DefaultPreprocConfig(),
#'   nodeID = NODE_ID(), projectID = PROJ_ID())
#'
#' @param ... arguments to pass to the \code{SingleCellExperiment} constructor.
#' @param childNodes a \code{SimpleList} object containing the child nodes of the
#' constructed \code{ScandalDataSet} object. Default is an empty \code{SimpleList}
#' @param parentNode the parent node of the constructed \code{ScandalDataSet}
#' object. Default is \code{NULL}, meaning that the constructed object has no
#' parent.
#' @param preprocConfig a configuration object of class \code{PreprocConfig}
#' @param nodeID a unique identifier of the constructed \code{ScandalDataSet}
#' object that should represent the specific sample. If not supplied a unique ID
#' will be generated randomly however it is advised to set this field.
#' @param projectID an identifier common to all the nodes in the constructed
#' \code{ScandalDataSet} object.  If not supplied a unique ID
#' will be generated randomly however it is advised to set this field.
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

  message("In setValidity")

  if (length(childNodes(object)) > 0) {

    is_child_valid <- sapply(childNodes(object), function(c) is(c, "ScandalDataSet"))

    if (!(base::all(is_child_valid) == TRUE))
      return (sprintf("Every child node must be a ScandalDataSet object"))
  }

  return (TRUE)
})

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Getters and setters.
###

#'
#' @rdname ScandalDataSet
#'
#' @export
setMethod("logtpm", "ScandalDataSet",   function(object, ...) {
  return (assay(object, i = "logtpm", ...))
})

#'
#' @param value a value to replace the currently set value (applies to all setter methods).
#'
#' @rdname ScandalDataSet
#'
#' @export
setReplaceMethod("logtpm", c("ScandalDataSet", "ANY"),   function(object, ..., value) {
  assay(object, i = "logtpm", ...) <- value
  return (object)
})

#'
#' @rdname ScandalDataSet
#'
#' @export
setMethod("childNodes", "ScandalDataSet", function(object) {
  return(object@childNodes)
})

#'
#' @rdname ScandalDataSet
#'
#' @export
setReplaceMethod("childNodes", "ScandalDataSet", function(object, value) {
  object@childNodes <- value
  return(object)
})

#'
#' @rdname ScandalDataSet
#'
#' @export
setMethod("parentNode", "ScandalDataSet", function(object) {
  return(object@parentNode)
})

#'
#' @rdname ScandalDataSet
#'
#' @export
setReplaceMethod("parentNode", "ScandalDataSet", function(object, value) {
  object@parentNode <- value
  return(object)
})

#'
#' @rdname ScandalDataSet
#'
#' @export
setMethod("nodeID", "ScandalDataSet", function(object) {
  return(object@nodeID)
})

#'
#' @rdname ScandalDataSet
#'
#' @export
setMethod("projectID", "ScandalDataSet", function(object) {
  return(object@projectID)
})

#'
#' @rdname ScandalDataSet
#'
#' @export
setMethod("sampleNames", "ScandalDataSet", function(object) {
  return(unique(gsub("*-.*", "", colnames(object))))
})

#'
#' @rdname ScandalDataSet
#'
#' @export
setMethod("unprocessedData", "ScandalDataSet", function(object) {

  o <- object
  while(!is.null(parentNode(o)))
    o <- parentNode(o)

  sname <- base::unique(.cell2tumor(colnames(object)))

  # return only the cells that belong to this specific tumor
  return (o@unprocessedData[, .subset_cells(colnames(o@unprocessedData), sname), drop = FALSE])
})

#'
#' @rdname ScandalDataSet
#'
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
