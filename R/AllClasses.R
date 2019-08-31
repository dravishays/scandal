
#' @title S4 classes definition
#' @details Includes the generics definition file
#' @include AllGenerics.R
NULL # Do not remove me!!!

### =========================================================================
### PreprocConfig objects (start)
### -------------------------------------------------------------------------
###

setClassUnion("ListOrVector", c("list", "vector"))
setClassUnion("NumericOrVector", c("numeric", "vector"))

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
#' pc <- PreprocConfig(complexityCutoff = c(0, 10000),
#'                     expressionCutoff = 5,
#'                     housekeepingCutoff = 7,
#'                     logBase = 2,
#'                     scalingFactor = 10,
#'                     pseudoCount = 1,
#'                     typeMatrix = TRUE)
#'
#' logBase(pc) # Equals 2
#' logBase(pc) <- 10
#' logBase(pc) # Equals 10
#'
#' @aliases PreprocConfig
#'
#' @author Avishay Spitzer
#'
#' @export
setClass("PreprocConfig",
         slots = c(complexityCutoff = "ListOrVector",
                   expressionCutoff = "NumericOrVector",
                   housekeepingCutoff = "NumericOrVector",
                   logBase = "numeric",
                   scalingFactor = "numeric",
                   pseudoCount = "numeric",
                   typeMatrix = "logical"))

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Constructors
###

#'
#' @param complexityCutoff A numeric vector of length 2 representing the lower and
#' upper bounds of complexity (i.e. the number of detected genes per cell).
#' @param expressionCutoff A numeric representing the minimal log2 mean expression
#' per gene below which a gene is considered lowly expressed.
#' @param housekeepingCutoff A numeric representing the log2 mean expression of
#' house-keeping genes (i.e. genes that are highly expressed in all cells) per
#' cell below which a cell is considered low quality.
#' @param logBase A numeric representing the logarithm base for performing log
#' transformation on the data.
#' @param scalingFactor A numeric representing a scaling factor by which to divide
#' each data point before log transformation.
#' @param pseudoCount A numeric representing the pseudo count added when performing
#' log transformation to avoid taking the log of zero.
#' @param typeMatrix A logical indicating if the dataset should be represented using
#' the S4 Matrix class (instead of base R matrix) to reduce memory overhead using
#' sparse matrix representation.
#'
#' @describeIn PreprocConfig-class Constructs a new \code{PreprocConfig} object.
#'
#' @importFrom methods new
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

.valid_complexity_cutoff <- function(x) {
  if (!is.numeric(x))
    return (sprintf("complexity cutoff must be of numeric type"))
  if (length(x) != 2)
    return (sprintf("complexity cutoff must be a numeric vector of length equals to 2"))
  if (x[1] < 0)
    return (sprintf("lower bound of complexity cutoff must be greater than or equal to 0"))
  if (x[2] <= x[1])
    return (sprintf("upper bound of complexity cutoff must be greater than lower bound"))

  return (NULL)
}

setValidity("PreprocConfig", function(object) {

  if (is.list(object@complexityCutoff)) {
    for (co in object@complexityCutoff) {
      res <- .valid_complexity_cutoff(co)

      if (!is.null(res))
        return (res)
    }
  } else {
    res <- .valid_complexity_cutoff(object@complexityCutoff)

    if (!is.null(res))
      return (res)
  }

  if (length(object@expressionCutoff) > 1){
    if (is.null(names(object@expressionCutoff)))
      return (sprintf("expression cutoff vector must be a named vector"))
  }

  if (any(object@expressionCutoff <= 0))
    return (sprintf("expression cutoff must be greater than or equal to 0"))

  if (length(object@housekeepingCutoff) > 1){
    if (is.null(names(object@housekeepingCutoff)))
      return (sprintf("housekeeping cutoff vector must be a named vector"))
  }

  if (any(object@housekeepingCutoff <= 0))
    return (sprintf("housekeeping cutoff must be greater than or equal to 0"))

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
setMethod("complexityCutoff", "PreprocConfig", function(x, sid = NULL) {
  if (is.list(x@complexityCutoff) & !is.null(sid))
    return (x@complexityCutoff[[sid]])

  return(x@complexityCutoff)
})

#'
#' @param x a \code{PreprocConfig} object.
#' @param value a value to replace the currently set value.
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
setMethod("expressionCutoff", "PreprocConfig", function(x, sid = NULL) {
  if (length(x@expressionCutoff) > 1 & !is.null(sid))
    return (x@expressionCutoff[sid])

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
setMethod("housekeepingCutoff", "PreprocConfig", function(x, sid = NULL) {
  if (length(x@housekeepingCutoff) > 1 & !is.null(sid))
    return (x@housekeepingCutoff[sid])

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

.valid_index <- function(x, i) {

  if (is.character(i))
    return(all(i %in% names(x)))
  else if (is.numeric(i)) {
    return (all(!(is.null(x[i]) | is.na(x[i]))))
  }

  return (FALSE)
}

setMethod("[", "PreprocConfig", function(x, i, ...) {

  if (missing(i))
    return(x)

  if (is.list(x@complexityCutoff)) {

    stopifnot(.valid_index(x@complexityCutoff, i))

    if (length(i) > 1)
      x@complexityCutoff <- x@complexityCutoff[i]
    else
      x@complexityCutoff <- x@complexityCutoff[[i]]
  }

  if (length(x@expressionCutoff) > 1) {

    stopifnot(.valid_index(x@expressionCutoff, i))

    x@expressionCutoff <- x@expressionCutoff[i]
  }

  if (length(x@housekeepingCutoff) > 1) {

    stopifnot(.valid_index(x@housekeepingCutoff, i))

    x@housekeepingCutoff <- x@housekeepingCutoff[i]
  }

  return (x)
})

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Display
###

#'
#' @rdname PreprocConfig-class
#'
#' @export
setMethod("show", "PreprocConfig", function(object) {

  if (is.list(complexityCutoff(object))) {

    str <- sprintf("Complexity cutoff:")
    for (n in names(complexityCutoff(object))) {
      co <- complexityCutoff(object)[[n]]
      str <- paste0(str, sprintf(" %s=(%d, %d)", n, co[1], co[2]))
    }

    cat(paste0(str, "\n"))

  } else
    cat(sprintf("Complexity cutoff: (%d, %d)\n", complexityCutoff(object)[1], complexityCutoff(object)[2]))

  if (length(expressionCutoff(object)) == 1)
    cat("Expression cutoff:", expressionCutoff(object), "\n")
  else {

    str <- sprintf("Expression cutoff:")
    for (n in names(expressionCutoff(object))) {
      co <- expressionCutoff(object)[n]
      str <- paste0(str, sprintf(" %s=%.2f", n, co))
    }

    cat(paste0(str, "\n"))
  }

  if (length(housekeepingCutoff(object)) == 1)
    cat("Housekeeping cutoff:", housekeepingCutoff(object), "\n")
  else {

    str <- sprintf("Housekeeping cutoff:")
    for (n in names(housekeepingCutoff(object))) {
      co <- housekeepingCutoff(object)[n]
      str <- paste0(str, sprintf(" %s=%.2f", n, co))
    }

    cat(paste0(str, "\n"))
  }

  cat("Log base:", logBase(object), "\n")
  cat("Scaling factor:", scalingFactor(object), "\n")
  cat("Pseudo count:", pseudoCount(object), "\n")
  cat("Matrix type:", typeMatrix(object), "\n")
})

#' @importFrom methods is
is_config_object <- function(object) { return (!is.null(object) & is(object, "PreprocConfig")) }

### -------------------------------------------------------------------------
### PreprocConfig objects (end)
### =========================================================================
###

### =========================================================================
### QCResults objects (start)
### -------------------------------------------------------------------------
###

#'
#' @title QCResults class
#'
#' @description An S4 class for storing quality control results.
#'
#' @slot preprocConfig A object of class \linkS4class{PreprocConfig}.
#' @slot cellIDs A character vector of IDs of cells that passed QC.
#' @slot geneIDs A character vector of IDs of genes that passe QC.
#' @slot statsQC A \linkS4class{DataFrame} containing statistics gathered through
#' quality control process.
#'
#' @section Methods:
#' \describe{
#'   \item{\code{preprocConfig}}{Getter for the \linkS4class{PreprocConfig} object
#'   used when creating the quality control result.}
#'   \item{\code{cellIDs}}{Getter for the cell IDs that passed QC.}
#'   \item{\code{geneIDs}}{Getter for the gene IDs that passed QC}
#'   \item{\code{statsQC}}{Getter for the quality control statistics}
#' }
#'
#' @examples
#'
#' @aliases QCResults
#'
#' @author Avishay Spitzer
#'
#' @export
setClass("QCResults",
         slots = c(preprocConfig = "PreprocConfig",
                   cellIDs = "vector",
                   geneIDs = "vector",
                   statsQC = "DataFrame"))

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Constructors
###

#'
#' @describeIn QCResults-class Constructs a new \code{QCResults} object.
#'
#' @importFrom methods new
#' @importClassesFrom S4Vectors DataFrame
#' @importFrom S4Vectors DataFrame
#'
#' @export
QCResults <- function(preprocConfig, cellIDs = c(), geneIDs = c(), statsQC = DataFrame()) {
  qc <- new("QCResults",
            preprocConfig = preprocConfig,
            cellIDs = cellIDs,
            geneIDs = geneIDs,
            statsQC = statsQC)

  return (qc)
}

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Getters and setters.
###

#'
#' @rdname QCResults-class
#'
#' @export
setMethod("preprocConfig", "QCResults", function(object) {
  return(object@preprocConfig)
})

#'
#' @rdname QCResults-class
#'
#' @export
setMethod("cellIDs", "QCResults", function(object) {
  return(object@cellIDs)
})

#'
#' @rdname QCResults-class
#'
#' @export
setReplaceMethod("cellIDs", "QCResults", function(object, value) {
  object@cellIDs <- value
  return(object)
})

#'
#' @rdname QCResults-class
#'
#' @export
setMethod("geneIDs", "QCResults", function(object) {
  return(object@geneIDs)
})

#'
#' @rdname QCResults-class
#'
#' @export
setReplaceMethod("geneIDs", "QCResults", function(object, value) {
  object@geneIDs <- value
  return(object)
})

#'
#' @rdname QCResults-class
#'
#' @export
setMethod("statsQC", "QCResults", function(object) {
  return(object@statsQC)
})

#'
#' @rdname QCResults-class
#'
#' @export
setReplaceMethod("statsQC", "QCResults", function(object, value) {
  object@statsQC <- value
  return(object)
})

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Display
###

#'
#' @rdname QCResults-class
#'
#' @export
setMethod("show", "QCResults", function(object) {
  cat("cellIDs:", length(cellIDs(object)), "IDs\n")
  cat("geneIDs:", length(geneIDs(object)), "IDs\n")
  show(statsQC(object))
})

#' @importFrom methods is
is_qc_results_object <- function(object) { return (!is.null(object) & is(object, "QCResults")) }

### -------------------------------------------------------------------------
### QCResults objects (end)
### =========================================================================
###

### =========================================================================
### ScandalDataSet objects (start)
### -------------------------------------------------------------------------
###

#' @import Matrix
setClassUnion("MatrixOrNULL", c("Matrix", "matrix", "NULL"))

#'
#' @title ScandalDataSet class
#'
#' @description An S4 class for storing single-cell seqeuncing data, reduced
#' dimensions representations of the data reuqired in the analysis process such as
#' t-SNE and UMAP coordinates and the end-product of the analysis which are the
#' transcriptional programs.
#'
#' @details The S4 class \code{ScandalDataSet} inherits from and extends Bioconductor's
#' base class for single-cell related applications, the \linkS4class{SingleCellExperiment}
#' class.
#' The idea behind \link{scandal} is that in order to detect intra-tumor heterogeneity
#' one needs inspect each tumor individually to collect the different transcriptomic
#' programs that can be found within each tumor and then assess these programs at the
#' level of the entire dataset to define the programs that generalize best
#' (meta-programs).
#' Besides the functionality supplied by its superclasses, \code{ScandalDataSet}
#' supplies methods to keep
#'
#' @slot unprocessedData A read-only matrix that contains the unprocessed data that
#' allows re-accessing this data without the need to read it from file. Sparse matrix
#' representation as well as maintaining a single copy for the entire objects tree
#' decreases the memory overhead of this approach.
#' @slot preprocConfig A configuration object of class \linkS4class{PreprocConfig}.
#' @slot qualityControl A \linkS4class{SimpleList} object containing objects of class
#' \code{QCResults} representing the quality control results, i.e. the cells and
#' genes in each individual node that passed qaulity control and are available for
#' downstream analysis.
#' @slot nodeID A unique character identifier of the constructed \code{ScandalDataSet}
#' object that should represent the specific sample.
#' @slot projectID A character identifier common to all the nodes in the constructed
#' \code{ScandalDataSet} object.
#' @slot cell2SampleMap A **function** that maps a vector of cell IDs to a vector of
#' sample IDs to which the cells belong.
#'
#' @section Constructor:
#' Constructs a \code{ScandalDataSet} object.
#'
#' @section Methods:
#' \describe{
#'   \item{\code{logtpm}}{Getter/setter for the logtpm assay}
#'   \item{\code{unprocessedData}}{Getter for the unprocessedData (read-only)}
#'   \item{\code{preprocConfig}}{Getter for the preprocConfig (read-only)}
#'   \item{\code{qualityControl}}{Getter/Setter for the qualityControl list}
#'   \item{\code{nodeID}}{Getter/setter for the nodeID}
#'   \item{\code{projectID}}{Getter/setter for the projectID}
#'   \item{\code{sampleIDs}}{Returns a character vector containing the IDs of all
#'   samples in the dataset.}
#'   \item{\code{inspectSamples}}{Returns a ScandalDataSet object representing a
#'   set of specific samples.}
#'   \item{\code{cell2SampleMap}}{Getter for the cell2SampleMap function (read-only)}
#' }
#'
#' @seealso \linkS4class{SummarizedExperiment}, \linkS4class{SingleCellExperiment}, \link{scandal_preprocess}
#'
#' @examples
#' # Building a mock dataset with 30 cells and 100 genes
#' ngenes <- 100
#' ncells <- 30
#' dataset <- matrix(sample(0:1e4, ngenes * ncells, replace = FALSE), nrow = ngenes, ncol = ncells)
#' rownames(dataset) <- sapply(seq_len(ngenes), function(x) paste0("GENE", x))
#' colnames(dataset) <- c(sapply(seq_len(ncells / 2), function(x) paste0("TUMOR1-Cell", x)),
#'                        sapply(seq(from = ncells / 2 + 1, to = ncells), function(x) paste0("TUMOR2-Cell", (x - ncells/2))))
#'
#' # Declare a global confguration object for the top-level ScandalDataSet object and
#' # a named list of configuration objects for each single tumor. Note that the names
#' # of the elements in the named list correspond to the names of the tumors that appear
#' # in the column names
#' global_config <- PreprocConfig(complexityCutoff = c(0, 10000), expressionCutoff = 1, housekeepingCutoff = 1, logBase = 2, scalingFactor = 1, pseudoCount = 1, typeMatrix = TRUE)
#' tumor_config <- list(TUMOR1 = PreprocConfig(complexityCutoff = c(0, 10000), expressionCutoff = 1, housekeepingCutoff = 1, logBase = 2, scalingFactor = 1, pseudoCount = 1, typeMatrix = TRUE),
#'                      TUMOR2 = PreprocConfig(complexityCutoff = c(0, 10000), expressionCutoff = 1, housekeepingCutoff = 1, logBase = 2, scalingFactor = 1, pseudoCount = 1, typeMatrix = TRUE))
#'
#' # Instantiate a new ScandalDataSet object
#' sds <- ScandalDataSet(assays = list(tpm = dataset), preprocConfig = global_config, nodeID = "Example1", projectID = "Project1")
#'
#' sds # Prints a user-readable summary of sds
#'
#' all(colnames(sds) == colnames(dataset)) # TRUE
#' all(rownames(sds) == rownames(dataset)) # TRUE
#' sampleIDs(sds) # Return a vector (TUMOR1, TUMOR2)
#' qualityControl(sds) # Empty list
#' nodeID(sds) # Returns "Example1"
#' projectID(sds) # Returns "Project1"
#'
#' @rdname ScandalDataSet
#'
#' @author Avishay Spitzer
#'
#' @import SingleCellExperiment
#' @import SummarizedExperiment
#' @importClassesFrom S4Vectors SimpleList
#'
#' @export
#' @exportClass ScandalDataSet
setClass("ScandalDataSet",
         slots = c(unprocessedData = "MatrixOrNULL",
                   preprocConfig = "PreprocConfig",
                   qualityControl = "SimpleList",
                   nodeID = "character",
                   projectID = "character",
                   cell2SampleMap = "function"),
         contains = "SingleCellExperiment"
)

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Constructor
###

#  Will be added by Roxygen to the class documentation

#' @usage
#' ## Constructor
#' ScandalDataSet(..., preprocConfig = DefaultPreprocConfig(),
#'   nodeID = NODE_ID(), projectID = PROJ_ID())
#'
#' @param ... arguments to pass to the \linkS4class{SingleCellExperiment} constructor.
#' @param preprocConfig a configuration object of class \linkS4class{PreprocConfig}
#' @param nodeID a unique identifier of the constructed \code{ScandalDataSet}
#' object that should represent the specific sample. If not supplied a unique ID
#' will be generated randomly however it is advised to set this field.
#' @param projectID an identifier common to all the nodes in the constructed
#' \code{ScandalDataSet} object.  If not supplied a unique ID
#' will be generated randomly however it is advised to set this field.
#' @param cell2SampleMap a **function** that maps a vector of cell IDs to a vector of
#' sample IDs to which the cells belong. The default function assumes that the cell ID
#' is a string separated by "-" and that the node ID is contained in the substring
#' until the first "-" character.
#'
#' @importClassesFrom S4Vectors DataFrame SimpleList
#' @importFrom S4Vectors SimpleList
#' @importFrom S4Vectors DataFrame
#' @importFrom methods new is as
#'
#' @export
ScandalDataSet <- function(..., preprocConfig = DefaultPreprocConfig(), nodeID = NODE_ID(), projectID = PROJ_ID(), cell2SampleMap = DEFAULT_CELL_2_NODE_MAP) {

  sce <- SingleCellExperiment(...)

  if(!is(sce, "SingleCellExperiment")) {
    sce <- as(sce, "SingleCellExperiment")
  }

  object <- new("ScandalDataSet", sce,
                unprocessedData = assay(sce),
                preprocConfig = preprocConfig,
                qualityControl = SimpleList(),
                nodeID = nodeID,
                projectID = projectID,
                cell2SampleMap = cell2SampleMap)

  int_colData(object)$Scandal <- DataFrame(row.names = colnames(object))
  int_elementMetadata(object)$Scandal <- DataFrame(row.names = rownames(object))
  int_metadata(object)$Scandal <- list()

  int_metadata(object)$Scandal[["Version"]] <- 1.0

  return (object)
}

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Validity
###

setValidity("ScandalDataSet", function(object) {

  if (length(qualityControl(object)) > 0) {

    is_qc_valid <- sapply(qualityControl(object), function(c) is(c, "QCResults"))

    if (!(base::all(is_qc_valid) == TRUE))
      return (sprintf("Every QC node must be a QCResults object"))
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
#' @param object a \code{ScandalDataSet} object.
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
setMethod("qualityControl", "ScandalDataSet", function(object) {
  return(object@qualityControl)
})

#'
#' @rdname ScandalDataSet
#'
#' @export
setReplaceMethod("qualityControl", "ScandalDataSet", function(object, value) {
  object@qualityControl <- value
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
setReplaceMethod("nodeID", "ScandalDataSet", function(object, value) {
  object@nodeID <- value
  return(object)
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
setMethod("sampleIDs", "ScandalDataSet", function(object) {
  return(unique(cell2SampleMap(object)(colnames(object))))
})

#'
#' @rdname ScandalDataSet
#'
#' @export
setMethod("unprocessedData", "ScandalDataSet", function(object) {
  return (object@unprocessedData)
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

#'
#' @rdname ScandalDataSet
#'
#' @export
setMethod("cell2SampleMap", "ScandalDataSet", function(object) {
  return (object@cell2SampleMap)
})

#'
#' @rdname ScandalDataSet
#'
#' @importClassesFrom S4Vectors DataFrame SimpleList
#' @importFrom S4Vectors SimpleList
#' @importFrom S4Vectors DataFrame
#'
#' @export
setMethod("inspectSamples", "ScandalDataSet", function(object, sampleIDs, nodeID = NODE_ID()) {

  stopifnot(!is.null(sampleIDs), is.character(sampleIDs), base::all(sampleIDs %in% sampleIDs(object)) == TRUE)

  qc_list <- qualityControl(object)[sampleIDs]

  stopifnot(!is.null(qc_list), base::all(sampleIDs %in% names(qc_list)) == TRUE)

  res <- object[unique(unname(unlist(sapply(qc_list, function(qc) geneIDs(qc))))),
                unname(unlist(sapply(qc_list, function(qc) cellIDs(qc))))]

  reducedDims(res) <- SimpleList()
  res@nodeID <- nodeID
  res@preprocConfig <- preprocConfig(object)[sampleIDs]
  res@qualityControl <- qc_list
  names(res@qualityControl) <- sampleIDs
  res@unprocessedData <- res@unprocessedData[, .subset_cells(colnames(res@unprocessedData), sampleIDs, cell2SampleMap(object))]

  return(res)
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

#'
#' @rdname ScandalDataSet
#'
#' @importFrom methods callNextMethod
#'
#' @export
setMethod("show", "ScandalDataSet", function(object) {
  callNextMethod()
  scat("sampleIDs(%d): %s\n", sampleIDs(object))
  cat("nodeID:", nodeID(object), "\n")
  cat("projectID:", projectID(object), "\n")
  scat("qualityControl(%d): %s\n", names(qualityControl(object)))
  cat("unprocessedData:", class(unprocessedData(object)), "with", NROW(unprocessedData(object)), "rows and", NCOL(unprocessedData(object)), "columns\n")
  show(preprocConfig(object))
})

#' @importFrom methods is
is_scandal_object <- function(object) { return (!is.null(object) & is(object, "ScandalDataSet")) }

#' @importFrom methods is
is_valid_assay <- function(x) { return (!(is.null(x)) & (is(x, "Matrix") | is.matrix(x))) }

# Generates a random experiment ID by sampling a random integer
NODE_ID <- function() { paste0("NODE", base::sample(1:1e9, 1, replace = FALSE)) }
PROJ_ID <- function() { paste0("PROJ", base::sample(1:1e9, 1, replace = FALSE)) }

DEFAULT_CELL_2_NODE_MAP <- function(cell_ids) {

  stopifnot(!is.null(cell_ids), is.vector(cell_ids), is.character(cell_ids))

  return (base::gsub("-.*", "", cell_ids))
}

### -------------------------------------------------------------------------
### ScandalDataSet objects (end)
### =========================================================================
###

### =========================================================================
### ScandalResults objects (start)
### -------------------------------------------------------------------------
###

#' @export
#' @exportClass ScandalResults
setClass("ScandalResults",
         slots = c(samples = "SimpleList", # List of ScandalDataSet objects
                   wsClusteringData = "SimpleList", # List of within-sample clustering data (e.g. NMF objects)
                   wsPrograms = "SimpleList", # List of lists representing within-sample programs (vectors of gene symbols)
                   wsScores = "SimpleList", # List of matrices representing the program scores within each sample
                   wsScoreSDs = "SimpleList", # List of vectors representing the standard deviation of program scores within each sample
                   bsScores = "matrix",
                   variablePrograms = "SimpleList",
                   programClusters = "SimpleList",
                   geneScores = "matrix",
                   mpScores = "SimpleList",
                   thresholdSD = "numeric", # Minimal SD below which programs are filtered out
                   nodeID = "character",
                   projectID = "character"),
         contains = "DataFrame"
)

### -------------------------------------------------------------------------
### ScandalResults objects (end)
### =========================================================================
###
