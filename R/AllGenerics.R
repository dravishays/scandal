
### =========================================================================
### Genercis for class ScandalDataSet
### -------------------------------------------------------------------------
###

#'
#' @title ScandalDataSet generics
#'
#' @rdname ScandalDataSet-generics
#'
#' @export
setGeneric("logtpm", function(object, ...) standardGeneric("logtpm"))

#' @rdname ScandalDataSet-generics
#'
#' @export
setGeneric("logtpm<-", function(object, ..., value) standardGeneric("logtpm<-"))

#' @rdname ScandalDataSet-generics
#'
#' @export
setGeneric("qualityControl", function(object, ...) standardGeneric("qualityControl"))

#' @rdname ScandalDataSet-generics
#'
#' @export
setGeneric("qualityControl<-", function(object, ..., value) standardGeneric("qualityControl<-"))

#' @rdname ScandalDataSet-generics
#'
#' @export
setGeneric("nodeID", function(object, ...) standardGeneric("nodeID"))

#' @rdname ScandalDataSet-generics
#'
#' @export
setGeneric("nodeID<-", function(object, ..., value) standardGeneric("nodeID<-"))

#' @rdname ScandalDataSet-generics
#'
#' @export
setGeneric("projectID", function(object, ...) standardGeneric("projectID"))

#' @rdname ScandalDataSet-generics
#'
#' @export
setGeneric("sampleIDs", function(object, ..., return_sorted = FALSE) standardGeneric("sampleIDs"))

#' @rdname ScandalDataSet-generics
#'
#' @export
setGeneric("unprocessedData", function(object, ...) standardGeneric("unprocessedData"))

#' @rdname ScandalDataSet-generics
#'
#' @export
setGeneric("preprocConfig", function(object, ...) standardGeneric("preprocConfig"))

#' @rdname ScandalDataSet-generics
#'
#' @export
setGeneric("inspectSamples", function(object, nodeID, ...) standardGeneric("inspectSamples"))

#' @rdname ScandalDataSet-generics
#'
#' @export
setGeneric("cell2SampleMap", function(object, ...) standardGeneric("cell2SampleMap"))

### =========================================================================
### Genercis for class ConfigParam
### -------------------------------------------------------------------------
###

#'
#' @title PreprocConfig generics
#'
#' @rdname PreprocConfig-generics
#'
#' @export
setGeneric("complexityCutoff", function(x, sid = NULL, ...) standardGeneric("complexityCutoff"))

#' @rdname PreprocConfig-generics
#'
#' @export
setGeneric("complexityCutoff<-", function(x, ..., value) standardGeneric("complexityCutoff<-"))

#' @rdname PreprocConfig-generics
#'
#' @export
setGeneric("expressionCutoff", function(x, sid = NULL, ...) standardGeneric("expressionCutoff"))

#' @rdname PreprocConfig-generics
#'
#' @export
setGeneric("expressionCutoff<-", function(x, ..., value) standardGeneric("expressionCutoff<-"))

#' @rdname PreprocConfig-generics
#'
#' @export
setGeneric("housekeepingCutoff", function(x, sid = NULL, ...) standardGeneric("housekeepingCutoff"))

#' @rdname PreprocConfig-generics
#'
#' @export
setGeneric("housekeepingCutoff<-", function(x, ..., value) standardGeneric("housekeepingCutoff<-"))

#' @rdname PreprocConfig-generics
#'
#' @export
setGeneric("scalingFactor", function(x, ...) standardGeneric("scalingFactor"))

#' @rdname PreprocConfig-generics
#'
#' @export
setGeneric("scalingFactor<-", function(x, ..., value) standardGeneric("scalingFactor<-"))

#' @rdname PreprocConfig-generics
#'
#' @export
setGeneric("logBase", function(x, ...) standardGeneric("logBase"))

#' @rdname PreprocConfig-generics
#'
#' @export
setGeneric("logBase<-", function(x, ..., value) standardGeneric("logBase<-"))

#' @rdname PreprocConfig-generics
#'
#' @export
setGeneric("pseudoCount", function(x, ...) standardGeneric("pseudoCount"))

#' @rdname PreprocConfig-generics
#'
#' @export
setGeneric("pseudoCount<-", function(x, ..., value) standardGeneric("pseudoCount<-"))

#' @rdname PreprocConfig-generics
#'
#' @export
setGeneric("typeMatrix", function(x, ...) standardGeneric("typeMatrix"))

#' @rdname PreprocConfig-generics
#'
#' @export
setGeneric("typeMatrix<-", function(x, ..., value) standardGeneric("typeMatrix<-"))

### =========================================================================
### Genercis for class QCResults
### -------------------------------------------------------------------------
###

#'
#' @title QCResults generics
#'
#' @rdname QCResults-generics
#'
#' @export
setGeneric("cellIDs", function(object, ...) standardGeneric("cellIDs"))

#' @rdname QCResults-generics
#'
#' @export
setGeneric("cellIDs<-", function(object, ..., value) standardGeneric("cellIDs<-"))

#' @rdname QCResults-generics
#'
#' @export
setGeneric("geneIDs", function(object, ...) standardGeneric("geneIDs"))

#' @rdname QCResults-generics
#'
#' @export
setGeneric("geneIDs<-", function(object, ..., value) standardGeneric("geneIDs<-"))

#' @rdname QCResults-generics
#'
#' @export
setGeneric("statsQC", function(object, ...) standardGeneric("statsQC"))

#' @rdname QCResults-generics
#'
#' @export
setGeneric("statsQC<-", function(object, ..., value) standardGeneric("statsQC<-"))

### =========================================================================
### Genercis for class ScandalMetaprograms
### -------------------------------------------------------------------------
###

#'
#' @title ScandalMetaprograms generics
#'
#' @rdname ScandalMetaprograms-generics
#'
#' @export
setGeneric("l2R", function(object, ..., as_tibble = FALSE) standardGeneric("l2R"))

#'
#' @title ScandalMetaprograms generics
#'
#' @rdname ScandalMetaprograms-generics
#'
#' @export
setGeneric("corL2R", function(object, ...) standardGeneric("corL2R"))

#'
#' @title ScandalMetaprograms generics
#'
#' @rdname ScandalMetaprograms-generics
#'
#' @export
setGeneric("consensusClusters", function(object, ..., best = FALSE) standardGeneric("consensusClusters"))

#'
#' @title ScandalMetaprograms generics
#'
#' @rdname ScandalMetaprograms-generics
#'
#' @export
setGeneric("mpL2R", function(object, ..., as_tibble = FALSE) standardGeneric("mpL2R"))

#'
#' @title ScandalMetaprograms generics
#'
#' @rdname ScandalMetaprograms-generics
#'
#' @export
setGeneric("metaPrograms", function(object, ..., as_tibble = FALSE) standardGeneric("metaPrograms"))

#'
#' @title ScandalMetaprograms generics
#'
#' @rdname ScandalMetaprograms-generics
#'
#' @export
setGeneric("mpScores", function(object, ..., as_tibble = FALSE) standardGeneric("mpScores"))

#'
#' @title ScandalMetaprograms generics
#'
#' @rdname ScandalMetaprograms-generics
#'
#' @export
setGeneric("mpAssigned", function(object, ..., return_mapped = TRUE) standardGeneric("mpAssigned"))

#'
#' @title ScandalMetaprograms generics
#'
#' @rdname ScandalMetaprograms-generics
#'
#' @export
setGeneric("mpMap", function(object, ...) standardGeneric("mpMap"))

#' @rdname ScandalMetaprograms-generics
#'
#' @export
setGeneric("mpMap<-", function(object, ..., value) standardGeneric("mpMap<-"))

#'
#' @title ScandalMetaprograms generics
#'
#' @rdname ScandalMetaprograms-generics
#'
#' @export
setGeneric("scoringStrategy", function(object, ...) standardGeneric("scoringStrategy"))

#'
#' @title ScandalMetaprograms generics
#'
#' @rdname ScandalMetaprograms-generics
#'
#' @export
setGeneric("scoreThreshold", function(object, ...) standardGeneric("scoreThreshold"))
