
#'
#' @title
#'
#' @description
#'
#' @param object
#'
#' @details
#'
#' @return
#'
#' @author Avishay Spitzer
#'
#' @export
preprocess <- function(object, config_param_list, forced_genes_set = NULL, use_housekeeping_filter = FALSE, aggregate_res = TRUE) {

  stopifnot(!is_scandal_object(object), is.logical(use_housekeeping_filter), is.logical(aggregate_res))

  if (is.null(config_param_list)) {
    warning("received empty config_param_list, will use default configuration")

    config_param_list <- .default_config(object)
  } else
    stopifnot(is.list(config_param_list),
              base::setequal(sampleNames(object), names(config_param_list)),
              base::all(base::lapply(config_param_list, function(x) is_config_object(x)) == TRUE))

  for(sname in sampleNames(object)) {
    sconf <- config_param_list[[sname]]
    sdata <- ScandalDataSet(assays = list(tpm = tpm(object)[, subset_cells(colnames(object), sname), drop = FALSE]), identifier = sname, configParam = sconf)

    if (is.null(sconf))
      stop(paste0("No configuration object supplied for ", sname))

    sdata <- .preprocess(sdata, forced_genes_set = forced_genes_set, use_housekeeping_filter = use_housekeeping_filter)

    # Set plate data
    plate <- gsub("-.*", "", substr(colnames(sdata), nchar(sname) + 2, nchar(colnames(sdata))))
    names(plate) <- colnames(sdata)
    colData(sdata)$Plate <- plate

    indSampleData(object)[[sname]] <- sdata
  }

  colData(object)$Tumor <- setNames(.cell2tumor(colnames(object)), colnames(object))
  colData(object)$CellType <- rep("Undetermined", ncol(object))
  names(colData(object)$CellType) <- colnames(object)

  if (isTRUE(aggregate_res)) {
    plot_mean_complexity(object, save_to_file = TRUE, filename = "complexity_per_tumor_pre_qc.png", title = "Complexity per tumor (pre QC)")
    object <- aggregate_processed_data(object, forced_genes_set)
    plot_mean_complexity(object, save_to_file = TRUE, filename = "complexity_per_tumor_post_qc.png", title = "Complexity per tumor (post QC)")
  }

  return (object)
}

.cell2tumor <- function(cell_ids) gsub("-.*", "", cell_ids)

.default_config <- function(object) {
  res <- rep(list(DefaultConfigParam()), length(sampleNames(object)))
  names(res) <- sampleNames(object)
  return (res)
}
