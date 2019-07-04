
#'
#' @title
#'
#' @description
#'
#' @param object
#' @param preproc_config_list
#' @param forced_genes_set
#' @param use_housekeeping_filter
#'
#' @details
#'
#' @return
#'
#' @author Avishay Spitzer
#'
#' @export
scandal_preprocess <- function(object, preproc_config_list, forced_genes_set = NULL, use_housekeeping_filter = FALSE, verbose = FALSE) {

  stopifnot(!is_scandal_object(object),
            (is.null(forced_genes_set) | is.list(forced_genes_set)),
            is.logical(use_housekeeping_filter),
            is.logical(aggregate_res))

  stopifnot(!is.null(preproc_config_list),
            is.list(preproc_config_list),
            base::setequal(sampleNames(object), names(preproc_config_list)),
            base::all(base::lapply(preproc_config_list, function(x) is_config_object(x)) == TRUE))

  for(sname in sampleNames(object)) {
    sconf <- preproc_config_list[[sname]]
    sdata <- ScandalDataSet(assays = list(tpm = tpm(object)[, subset_cells(colnames(object), sname), drop = FALSE]), identifier = sname, preprocConfig = sconf)

    parentNode(sdata) <- object

    childNodes(object)[[sname]] <- .scandal_preprocess(sdata, cell_ids = NULL, forced_genes_set = forced_genes_set, use_housekeeping_filter = use_housekeeping_filter, verbose = verbose)
  }

  object <- .scandal_preprocess(object, cell_ids = .aggregate_cell_ids(object), forced_genes_set = forced_genes_set, use_housekeeping_filter = use_housekeeping_filter, verbose = verbose)

  return (object)
}

#' @export
preprocess <- function(x, complexity_cutoff, housekeeping_cutoff, expression_cutoff, log_base, scaling_factor,
                       sample_id = "", cell_ids = NULL, forced_genes_set = NULL, use_housekeeping_filter = FALSE, verbose = FALSE) {

  stopifnot(!is_valid_assay(x),
            (is.vecrot(complexity_cutoff) & is.numeric(complexity_cutoff)),
            is.numeric(housekeeping_cutoff),
            is.numeric(expression_cutoff),
            is.numeric(log_base),
            is.numeric(scaling_factor),
            (is.null(forced_genes_set) | is.vector(forced_genes_set)),
            (is.null(cell_ids) | is.vector(cell_ids)),
            is.logical(use_housekeeping_filter),
            is.logical(aggregate_res))

  if (isTRUE(verbose))
    message("Preprocessing sample ", sample_id, "...")

  sparcity_before_qc <- length(which(x == 0)) / (dim(x)[1] * dim(x)[2]) * 100

  if (isTRUE(verbose))
    message("Detecting high quality cells...")

  if (is.null(cell_ids)) {
    x <- x[, filter_low_quality_cells(x, complexity_cutoff = complexity_cutoff, verbose = verbose)]

    if (isFALSE(use_housekeeping_filter)) {
      if (isTRUE(verbose))
        message("Skipped filtering-out cells according to housekeeping cutoff")
    }
    else {
      if (isTRUE(verbose))
        message("Detecting low quality cells based on mean expression of housekeeping genes...")

      x <- x[, filter_low_housekeeping_cells(x, housekeeping_cutoff = housekeeping_cutoff, verbose = verbose)]
    }
  } else {

    if (isTRUE(verbose))
      message("Skipping cell filtering as cell IDs list was supplied")

    if (!base::all(cell_ids %in% colnames(x)))
      warning("Some (or all) of the supplied cell IDs were not found in the dataset")

    x <- x[, colnames(x) %in% cell_ids]
  }

  if (isTRUE(verbose))
    message("Detecting highly expressed genes...")

  x <- x[filter_lowly_expressed_genes(x, expression_cutoff = expression_cutoff, forced_genes_set = forced_genes_set, verbose = verbose), ]

  if (isTRUE(verbose))
    message(paste0("Log transforming TPM matrix, (base - ", log_base, ", scaling factor - ", scaling_factor, ")"))

  x <- log_transform(x, log_base = log_base, scaling_factor = scaling_factor, verbose = verbose)

  if (isTRUE(verbose))
    message("Centering matrix...")

  x <- center_matrix(x, by = "row", method = "mean", verbose = verbose)

  sparcity_after_qc <- length(which(x == 0)) / (dim(x)[1] * dim(x)[2]) * 100

  if (isTRUE(verbose))
    message(sprintf("Sparcity before QC: %.2f%%, after QC: %.2f%%", sparcity_before_qc, sparcity_after_qc))

  if (isTRUE(verbose))
    message("Preprocessing done!")

  return (x)
}

#' @export
log_transform <- function(x, log_base = 2, scaling_factor = 1, pseudo_count = 1, verbose = FALSE) {

  stopifnot(!is.null(x), log_base > 0, scaling_factor > 0, pseudo_count > 0)

  x <- base::log( (x / scaling_factor) + pseudo_count, base = log_base)

  return(x)
}

#' @export
reverse_log_transform <- function(x, log_base = 2, scaling_factor = 1, pseudo_count = 1, verbose = FALSE) {

  stopifnot(!is.null(x), log_base > 0, scaling_factor > 0, pseudo_count > 0)

  x <- (log_base^x * scaling_factor) - pseudo_count

  return(x)
}

#' @export
compute_complexity <- function(x, return_sorted = FALSE, cell_subset = NULL, verbose = FALSE) {

  stopifnot(!is_valid_assay(x), is.logical(return_sorted), (is.null(cell_subset) | is.vector(cell_subset)))

  if (!is.null(cell_subset))
    x <- x[, cell_subset]

  c <- base::apply(x, 2, function(y) base::sum(y != 0))

  if(isTRUE(return_sorted)) {
    c <- base::sort(c)
  }

  return(c)
}

#' @export
compute_central_tendency <- function(x, by = "row", method = "mean", log_transform_res = FALSE, genes_subset = NULL, verbose = FALSE) {

  stopifnot(!is_valid_assay(x), by %in% c("row", "col"), is.logical(log_transform_res), (is.null(genes_subset) | is.vector(genes_subset)))

  if (!is.null(genes_subset))
    x <- x[rownames(x) %in% genes_subset, ]

  func <- base::match.fun(method)

  x <- base::apply(x, ifelse(by == "row", 1, 2), func)

  if (isTRUE(log_transform_res))
    x <- base::log2(x + 1)

  return (x)
}

#' @export
filter_low_quality_cells <- function(x, complexity_cutoff, verbose = FALSE) {

  stopifnot(!is_valid_assay(x), is.numeric(complexity_cutoff), complexity_cutoff > 0)

  d <- compute_complexity(x, return_sorted = FALSE, cell_subset = NULL, verbose = verbose)

  passQC <- d[(d >= complexity_cutoff[1]) & (d <= complexity_cutoff[2]), drop = FALSE]

  if (isTRUE(verbose))
    message(sprintf("%d cells pre-QC, cell cutoff - lower bound %d [cells], upper bound %d [cells], %d cells passed QC, %2.0f%% of cells dropped",
                    length(d),
                    complexity_cutoff[1],
                    complexity_cutoff[2],
                    length(passQC),
                    (1 - length(passQC)/length(d)) * 100))

  return (names(passQC))
}

#' @export
filter_low_housekeeping_cells <- function(x, housekeeping_cutoff, verbose = FALSE) {

  stopifnot(!is_valid_assay(x), is.numeric(housekeeping_cutoff), housekeeping_cutoff > 0)

  hk_mean_exp <- compute_central_tendency(x, by = "col", method = "mean", log_transform_res = TRUE, genes_subset = HOUSEKEEPING_GENES_LIST, verbose = verbose)

  passQC <- hk_mean_exp[hk_mean_exp >= housekeeping_cutoff]

  if (isTRUE(verbose))
    message(sprintf("Housekeeping cutoff %d, %d cells pre-QC, %d cells passed QC, %2.0f%% of cells dropped",
                    housekeeping_cutoff,
                    length(hk_mean_exp),
                    length(passQC),
                    (1 - length(passQC)/length(hk_mean_exp)) * 100))

  return (names(passQC))
}

#' @export
filter_lowly_expressed_genes <- function(x, expression_cutoff, forced_genes_set = NULL, verbose = FALSE) {

  stopifnot(!is_valid_assay(x), is.numeric(expression_cutoff), expression_cutoff > 0, (is.null(forced_genes_set) | is.vector(forced_genes_set)))

  gene_counts <- compute_central_tendency(x, by = "row", method = "mean", log_transform_res = TRUE, genes_subset = forced_genes_set, verbose = verbose)

  passQC <- gene_counts[gene_counts >= expression_cutoff]

  if (!is.null(forced_genes_set))
    passQC_forced_genes <- gene_counts[names(which(gene_counts[forced_genes_set] > 0))]
  else
    passQC_forced_genes <- 0

  if (isTRUE(verbose))
    message(sprintf("%d genes pre-QC, gene cutoff %.2f, %d genes passed QC, %2.0f%% of genes dropped, mean expression of %d of %d forced genes is > 0",
                    length(gene_counts),
                    expression_cutoff,
                    length(passQC),
                    (1 - length(passQC)/length(gene_counts)) * 100,
                    ifelse(!is.null(forced_genes_set), length(passQC_forced_genes), 0),
                    ifelse(!is.null(forced_genes_set), length(forced_genes_set), 0)))

  return (unique(c(names(passQC), names(passQC_forced_genes))))
}

#'
#' @title Centers a matrix
#'
#' @description Centers each row or column of the given matrix around the mean or median
#'
#' @param x a numeric matrix
#' @param by
#' @param method the
#'
#' @details
#'
#' @return
#'
#' @author Avishay Spitzer
#'
#' @export
center_matrix <- function(x, by = "row", method = "mean", verbose = FALSE) {

  stopifnot(!is_valid_assay(x), by %in% c("row", "col"), method %in% c("mean", "median"))

  if (by == "row")
    x <- t(t(x) - compute_central_tendency(x, by = "row", method = method, log_transform_res = FALSE, genes_subset = NULL, verbose = verbose))
  else
    x <- x - compute_central_tendency(x, by = by, method = method, log_transform_res = FALSE, genes_subset = NULL, verbose = verbose)

  return (x)
}

subset_cells <- function(cell_names, sample_name) cell_names[which(.cell2tumor(cell_names) %in% sample_name)]

.cell2tumor <- function(cell_ids) gsub("-.*", "", cell_ids)

.scandal_preprocess <- function(object, cell_ids = NULL, forced_genes_set = NULL, use_housekeeping_filter = FALSE, verbose = FALSE) {

  # Extract the preprocessing configuration object
  preproc_config <- preprocConfig(object)

  # Call the matrix preprocessing function
  x <- preprocess(assay(object),
                  complexity_cutoff = complexityCutoff(preproc_config),
                  housekeeping_cutoff = housekeepingCutoff(preproc_config),
                  expression_cutoff = expressionCutoff(preproc_config),
                  log_base = logBase(preproc_config),
                  scaling_factor = scalingFactor(preproc_config),
                  forced_genes_set = forced_genes_set, use_housekeeping_filter = use_housekeeping_filter)

  # subset the object according to the result of the preprocessing function, basically dropping the low quality cells and lowly expressed genes
  object <- object[rownames(x), colnames(x)]

  # Add the new log TPM assay to the object
  logtpm(object) <- x

  # Drop the tpm assay
  assays(object) <- assays(object)[c("logtpm")]

  return (object)
}

.aggregate_cell_ids <- function(object) {

  cell_ids <- c()
  for (c in childNodes(object))
    cell_ids <- c(cell_ids, colnames(c))

  return (cell_ids)
}
