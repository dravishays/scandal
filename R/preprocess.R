
#'
#' @title Load dataset from file
#'
#' @description Loads a dataset from file
#'
#' @param filename the name of the file containing the dataset (in tab-delimited
#' format).
#' @param drop_cols number of columns that should be dropped from the dataset,
#' i.e. if drop_cols == 3 then columns 1:3 will be dropped. Default is 1.
#' @param rownames_col the column that contains the rownames, i.e. gene symbols
#' or identifiers. Default is 1.
#' @param excluded_samples a vector containing the names of the samples that
#' should be excluded from the returned dataset. Defaults to NULL.
#' @param as_Matrix logical indicating whether the loaded dataset should be returned
#' as an S4 Matrix class (supports sparse representation) or base R matrix type.
#' Default is TRUE.
#'
#' @details
#'
#' @return
#'
#' @examples
#'
#' @author Avishay Spitzer
#'
#' @export
load_dataset <- function(filename, drop_cols = 1, rownames_col = 1, excluded_samples = NULL, as_Matrix = TRUE, verbose = FALSE) {

  stopifnot(is.character(filename), base::file.exists(filename))
  stopifnot(is.integer(drop_cols), is.integer(rownames_col), drop_cols > 0, rownames_col > 0)
  stopifnot(is.null(excluded_samples) | (is.vector(excluded_samples) & is.character(excluded_samples)))
  stopifnot(is.logical(as_Matrix))

  if (isTRUE(verbose))
    message("Loading dataset from ", filename)

  dataset <- utils::read.delim(filename, header = TRUE, sep = '\t', stringsAsFactors = FALSE)

  #tpm_data[, rownames_col][c(11570, 11573)] <- c("MARCH1", "MARCH2")
  rownames(dataset) <- dataset[, rownames_col]

  dataset <- dataset[, -seq_len(drop_cols)]

  colnames(dataset) <- base::gsub("\\.", "-", colnames(dataset))

  if (isTRUE(verbose))
    message("Loaded dataset with ", nrow(dataset), " rows and ", ncol(dataset), " columns")

  if (!is.null(excluded_samples)) {

    dataset <- dataset[, !(.cell2tumor(colnames(dataset)) %in% excluded_samples)]

    if (isTRUE(verbose))
      message("Excluding ", length(excluded_samples), " samples from dataset which now contains ", nrow(dataset), " rows and ", ncol(dataset), " columns")
  }

  if (isTRUE(as_Matrix))
    dataset <- methods::as(dataset, "Matrix")
  else
    dataset <- as.matrix(dataset)

  return (dataset)
}

#'
#' @title Scandal object preprocessing
#'
#' @description Performs preprocessing of ScandalDataSet objects including breaking
#' up the dataset into objects representing the specific samples, filtering out low
#' quality cells and lowly expressed genes and log transforming.
#'
#' @param object a ScandalDataSet object (the underlying object).
#' @param preproc_config_list a named list of containing a \linkS4class{PreprocConfig}
#' object for each sample. The names should correspond to the sample name.
#' @param forced_genes_set a vector of genes that should be included in the final
#' processed object even if their expression is low with the exception of forced
#' genes with absolute count equals to zero which will be filtered out. Default is NULL.
#' @param use_housekeeping_filter should cells with low expression of house-keeping genes
#' be filtered out. Default is FALSE.
#'
#' @details This function is the entry point for preprocessing a \linkS4class{ScandalDataSet}
#' object to allow downstream analysis.
#'
#' The main steps of preprocessing are as follows:
#' \enumerate{
#'   \item Create a \linkS4class{ScandalDataSet} object for each individual (samples
#'   are identified by the names of the \code{preproc_config_list} components).
#'   \item Filtering out low quality cells (cells with low complexity) by summing-up for
#'   each cell (column) the number of genes with count greater than zero and removing
#'   the cells outside the complexity cutoff range configured for the specific sample
#'   in \code{preproc_config_list}.
#'   \item A possible step of filtering out cells with low expression of house-keeping
#'   genes, i.e. genes that are normally highly expressed in most cells (for example,
#'   genes that encode ribosomes). Cells with mean expression of HK genes less than the
#'   housekeeping cutoff configured  for the specific sample in \code{preproc_config_list}
#'   will be removed.
#'   \item Filtering out lowly expressed genes i.e. genes with log2 mean expression less
#'   than the expression cutoff range configured for the specific sample in
#'   \code{preproc_config_list}.
#'   \item Log-transforming the expression data.
#'   \item The \linkS4class{ScandalDataSet} object is added to the \code{childNodes} slot
#'   of underlying \code{object}
#'   \item The underlying \code{object} is preprocessed in the same way described above
#'   beside cell filtering (the cells that passed qaulity control are aggregated from
#'   the childNodes objects) to allow plotting all cells together in downstream analysis.
#' }
#'
#' @return A processed ScandalDataSet object ready for analysis.
#'
#' @examples
#'
#' @author Avishay Spitzer
#'
#' @export
scandal_preprocess <- function(object, preproc_config_list, forced_genes_set = NULL, use_housekeeping_filter = FALSE, verbose = FALSE) {

  stopifnot(!is_scandal_object(object),
            (is.null(forced_genes_set) | is.vector(forced_genes_set)),
            is.logical(use_housekeeping_filter))

  stopifnot(!is.null(preproc_config_list),
            is.list(preproc_config_list),
            base::setequal(nodeIDs(object), names(preproc_config_list)),
            base::all(base::lapply(preproc_config_list, function(x) is_config_object(x)) == TRUE))

  object <- .scandal_preprocess(object, preproc_config_list = preproc_config_list, forced_genes_set = forced_genes_set, use_housekeeping_filter = use_housekeeping_filter, verbose = verbose)

  return (object)
}

#'
#' @title Matrix preprocessing
#'
#' @description Performs preprocessing of a matrix object including breaking
#' up the dataset into objects representing the specific samples, filtering out low
#' quality cells and lowly expressed genes and log transforming.
#'
#' @param x a numeric matrix.
#' @param complexity_cutoff a numeric vector of length 2 representing the lower and
#' upper bounds of complexity (i.e. the number of detected genes per cell).
#' @param expression_cutoff a numeric representing the minimal log2 mean expression
#' per gene below which a gene is considered lowly expressed.
#' @param housekeeping_cutoff a numeric representing the log2 mean expression of
#' house-keeping genes (i.e. genes that are highly expressed in all cells) per
#' cell below which a cell is considered low quality.
#' @param log_base a numeric representing the logarithm base for performing log
#' transformation on the data.
#' @param scaling_factor a numeric representing a scaling factor by which to divide
#' each data point before log transformation.
#' @param pseudo_count a numeric representing the pseudo count added when performing
#' log transformation to avoid taking the log of zero.
#' @param cell_ids a charactyer vector containing IDs of cells that already passed QC.
#' enables bypassing the low-quality cells filtering step.
#' @param gene_ids a charactyer vector containing IDs of genes that already passed QC.
#' enables bypassing the lowly-expressed genes filtering step.
#' @param forced_genes_set a vector of genes that should be included in the final
#' processed object even if their expression is low with the exception of forced
#' genes with absolute count equals to zero which will be filtered out. Default is NULL.
#' @param sample_id an ID of the sample being processed. Used only for printing and hence
#' is not a mandatory parameter. Default is NULL.
#' @param forced_genes_set a vector of genes that should be included in the final
#' processed object even if their expression is low with the exception of forced
#' genes with absolute count equals to zero which will be filtered out. Default is NULL.
#' @param use_housekeeping_filter should cells with low expression of housekeeping genes
#' should be filtered out. Default is FALSE.
#'
#' @details This function is performs the matrix preprocessing steps and is used for
#' preprocessing \linkS4class{ScandalDataSet} objects to allow downstream analysis.
#'
#' The main steps of preprocessing are as follows:
#' \enumerate{
#'   \item Filtering out low quality cells (cells with low complexity) by summing-up for
#'   each cell (column) the number of genes with count greater than zero and removing
#'   the cells outside the complexity cutoff range configured for the specific sample
#'   in \code{complexity_cutoff}.
#'   \item A possible step of filtering out cells with low expression of house-keeping
#'   genes, i.e. genes that are normally highly expressed in most cells (for example,
#'   genes that encode ribosomes). Cells with mean expression of HK genes less than the
#'   housekeeping cutoff configured  for the specific sample in \code{housekeeping_cutoff}
#'   will be removed.
#'   \item Filtering out lowly expressed genes i.e. genes with log2 mean expression less
#'   than the expression cutoff range configured for the specific sample in
#'   \code{expression_cutoff}.
#'   \item Log-transforming the expression data.
#' }
#'
#' @return A processed matrix ready for downstream analysis.
#'
#' @note The function assumes that each column represents a cell and each row represents a gene.
#'
#' @seealso \code{\link{scandal_preprocess}}
#'
#' @author Avishay Spitzer
#'
#' @export
preprocess <- function(x, complexity_cutoff, expression_cutoff, housekeeping_cutoff, log_base, scaling_factor, pseudo_count,
                       sample_id = NULL, cell_ids = NULL, gene_ids = NULL, forced_genes_set = NULL, use_housekeeping_filter = FALSE, verbose = FALSE) {

  stopifnot(is_valid_assay(x),
            (is.vector(complexity_cutoff) & is.numeric(complexity_cutoff)),
            is.numeric(housekeeping_cutoff),
            is.numeric(expression_cutoff),
            is.numeric(log_base),
            is.numeric(scaling_factor),
            is.numeric(pseudo_count),
            (is.null(forced_genes_set) | is.vector(forced_genes_set)),
            (is.null(cell_ids) | is.vector(cell_ids)),
            (is.null(gene_ids) | is.vector(gene_ids)),
            is.logical(use_housekeeping_filter))

  if (isTRUE(verbose))
    message("Preprocessing sample ", sample_id, "...")

  sparcity_before_qc <- length(which(x == 0)) / (dim(x)[1] * dim(x)[2]) * 100

  if (is.null(cell_ids)) {

    if (isTRUE(verbose))
      message("Detecting high quality cells...")

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

  if (is.null(gene_ids)) {

    if (isTRUE(verbose))
      message("Detecting highly expressed genes...")

    x <- x[filter_lowly_expressed_genes(x, expression_cutoff = expression_cutoff, forced_genes_set = forced_genes_set, verbose = verbose), ]

  } else {

    if (isTRUE(verbose))
      message("Skipping gene filtering as gene IDs list was supplied")

    if (!base::all(gene_ids %in% rownames(x)))
      warning("Some (or all) of the supplied gene IDs were not found in the dataset")

    x <- x[rownames(x) %in% gene_ids,]
  }

  if (isTRUE(verbose))
    message(paste0("Log transforming matrix, (base - ", log_base, ", scaling factor - ", scaling_factor, ", pseudo count - ", pseudo_count, ")"))

  x <- log_transform(x, log_base = log_base, scaling_factor = scaling_factor, pseudo_count = pseudo_count, verbose = verbose)

  # if (isTRUE(verbose))
  #   message("Centering matrix...")
  #
  # x <- center_matrix(x, by = "row", method = "mean", verbose = verbose)

  sparcity_after_qc <- length(which(x == 0)) / (dim(x)[1] * dim(x)[2]) * 100

  if (isTRUE(verbose))
    message(sprintf("Sparcity before QC: %.2f%%, after QC: %.2f%%", sparcity_before_qc, sparcity_after_qc))

  if (isTRUE(verbose))
    message("Preprocessing done!")

  return (x)
}

#'
#' @title Log-transforms a numeric object
#'
#' @description Performs log transformation on a numeric object. The data can be
#' scaled first if requested by dividing X by the provided scaling factor.
#'
#' @param x a numeric object (i.e. vector or matrix)
#' @param log_base a positive number representing the base with respect to which
#' the logarithm will be computed. Default is 2
#' @param scaling_factor a positive number by which the data will be scaled, i.e.
#' divided prior to log computation. Default is 1, i.e. no scaling
#' @param pseudo_count a positive number which will be added to the possibly scaled
#' x prior to log computation to avoid taking the logarithm of zero. Default is 1
#'
#' @return A log-transformed numeric object.
#'
#' @author Avishay Spitzer
#'
#' @export
log_transform <- function(x, log_base = 2, scaling_factor = 1, pseudo_count = 1, verbose = FALSE) {

  stopifnot(!is.null(x), is.numeric(x),
            (is.numeric(log_base) & log_base > 0),
            (is.numeric(scaling_factor) & scaling_factor > 0),
            (is.numeric(pseudo_count) & pseudo_count > 0))

  x <- log( (x / scaling_factor) + pseudo_count, base = log_base)

  return (x)
}

#'
#' @describeIn log_transform Reverses the log transformation, i.e \deqn{(log_base^x * scaling_factor) - pseudo_count}
#'
#' @export
reverse_log_transform <- function(x, log_base = 2, scaling_factor = 1, pseudo_count = 1, verbose = FALSE) {

  stopifnot(!is.null(x), is.numeric(x),
            (is.numeric(log_base) & log_base > 0),
            (is.numeric(scaling_factor) & scaling_factor > 0),
            (is.numeric(pseudo_count) & pseudo_count > 0))

  x <- (log_base^x * scaling_factor) - pseudo_count

  return (x)
}

#'
#' @title Compute cell complexity
#'
#' @description Computes the complexity, i.e. the number of genes with count
#' greater than zero for each cell (column).
#'
#' @param x a numeric matrix or Matrix object
#' @param return_sorted should the result be returned sorted. Default is FALSE
#' @param cell_subset a vector of cell IDs for which the complexity should be
#' calculated. Default is NULL
#'
#' @return A named vector of complexity per cell ID.
#'
#' @examples
#' m <- matrix(c(10, 0, 2, 3, 0, 0, 1, 15, 3), nrow = 3, ncol = 3)
#'
#' c <- compute_complexity(m) # 2, 1, 3
#'
#' c <- compute_complexity(m, cell_subset = c(1, 3)) # 2, 3
#'
#' @author Avishay Spitzer
#'
#' @export
compute_complexity <- function(x, return_sorted = FALSE, cell_subset = NULL, verbose = FALSE) {

  stopifnot(is_valid_assay(x), is.logical(return_sorted), (is.null(cell_subset) | is.vector(cell_subset)))

  if (!is.null(cell_subset))
    x <- x[, cell_subset]

  c <- base::apply(x, 2, function(y) base::sum(y != 0))

  if(isTRUE(return_sorted)) {
    c <- base::sort(c)
  }

  return (c)
}

#'
#' @title Centers a matrix
#'
#' @description Centers the mean/median of each row/column of the given matrix
#' around zero.
#'
#' @param x a numeric matrix or Matrix object.
#' @param by either "row" or "col". Default is "row".
#' @param method either "mean" or "median". Default is "mean".
#' @param scale logical indicating whether to scale the rows/columns of x, i.e.
#' divide each row/column by the standard deviation of the row/column. Default
#' is FALSE.
#'
#' @return A matrix with either mean or median of row/column centered around zero.
#'
#' @examples
#' # Center the mean of each row around zero
#' m <- matrix(runif(25, 0, 100), nrow = 5, ncol = 5) # Generate a 5x5 numeric matrix
#' m <- center_matrix(m, by = "row", method = "mean", scale = FALSE) # Center
#' all(rowMeans(m) == 0) # TRUE
#' all(colMeans(m) == 0) # FALSE
#'
#' # Center the median of each column around zero
#' m <- matrix(runif(25, 0, 100), nrow = 5, ncol = 5) # Generate a 5x5 numeric matrix
#' m <- center_matrix(m, by = "col", method = "median", scale = FALSE) # Center
#' all(rowMedians(m) == 0) # FALSE
#' all(colMedians(m) == 0) # TRUE
#'
#' @author Avishay Spitzer
#'
#' @export
center_matrix <- function(x, by = "row", method = "mean", scale = FALSE, verbose = FALSE) {

  stopifnot(is_valid_assay(x), by %in% c("row", "col"), method %in% c("mean", "median"), is.logical(scale))

  center <- .compute(x, by = by, method = method, log_transform_res = FALSE, genes_subset = NULL, verbose = verbose)

  margin <- ifelse(by == "row", 1, 2)

  x <- base::sweep(x, MARGIN = margin, STATS = center, FUN = "-", check.margin = FALSE)

  if (isTRUE(scale)) {
    scale <- .compute(x, by = by, method = "sd", log_transform_res = FALSE, genes_subset = NULL, verbose = verbose)

    x <- sweep(x, MARGIN = margin, STATS = scale, FUN = "/", check.margin=FALSE)
  }

  return (x)
}

# Not exported
filter_low_quality_cells <- function(x, complexity_cutoff, verbose = FALSE) {

  stopifnot(is_valid_assay(x), is.numeric(complexity_cutoff), length(complexity_cutoff) == 2, complexity_cutoff[1] >= 0, complexity_cutoff[2] > complexity_cutoff[1])

  d <- compute_complexity(x, return_sorted = FALSE, cell_subset = NULL, verbose = verbose)

  passQC <- d[(d >= complexity_cutoff[1]) & (d <= complexity_cutoff[2]), drop = FALSE]

  if (isTRUE(verbose))
    message(sprintf("%d cells pre-QC, cell cutoff - lower bound %d [genes], upper bound %d [genes], %d cells passed QC, %2.0f%% of cells dropped",
                    length(d),
                    complexity_cutoff[1],
                    complexity_cutoff[2],
                    length(passQC),
                    (1 - length(passQC)/length(d)) * 100))

  return (names(passQC))
}

# Not exported
filter_low_housekeeping_cells <- function(x, housekeeping_cutoff, verbose = FALSE) {

  stopifnot(is_valid_assay(x), is.numeric(housekeeping_cutoff), housekeeping_cutoff > 0)

  hk_mean_exp <- .compute(x, by = "col", method = "mean", log_transform_res = TRUE, genes_subset = HOUSEKEEPING_GENES_LIST, verbose = verbose)

  passQC <- hk_mean_exp[hk_mean_exp >= housekeeping_cutoff]

  if (isTRUE(verbose))
    message(sprintf("Housekeeping cutoff %d, %d cells pre-QC, %d cells passed QC, %2.0f%% of cells dropped",
                    housekeeping_cutoff,
                    length(hk_mean_exp),
                    length(passQC),
                    (1 - length(passQC)/length(hk_mean_exp)) * 100))

  return (names(passQC))
}

# Not exported
filter_lowly_expressed_genes <- function(x, expression_cutoff, forced_genes_set = NULL, verbose = FALSE) {

  stopifnot(is_valid_assay(x), is.numeric(expression_cutoff), expression_cutoff > 0, (is.null(forced_genes_set) | is.vector(forced_genes_set)))

  gene_counts <- .compute(x, by = "row", method = "mean", log_transform_res = TRUE, genes_subset = NULL, verbose = verbose)

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
#' @title Quality control statistics
#'
#' @description This function returns a \code{DataFrame} containing statistics
#' gathered while performing the preprocessing procedure including cell and gene
#' drop rates.
#'
#' @param object a \code{ScandalDataSet} object
#'
#' @export
scandal_quality_control_stats <- function(object) {
  stopifnot(!is_scandal_object(object))

  res <- do.call(rbind, lapply(qualityControl(object), function(x) statsQC(x)))

  return(res)
}

#' @export
scandal_plot_qc_metrics <- function(object, preproc_config_list, show_plot = TRUE, save_to_file = TRUE) {

  stopifnot(!is_scandal_object(object), is.logical(show_plot), is.logical(save_to_file))

  stopifnot(!is.null(preproc_config_list),
            is.list(preproc_config_list),
            base::setequal(nodeIDs(object), names(preproc_config_list)),
            base::all(base::lapply(preproc_config_list, function(x) is_config_object(x)) == TRUE))

  project_id <- projectID(object)

  for(nid in nodeIDs(object)) {
    nconf <- preproc_config_list[[nid]]
    ndata <- assay(object)[, .subset_cells(colnames(object), nid), drop = FALSE]

    plot_cell_complexity_distribution(ndata, complexity_cutoff = complexityCutoff(nconf), node_id = nid, project_id = project_id, show_plot = show_plot, save_to_file = save_to_file)

    plot_mean_housekeeping_expression(ndata, housekeeping_cutoff = housekeepingCutoff(nconf), node_id = nid, project_id = project_id, show_plot = show_plot, save_to_file = save_to_file)

    plot_mean_expression_frequency(ndata, expression_cutoff = expressionCutoff(nconf), node_id = nid, project_id = project_id, show_plot = show_plot, save_to_file = save_to_file)
  }
}

#' @export
plot_cell_complexity_distribution <- function(x, complexity_cutoff, node_id, project_id, show_plot = TRUE, save_to_file = TRUE) {

  stopifnot(is_valid_assay(x), is.numeric(complexity_cutoff), is.character(node_id), is.character(project_id), is.logical(show_plot), is.logical(save_to_file))

  complexity <- compute_complexity(x)

  p <- generate_scatter_plot(y_data = complexity, title = paste0(node_id, " - Cell complexity distribution"), xlab = "Cells", ylab = "Complexity", plot_ordered = TRUE)

  p <- p +
    ggplot2::geom_hline(yintercept = complexity_cutoff[1], linetype = "dashed", color = "red") +
    ggplot2::geom_hline(yintercept = complexity_cutoff[2], linetype = "dashed", color = "red")

  scandal_plot(p, show_plot = show_plot, project_dir = project_id, save_to_file = save_to_file, filename = paste0(node_id, "_complexity.png"), dirname = "QC")
}

#' @export
plot_mean_housekeeping_expression <- function(x, housekeeping_cutoff, node_id, project_id, show_plot = TRUE, save_to_file = TRUE) {

  stopifnot(is_valid_assay(x), is.numeric(housekeeping_cutoff), is.character(node_id), is.character(project_id), is.logical(show_plot), is.logical(save_to_file))

  hk_mean_exp <- .compute(x, by = "col", method = "mean", log_transform_res = TRUE, genes_subset = HOUSEKEEPING_GENES_LIST)

  p <- generate_scatter_plot(y_data = hk_mean_exp, title = paste0(node_id, " - Mean expression of housekeeping genes distribution"), xlab = "Cells", ylab = "Mean expression [log2]", plot_ordered = TRUE)

  p <- p + ggplot2::geom_hline(yintercept = housekeeping_cutoff, linetype = "dashed", color = "red")

  scandal_plot(p, show_plot = show_plot, save_to_file = save_to_file, project_dir = project_id, filename = paste0(node_id, "_mean_hk_expression.png"), dirname = "QC")
}

#' @export
plot_mean_expression_frequency <- function(x, expression_cutoff, node_id, project_id, show_plot = TRUE, save_to_file = TRUE) {

  stopifnot(is_valid_assay(x), is.numeric(expression_cutoff), is.character(node_id), is.character(project_id), is.logical(show_plot), is.logical(save_to_file))

  gene_exp <- .compute(x, by = "row", method = "mean", log_transform_res = TRUE)

  p <- generate_histogram_plot(data = gene_exp, title = paste0(node_id, " - Mean gene expression frequency"), xlab = "Mean expression [log2]", ylab = "Frequency")

  p <- p + ggplot2::geom_vline(xintercept = expression_cutoff, linetype = "dashed", color = "red")

  scandal_plot(p, show_plot = show_plot, save_to_file = save_to_file, project_dir = project_id, filename = paste0(node_id, "_mean_gene_expression.png"), dirname = "QC")
}

#' @export
scandal_plot_complexity_distribution <- function(object, show_plot = TRUE, save_to_file = TRUE) {

  stopifnot(!is_scandal_object(object), is.logical(show_plot), is.logical(save_to_file))

  complexity_pre  <- compute_complexity(unprocessedData(object))
  complexity_post <- compute_complexity(assay(object))

  p1 <- generate_whiskers_plot(data = complexity_pre, title = paste0(nodeID(object), " - Complexity distribution pre-QC"), labels = .cell2tumor(colnames(unprocessedData(object))), xlab = NULL, ylab = "Complexity")
  p2 <- generate_whiskers_plot(data = complexity_post, title = paste0(nodeID(object), " - Complexity distribution post-QC"), labels = .cell2tumor(colnames(object)), xlab = NULL, ylab = "Complexity")

  p <- p1 + p2

  scandal_plot(p, show_plot = show_plot, save_to_file = save_to_file, project_dir = projectID(object), filename = paste0(nodeID(object), "_complexity_per_tumor_whiskers.png"), dirname = "QC")
}

#'
#' @title
#'
#' @description
#'
#' @param object a \code{ScandalDataSet} object
#' @param node_id the ID of the node to be inspected
#'
#' @details
#'
#' @examples
#'
#' @author Avishay Spitzer
#'
#' @export
scandal_inspect_node <- function(object, node_id, verbose = FALSE) {

  stopifnot(!is_scandal_object(object))
  stopifnot(!is.null(node_id), is.character(node_id), length(node_id) == 1)

  if (!(node_id %in% nodeIDs(object)))
    stop("Node ", node_id, " not found")

  node <- inspectNode(object, nodeID = node_id)

  preproc_config <- preprocConfig(node)

  x <- unprocessedData(node)

  x <- preprocess(as.matrix(x),
                  complexity_cutoff = complexityCutoff(preproc_config),
                  expression_cutoff = expressionCutoff(preproc_config),
                  housekeeping_cutoff = housekeepingCutoff(preproc_config),
                  log_base = logBase(preproc_config),
                  scaling_factor = scalingFactor(preproc_config),
                  pseudo_count = pseudoCount(preproc_config),
                  sample_id = nodeID(node),
                  cell_ids = colnames(node),
                  gene_ids = rownames(node),
                  forced_genes_set = NULL, use_housekeeping_filter = FALSE, verbose = verbose)

  node <- .assign_assay(node, x)

  return (node)
}

.subset_cells <- function(cell_names, sample_name) cell_names[which(.cell2tumor(cell_names) %in% sample_name)]

.cell2tumor <- function(cell_ids) base::gsub("-.*", "", cell_ids)

.sparsity <- function(x) length(which(x == 0)) / (dim(x)[1] * dim(x)[2]) * 100

.qc_stats <- function(x_pre, x_post, nid) {

  cells_pre_qc <- ncol(x_pre)
  cells_post_qc <- ncol(x_post)
  genes_pre_qc <- nrow(x_pre)
  genes_post_qc <- nrow(x_post)

  DataFrame("Cells pre-QC" = cells_pre_qc,
            "Cells post-QC" = cells_post_qc,
            #"Cells dropped" = paste0((1 - (cells_post_qc / cells_pre_qc)) * 100, "%%"),
            "Cells dropped" = sprintf("%.2f%%", (1 - (cells_post_qc / cells_pre_qc)) * 100),
            "Genes pre-QC" = genes_pre_qc,
            "Genes post-QC" = genes_post_qc,
            #"Genes dropped" = paste0((1 - (genes_post_qc / genes_pre_qc)) * 100, "%%"),
            "Genes dropped" = sprintf("%.2f%%", (1 - (genes_post_qc / genes_pre_qc)) * 100),
            #"Sparsity pre QC" = paste0(.sparsity(x_pre), "%"),
            #"Sparsity post QC" = paste0(.sparsity(x_post), "%"),
            "Sparsity pre QC" = sprintf("%.2f%%",.sparsity(x_pre)),
            "Sparsity post QC" = sprintf("%.2f%%",.sparsity(x_post)),
            row.names = nid)
}

.scandal_preprocess <- function(object, preproc_config_list, forced_genes_set = NULL, use_housekeeping_filter = FALSE, verbose = FALSE) {

  for(nid in nodeIDs(object)) {
    nconf <- preproc_config_list[[nid]]

    x_pre <- as.matrix(assay(object)[, .subset_cells(colnames(object), nid), drop = FALSE])

    x <- preprocess(x_pre,
                    complexity_cutoff = complexityCutoff(nconf),
                    expression_cutoff = expressionCutoff(nconf),
                    housekeeping_cutoff = housekeepingCutoff(nconf),
                    log_base = logBase(nconf),
                    scaling_factor = scalingFactor(nconf),
                    pseudo_count = pseudoCount(nconf),
                    sample_id = nid, cell_ids = NULL,
                    forced_genes_set = forced_genes_set, use_housekeeping_filter = use_housekeeping_filter, verbose = verbose)

    stats_qc <- .qc_stats(x_pre, x, nid)

    qc <- QCResults(nconf, cellIDs = colnames(x), geneIDs = rownames(x), statsQC = stats_qc)

    qualityControl(object)[[nid]] <- qc
  }

  # Coerce to base R matrix
  x <- as.matrix(assay(object))

  # Extract the preprocessing configuration object
  preproc_config <- preprocConfig(object)

  # Call the matrix preprocessing function
  x <- preprocess(x,
                  complexity_cutoff = complexityCutoff(preproc_config),
                  expression_cutoff = expressionCutoff(preproc_config),
                  housekeeping_cutoff = housekeepingCutoff(preproc_config),
                  log_base = logBase(preproc_config),
                  scaling_factor = scalingFactor(preproc_config),
                  pseudo_count = pseudoCount(preproc_config),
                  sample_id = nodeID(object),
                  cell_ids = .aggregate_cell_ids(object),
                  gene_ids = .aggregate_gene_ids(object),
                  forced_genes_set = forced_genes_set, use_housekeeping_filter = use_housekeeping_filter, verbose = verbose)

  stats_qc <- .qc_stats(as.matrix(assay(object)), x, nodeID(object))

  qc <- QCResults(preproc_config, cellIDs = colnames(x), geneIDs = rownames(x), statsQC = stats_qc)

  qualityControl(object)[[nodeID(object)]] <- qc

  # re-assign the matrix after preprocessing
  object <- .assign_assay(object, x)

  return (object)
}

.aggregate_cell_ids <- function(object) {

  cell_ids <- lapply(qualityControl(object), function(x) cellIDs(x))
  cell_ids <- base::unname(base::unlist(cell_ids))
  cell_ids <- base::unique(cell_ids)

  return (cell_ids)
}

.aggregate_gene_ids <- function(object) {

  gene_ids <- lapply(qualityControl(object), function(x) geneIDs(x))
  gene_ids <- base::unname(base::unlist(gene_ids))
  gene_ids <- base::unique(gene_ids)

  return (gene_ids)
}

.assign_assay <- function(object, x) {

  preproc_config <- preprocConfig(object)

  # subset the object according to the result of the preprocessing function, basically dropping the low quality cells and lowly expressed genes
  object <- object[rownames(x), colnames(x)]

  # Coerce to Matrix object to benefit from sparse representation
  if(isTRUE(typeMatrix(preproc_config)))
    x <- as(x, "Matrix")

  # Add the new log TPM assay to the object
  logtpm(object) <- x

  # Drop the tpm assay
  assays(object) <- assays(object)[c("logtpm")]

  return (object)
}

# Internal function that can run different computations given a method on either rows or columns with the posibility to subset them
.compute <- function(x, by = "row", method = "mean", log_transform_res = FALSE, genes_subset = NULL, verbose = FALSE) {

  stopifnot(is_valid_assay(x), by %in% c("row", "col"), is.logical(log_transform_res), (is.null(genes_subset) | is.vector(genes_subset)))

  if (!is.null(genes_subset))
    x <- x[rownames(x) %in% genes_subset, ]

  func <- base::match.fun(method)

  x <- base::apply(x, ifelse(by == "row", 1, 2), func)

  if (isTRUE(log_transform_res))
    x <- base::log2(x + 1)

  return (x)
}
