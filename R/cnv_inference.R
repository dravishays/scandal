
### =========================================================================
### CNV matrix computation and plotting
### -------------------------------------------------------------------------
###

#'
#' @title CNV inference
#'
#' @description This function infers CNVs (chromosomal copy-number variations) from the
#' single-cell expression data. CNV inference is the main method of the scandal framework
#' for classifying malignant and non-malignant cells.
#'
#' @param object a \linkS4class{ScandalDataSet} object.
#' @param gene_positions_table a data frame containing all genes ordered by the position
#' of the gene on the chromosome and by the order of the chromosomes. The data frame should
#' contain the column names "Gene" and "CHR".
#' @param reference_cells a named vector of the cluster assignments of the reference cells.
#' The names should correspond to the cell IDs of the reference (non-malignant) cells. The
#' CNV matrix can be computed without a reference (with \code{reference=NULL}) but this is
#' not recommended as downstream comoutations using the inferred CNV matrix will be less
#' reliable.
#' @param max_genes maximal number of genes to use for computing the CNV matrix. Default
#' is 5000.
#' @param expression_limits a numeric vector with two elements representing the upper
#' and lower values with which to bound the centered expression matrix prior to
#' calculating the CNV matrix. This blunts the effect of noisy genes. Defaut is (-3, 3).
#' @param window number of genes to consider when calculating the running mean. Default
#' is a window of 100 genes.
#' @param scaling_factor a small constant by which to increase the calculated
#' (-BM, +BM) interval to compensate for possible noise. Default is 0.2.
#' @param initial_centering direction of centering the expression matrix (row-wise or
#' col-wise) prior to computing the CNV matrix. Accepts either strings "row" or "col",
#' default is "col".
#' @param base_metric a metric to use for calculating the (-BM, + BM) interval. Accepts
#' either strings "mean" or "median", default is "median".
#' @param verbose suppresses all messages from this function. Default is FALSE.
#'
#' @details The CNV algorithm is as follows:\cr
#' Preprocessing steps:
#' \enumerate{
#'   \item Compute mean expression for each gene (log2[mean(TPM) + 1])
#'   \item Keep the \code{max_genes} highest expressed genes
#'   \item Order the rows (genes) of the expression matrix according to chromosomal
#'   position
#'   \item Log-transform the expression matrix
#'   \item Mean-center of the expression matrix in the \code{initial_centering} direction
#'   \item Bound the expression matrix according to the \code{expression_limits}
#' }
#' \cr
#'
#' @return Returns the \linkS4class{ScandalDataSet} object with CNV matrix in the
#' "cnv" element of the reducedDim slot (accessible by reducedDim(object, "cnv")). Note
#' that the matrix is stored with cell IDs as row names and gene IDs as column names.
#'
#' @seealso The CNV inference method was defined and developed by **Dr. Itay Tirosh**
#' during his time at the *Broad Institute* and published in several high-impact papers
#' including the following paper from *Cell*:
#' https://www.cell.com/cell/fulltext/S0092-8674(17)31270-9.
#'
#' @author Avishay Spitzer
#'
#' @export
scandal_cnv_infer <- function(object, gene_positions_table, reference_cells,
                              max_genes = 5000, expression_limits = c(-3, 3), window = 100, scaling_factor = 0.2,
                              initial_centering = "col", base_metric = "median",
                              verbose = FALSE) {

  stopifnot(is_scandal_object(object))
  stopifnot(!is.null(gene_positions_table), is.data.frame(gene_positions_table), c("Gene", "CHR") %in% colnames(gene_positions_table))
  stopifnot(!is.null(reference_cells), is.vector(reference_cells),!is.null(names(reference_cells)), all(names(reference_cells) %in% colnames(object)) == TRUE,
            is.character(reference_cells) | is.factor(reference_cells))
  stopifnot(is.numeric(max_genes))
  stopifnot(!is.null(expression_limits), is.vector(expression_limits), is.numeric(expression_limits), length(expression_limits) == 2)
  stopifnot(is.numeric(window))
  stopifnot(is.numeric(scaling_factor))
  stopifnot(!is.null(initial_centering), initial_centering %in% c("col", "row"))
  stopifnot(!is.null(base_metric), base_metric %in% c("mean", "median"))
  stopifnot(is.logical(verbose))

  x <- as.matrix(unprocessedData(object))[, colnames(object)]

  cnv_matrix <- scandal_cnv_compute_matrix(x = x,
                                           gene_positions_table = gene_positions_table,
                                           reference_cells = reference_cells,
                                           max_genes = max_genes,
                                           expression_limits = expression_limits,
                                           window = window,
                                           scaling_factor = scaling_factor,
                                           initial_centering = initial_centering,
                                           base_metric = base_metric,
                                           verbose = verbose)

  if (is.null(cnv_matrix))
    stop("An error has occured while computing the CNV matrix")

  reducedDim(object, "cnv") <- t(cnv_matrix)

  return (object)
}

#'
#' @title CNV matrix computation
#'
#' @description This function
#'
#' @param x
#' @param gene_positions_table
#' @param reference_cells
#' @param max_genes
#' @param expression_limits
#' @param window
#' @param scaling_factor
#' @param initial_centering
#' @param base_metric
#' @param verbose
#'
#' @details
#'
#' @return
#'
#' @author Avishay Spitzer
#'
#' @export
scandal_cnv_compute_matrix <- function(x, gene_positions_table, reference_cells,
                                       max_genes = 5000, expression_limits = c(-3, 3), window = 100, scaling_factor = 0.2,
                                       initial_centering = "col", base_metric = "median",
                                       verbose = FALSE) {

  stopifnot(is_valid_assay(x))
  stopifnot(!is.null(gene_positions_table), is.data.frame(gene_positions_table), c("Gene", "CHR") %in% colnames(gene_positions_table))
  stopifnot(!is.null(reference_cells), is.vector(reference_cells),!is.null(names(reference_cells)), all(names(reference_cells) %in% colnames(x)) == TRUE,
            is.character(reference_cells) | is.factor(reference_cells))
  stopifnot(is.numeric(max_genes))
  stopifnot(!is.null(expression_limits), is.vector(expression_limits), is.numeric(expression_limits), length(expression_limits) == 2)
  stopifnot(is.numeric(window))
  stopifnot(is.numeric(scaling_factor))
  stopifnot(!is.null(initial_centering), initial_centering %in% c("col", "row"))
  stopifnot(!is.null(base_metric), base_metric %in% c("mean", "median"))
  stopifnot(is.logical(verbose))

  x <- as.matrix(x)

  cnv_matrix <- .scandal_compute_cnv_matrix(x = x,
                                            gene_pos_tbl = gene_positions_table,
                                            reference_cells = reference_cells,
                                            max_genes = max_genes,
                                            expression_limits = expression_limits,
                                            window = window,
                                            scaling_factor = scaling_factor,
                                            initial_centering = initial_centering,
                                            base_metric = base_metric,
                                            verbose = verbose)

  if (is.null(cnv_matrix))
    stop("An error has occured while computing the CNV matrix")

  return (cnv_matrix)
}

#' @importFrom ComplexHeatmap Heatmap rowAnnotation HeatmapAnnotation draw decorate_heatmap_body
#' @importFrom circlize colorRamp2
#' @importFrom grDevices dev.off jpeg
#' @export
scandal_cnv_plot <- function(object, cnv_matrix = NULL, gene_positions_table, reference_cells = NULL, show_reference = FALSE,
                             row_annotation = NULL, row_annotation_cols = NULL, low = "dodgerblue", mid = "white", high = "red",
                             save_to_file = TRUE, show_plot = TRUE, filename = NULL, dirname = ".", width = 1600, height = 1000, quality = 100, verbose = FALSE,
                             cluster_rows = FALSE, ...) {

  if (!is.null(object))
    cnv_matrix <- reducedDim(object, "cnv")

  if (is.null(cnv_matrix))
    stop("CNV matrix not found, did you forget to run scandal_infer_cnv first?")

  if (isTRUE(verbose))
    message(sprintf("Plotting CNV: data range (%.2f, %.2f), limit set to (-1, 1)", min(cnv_matrix), max(cnv_matrix)))

  cell_ids <- NULL
  if (!is.null(reference_cells)) {

    malignant_cells <- rownames(cnv_matrix)[!rownames(cnv_matrix) %in% reference_cells]

    if (isTRUE(show_reference))
      cell_ids <- c(malignant_cells, reference_cells)
    else
      cell_ids <- malignant_cells
  } else
    cell_ids <- rownames(cnv_matrix)

  plot_name <- paste0("cnv", sample(1:10^9, 1))

  ra <- NULL
  if (!is.null(row_annotation)) {
    rann <- as.data.frame(row_annotation)[cell_ids, ]

    if (is.null(row_annotation_cols))
      ra <- rowAnnotation(df = rann, show_annotation_name = FALSE)
    else
      ra <- rowAnnotation(df = rann, col = row_annotation_cols, show_annotation_name = FALSE)
  }

  h <- Heatmap(cnv_matrix[cell_ids, ], name = plot_name,
               col = colorRamp2(c(-1, 0, 1), c(low, mid, high)),
               show_row_names = FALSE, show_column_names = FALSE, cluster_columns = FALSE, cluster_rows = cluster_rows,
               left_annotation = ra,
               heatmap_legend_param = list(at = c(-1, -0.5, 0, 0.5, 1), title = "Inferred CNV\n(log2 ratio)"), ...)

  if (isTRUE(show_plot)) {
    draw(h)

    decorate_heatmap_body(plot_name, .draw_grid(cnv_matrix, gene_positions_table))
  }

  if (isTRUE(save_to_file)) {

    if (is.null(filename))
      filename <- paste0(plot_name, ".jpeg")

    file <- paste0(dirname, "/", filename)

    jpeg(filename = file, width = width, height = height, quality = quality)
    draw(h)
    decorate_heatmap_body(plot_name, .draw_grid(cnv_matrix, gene_positions_table))
    dev.off()
  }

  invisible(h)
}

### =========================================================================
### CNV scores computation, classification and plotting
### -------------------------------------------------------------------------
###

#'
#' @title Compute the cell CNV score
#'
#' @description Computes the CNV score for each cell which is composed of the CNV
#' signal (the sum of squares of computed CNVs) and CNV correlation (the pearson's
#' correlation between the cell CNV profile and the tumor's CNV profile).
#'
#' @param object
#' @param cnv_matrix
#' @param gene_positions_table
#' @param hotspot_threshold
#' @param verbose
#'
#' @return A \link{tibble} with three columns:
#' \enumerate{
#'   \item CellID
#'   \item Signal
#'   \item Correlation
#' }
#'
#' @author Avishay Spitzer
#'
#' @importFrom stats cor quantile
#' @importFrom tibble tibble
#' @export
scandal_cnv_compute_scores <- function(object, cnv_matrix = NULL, gene_positions_table, hotspot_threshold = .9, verbose = FALSE) {

  stopifnot(!(is.null(object) & is.null(cnv_matrix)))
  stopifnot(is.null(object) | is_scandal_object(object))

  if (!is.null(object))
    cnv_matrix <- reducedDim(object, "cnv")

  if (is.null(cnv_matrix))
    stop("CNV matrix not found, did you forget to run scandal_infer_cnv first?")

  cnv_matrix <- t(cnv_matrix)

  # Square the CNV matrix
  tmp <- cnv_matrix^2

  # Calculate the hotspots - areas with high signal
  pos_mean <- rowMeans(tmp)
  hotspots <- names(pos_mean[pos_mean >= quantile(pos_mean, hotspot_threshold)])
  tmp <- tmp[hotspots, ]

  # Compute the mean squared CNV signal per
  cnv_signal <- colMeans(tmp)

  # Compute the tumor CNV profile
  tumor_cnv0 <- rowMeans(cnv_matrix)

  # Compute the correlation between each cell and the tumor CNV profile
  cnv_correlation <- apply(cnv_matrix, 2, cor, y = tumor_cnv0, method = "pearson")

  # Return the result as a data frame with two variables
  res <- tibble(CellID = colnames(cnv_matrix), Signal = cnv_signal, Correlation = cnv_correlation)

  return (res)
}

#'
#' @title Classify CNV scores
#'
#' @description This function classifies the CNV scores according to the supplied signal and correlation thresholds.
#'
#' @param cnv_scores
#' @param signal_threshold
#' @param correlation_threshold
#' @param verbose
#'
#' @return
#'
#' @seealso \link{scandal_cnv_compute_scores}
#'
#' @author Avishay Spitzer
#'
#' @export
scandal_cnv_classify_scores <- function(cnv_scores, signal_threshold = 0.05, correlation_threshold = 0.05, verbose = FALSE) {

  cnv_detected <- rep("Non-classifiable", nrow(cnv_scores))

  cnv_detected[which(cnv_scores$Signal >= signal_threshold & cnv_scores$Correlation >= correlation_threshold)] <- "Detected"
  cnv_detected[which(cnv_scores$Signal <  signal_threshold & cnv_scores$Correlation <  correlation_threshold)] <- "Not detected"
  cnv_detected[which(cnv_scores$Signal <  signal_threshold & cnv_scores$Correlation >= correlation_threshold)] <- "Low signal"
  cnv_detected[which(cnv_scores$Signal >= signal_threshold & cnv_scores$Correlation <  correlation_threshold)] <- "Low correlation"

  cnv_scores$CNVDetected <- cnv_detected

  return (cnv_scores)
}

#' @importFrom dplyr %>% group_by summarise n mutate filter select pull left_join
#' @importFrom tibble tibble
#' @importFrom stats setNames
#' @export
scandal_cnv_classify_cells <- function(cnv_scores, clusters, cnv_matrix = NULL, cell_class_init = NULL, min_cluster_cc_freq = .5, return_all = FALSE, verbose = FALSE) {

  # Simple classification - leaving the low signal and low correlation cells as "Unresolved"
  if (is.null(clusters))
    data <- .classify_cells_simple(data = cnv_scores)
  else # Complex classifiaction using cluster correlation
    data <- .classify_cells_cluster_correlation(data = cnv_scores,
                                                clusters = clusters,
                                                cnv_matrix = cnv_matrix,
                                                cell_class_init = cell_class_init,
                                                min_cluster_cc_freq = min_cluster_cc_freq,
                                                verbose = verbose)

  if (isTRUE(return_all))
    return (data)
  else
    return (data$Malignant)
}

#' @export
scandal_cnv_get_malignant_cells <- function(cnv_scores) {

  res <- (cnv_scores %>%
            select(CellID, Malignant) %>%
            filter(Malignant == "Malignant"))$CellID

  return (res)
}

#' @export
scandal_cnv_scores_plot <- function(cnv_scores, signal_threshold = 0.05, correlation_threshold = 0.5, title = NULL, labels = NULL, name = NULL, verbose = FALSE) {

  #cnv_scores <- scandal_cnv_classify_scores(cnv_scores, signal_threshold = signal_threshold, correlation_threshold = correlation_threshold)

  #cnv_detected <- cnv_scores$CNVDetected

  #non_classifiable <- (length(which(cnv_detected == "Low signal" | cnv_detected == "Low correlation")) / length(cnv_detected)) * 100

  p <- scandal_scatter_plot(x = cnv_scores$Correlation, y = cnv_scores$Signal, labels = labels, color_legend_name = name,
                            title = title, xlab = "CNV Correlation", ylab = "CNV Signal", plot_ordered = FALSE) +
        geom_vline(xintercept = correlation_threshold) +
        geom_hline(yintercept = signal_threshold) #+
    #labs(caption = sprintf("%.2f%% of cells are non-classifiable", non_classifiable))

  return (p)
}

### =========================================================================
### Exported utility function
### -------------------------------------------------------------------------
###

#' @importFrom utils read.delim
#' @export
load_gene_pos_file <- function(filename = "gencode_v19_gene_pos.txt") {

  gencode_v19_gene_pos <- read.delim(filename, header = TRUE, sep = '\t')
  colnames(gencode_v19_gene_pos) <- c("Gene", "CHR", "Start", "End")
  gencode_v19_gene_pos$Gene <- as.character(gencode_v19_gene_pos$Gene)
  gencode_v19_gene_pos$CHR <- as.character(gencode_v19_gene_pos$CHR)

  return (gencode_v19_gene_pos)
}

#' @export
count_genes_per_chromosome <- function(cnv_matrix, gene_pos_tbl) {

  chr_pos_tbl <- table(gene_pos_tbl[gene_pos_tbl$Gene %in% colnames(cnv_matrix), "CHR"])[CHRs]

  return (chr_pos_tbl)
}

### =========================================================================
### Internal functions
### -------------------------------------------------------------------------
###

#CHRs <- c("chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8", "chr9", "chr10", "chr11", "chr12", "chr13", "chr14", "chr15", "chr16", "chr17", "chr18", "chr19", "chr20", "chr21", "chr22", "chrX", "chrY")
CHRs <- c("chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8", "chr9", "chr10", "chr11", "chr12", "chr13", "chr14", "chr15", "chr16", "chr17", "chr18", "chr19", "chr20", "chr21", "chr22", "chrX")

.scandal_compute_cnv_matrix <- function(x, gene_pos_tbl, reference_cells, max_genes, expression_limits, window, scaling_factor, initial_centering, base_metric, verbose) {

  cnv_matrix <- .prepare_cnv_matrix(x,
                                    gene_pos_tbl = gene_pos_tbl,
                                    max_genes = max_genes,
                                    expression_limits = expression_limits,
                                    initial_centering = initial_centering,
                                    verbose = verbose)

  gene_pos_tbl <- gene_pos_tbl[gene_pos_tbl$Gene %in% rownames(cnv_matrix), ]

  if (!is.null(reference_cells)) {
    m_ref <- cnv_matrix[, names(reference_cells)]
    m_mal <- cnv_matrix[, !(colnames(cnv_matrix) %in% names(reference_cells))]
  } else
    m_mal <- cnv_matrix

  cnv_ref <- NULL
  cnv_mal <- NULL
  genes_per_chr <- replicate(length(CHRs), 0)
  names(genes_per_chr) <- CHRs

  for (c in CHRs) {

    if (isTRUE(verbose))
      message(paste0("Calculate CNV matrix for ", c))

    chr_pos_tbl <- gene_pos_tbl[as.character(gene_pos_tbl$CHR) == c, ]

    if (!is.null(reference_cells))
      m1 <- m_ref[rownames(m_ref) %in% chr_pos_tbl$Gene, ]

    m2 <- m_mal[rownames(m_mal) %in% chr_pos_tbl$Gene, ]

    if (!is.null(reference_cells))
      m1 <- .calc_cnv_matrix(m1, window = window, verbose = verbose)

    m2 <- .calc_cnv_matrix(m2, window = window, verbose = verbose)

    if (!is.null(reference_cells))
      cnv_ref <- rbind(cnv_ref, m1)

    cnv_mal <- rbind(cnv_mal, m2)

    genes_per_chr[c] <- nrow(m2)
  }

  if (!is.null(reference_cells))
    cnv_matrix <- cbind(cnv_ref, cnv_mal)
  else
    cnv_matrix <- cnv_mal

  if (isTRUE(verbose))
    message(paste0("Median-centering CNV matrix cellwise"))

  cnv_matrix <- center_matrix(cnv_matrix, by = "col", method = "median", scale = FALSE, verbose = verbose)

  if (!is.null(reference_cells)) {

    if (isTRUE(verbose))
      message(paste0("Calculating CNV score"))

    cnv_matrix <- .assign_cnv_score(cnv_matrix, ref_group = reference_cells, scaling_factor = scaling_factor, base_metric = base_metric, verbose = verbose)
  }

  cnv_matrix <- cnv_matrix[, colnames(x)]

  if (isTRUE(verbose))
    message(paste0("CNV matrix created!"))

  return (cnv_matrix)
}

.prepare_cnv_matrix <- function(cnv_matrix, gene_pos_tbl, max_genes, expression_limits, initial_centering, verbose) {

  if (isTRUE(verbose))
    message("Preprocessing CNV matrix...")

  if (isTRUE(verbose))
    message(paste0("Expression matrix contains ", ncol(cnv_matrix), " high quality cells"))

  mean_exp <- log2(rowMeans(cnv_matrix) + 1)

  mean_exp <- mean_exp[order(mean_exp, decreasing = TRUE)]
  passQC <- head(mean_exp, n = min(max_genes, length(mean_exp)))

  cnv_matrix <- cnv_matrix[names(passQC), ]

  if (isTRUE(verbose))
    message(paste0("Keeping ", min(max_genes, length(mean_exp)), " highest expressed genes"))

  gene_pos_tbl <- gene_pos_tbl[gene_pos_tbl$Gene %in% rownames(cnv_matrix), ]

  # Reorder the expression matrix according to the positions of the genes in each chromosome
  if (isTRUE(verbose))
    message(paste0("Reordering genes according to chromosomal location"))

  cnv_matrix <- cnv_matrix[gene_pos_tbl$Gene, ]

  # Log-transform the matrix
  if (isTRUE(verbose))
    message(paste0("Log2 expression matrix [TPM/10]"))

  cnv_matrix <- log_transform(cnv_matrix)

  # Center the matrix
  if (isTRUE(verbose))
    message(sprintf("Mean-centering expression matrix %s-wise", initial_centering))

  cnv_matrix <- center_matrix(cnv_matrix, by = initial_centering, method = "mean", scale = FALSE, verbose = verbose)

  # Bound the expression matrix to reduce the effect of specific highly expressed genes on the mean expression
  if (!is.null(expression_limits)) {

    if (isTRUE(verbose))
      message(paste0("Expresion matrix is in the range (", min(cnv_matrix), ", ", max(cnv_matrix), ")"))

    if (isTRUE(verbose))
      message(paste0("Limiting matrix to values betwwen (", expression_limits[1], ", ", expression_limits[2], ")"))

    cnv_matrix <-.bound_matrix(cnv_matrix, expression_limits = expression_limits, verbose = verbose)

    if (max(cnv_matrix) > expression_limits[2] | min(cnv_matrix) < expression_limits[1])
      stop("Matrix limiting failed")
  }

  if (isTRUE(verbose))
    message(paste0("Preprocessing done!"))

  return (cnv_matrix)
}

.bound_matrix <- function(cnv_matrix, expression_limits, verbose) {

  cnv_matrix[cnv_matrix > expression_limits[2]] <- expression_limits[2]
  cnv_matrix[cnv_matrix < expression_limits[1]] <- expression_limits[1]

  return (cnv_matrix)
}

#' @importFrom caTools runmean
.calc_cnv_matrix = function(m, window, verbose) {

  if (window > nrow(m)) {

    if (isTRUE(verbose))
      message(sprintf("Window too large - setting to %d", nrow(m)))

    window <- nrow(m)
  }

  m <- (2^m)*10 - 1

  movavg_m <- runmean(m, k = window, endrule = "mean")

  movavg_m <- log2(movavg_m + 1)

  rownames(movavg_m) <- rownames(m)
  colnames(movavg_m) <- colnames(m)

  return (movavg_m)
}

#' @importFrom stats median
.assign_cnv_score <- function(m, ref_group, scaling_factor, base_metric, verbose) {

  calc_score <- function(x) {

    clusters <- table(ref_group)

    base_max <- NA
    base_min <- NA

    for (cluster in names(clusters)) {
      if (base_metric == "mean")
        ref_metric = log2(mean(2^x[names(which(ref_group == cluster))] - 1) + 1)
      else if (base_metric == "median")
        ref_metric = log2(median(2^x[names(which(ref_group == cluster))] - 1) + 1)
      else
        stop("Unrecognized base statistic ", base_metric)

      if (is.na(base_max) || ref_metric > base_max)
        base_max <- ref_metric
      if (is.na(base_min) || ref_metric < base_min)
        base_min <- ref_metric
    }

    row_init <- rep(0, length(x))

    base_max <- base_max + scaling_factor
    base_min <- base_min - scaling_factor

    above_max <- which(x > base_max)
    below_min <- which(x < base_min)

    row_init[above_max] <- x[above_max] - base_max
    row_init[below_min] <- x[below_min] - base_min

    return(row_init)
  }

  res <- t(apply(m, 1, calc_score))

  colnames(res) <- colnames(m)

  return(res)
}

#' @importFrom grid grid.lines gpar
.draw_grid <- function(cnv_matrix, gene_positions_table) {

  gpc <- count_genes_per_chromosome(cnv_matrix, gene_positions_table)

  sgpc <- sum(gpc)

  rgpcs <- 0

  for (i in 1:(length(gpc) - 1)) {

    cgpc <- rgpcs + gpc[i]

    grid.lines(x = c(cgpc/sgpc, cgpc/sgpc), c(0, 1), gp = gpar(lty = 1, lwd = 1))

    rgpcs <- rgpcs + gpc[i]
  }

  grid.lines(x = c(0, 1), c(0, 0), gp = gpar(lty = 1, lwd = 1))
  grid.lines(x = c(0, 1), c(1, 1), gp = gpar(lty = 1, lwd = 1))
  grid.lines(x = c(0, 0), c(0, 1), gp = gpar(lty = 1, lwd = 1))
  grid.lines(x = c(1, 1), c(0, 1), gp = gpar(lty = 1, lwd = 1))
}

.classify_cells_simple <- function(data) {

  malignant <- setNames(rep("Unresolved", nrow(data)), nm = rownames(data))

  malignant[data$CNVDetected == "Detected"] <- "Malignant"
  malignant[data$CNVDetected == "Not detected"] <- "Nonmalignant"

  data$Malignant <- malignant

  return (data)
}

.classify_cells_cluster_correlation <- function(data, clusters, cnv_matrix, cell_class_init, min_cluster_cc_freq, verbose) {

  data$Cluster <- clusters[data$CellID]

  if (is.null(cell_class_init)) {

    if (isTRUE(verbose))
      message("Initializing cell Classification")

    data <- .init_cell_class(data)
  } else {

    if (isTRUE(verbose))
      message("Skipped cell classification initialization as it was supplied")

    data$CellClass <- cell_class_init
  }

  if (isTRUE(verbose))
    message("Classifying clusters")

  data <- .classify_clusters(data, min_cluster_cc_freq)

  if (isTRUE(verbose))
    message("Classifying NiMC (Normal in Malignant Cluster)")

  data <- .classify_nimc(data)

  if (isTRUE(verbose))
    message("Classifying MiNC (Malignant in Normal Cluster)")

  data <- .classify_minc(data)

  # Classify the low signal and low correlation scores
  if (!is.null(cnv_matrix)) {

    if (isTRUE(verbose))
      message("Classifying intermediate scores (low signal/correlation)")

    data <- .classify_intermediate_scores(data, cnv_matrix)
  }

  if (isTRUE(verbose))
    message("Classifying malignant cells")

  data <- .classify_cells(data)

  return (data)
}

.init_cell_class <- function(data) {

  cell_class <- setNames(rep("Unresolved", nrow(data)), data$CellID)

  malignant_cells <- data %>% filter(CNVDetected == "Detected") %>% pull(CellID)

  nonmalignant_cells <- data %>% filter(CNVDetected == "Not detected") %>% pull(CellID)

  cell_class[malignant_cells] <- "Malignant"
  cell_class[nonmalignant_cells] <- "Nonmalignant"

  data$CellClass <- cell_class

  return (data)
}

.classify_clusters <- function(data, min_cluster_cc_freq) {

  cluster_cell_class_freq <- data %>%
    group_by(Cluster, CellClass) %>%
    summarise (n = n()) %>%
    mutate(Freq = n / sum(n))

  cluster_class <- setNames(rep("Unresolved", length(unique(data$Cluster))), nm = unique(data$Cluster))

  malignant_clusters <- filter(cluster_cell_class_freq, CellClass == "Malignant", Freq >= min_cluster_cc_freq) %>% pull(Cluster)
  cluster_class[malignant_clusters] <- "Malignant"

  nonmalignant_clusters <- filter(cluster_cell_class_freq, CellClass == "Nonmalignant", Freq >= min_cluster_cc_freq) %>% pull(Cluster)
  cluster_class[nonmalignant_clusters] <- "Nonmalignant"

  cluster_class <- tibble(Cluster = names(cluster_class), ClusterClass = cluster_class)

  data <- data %>%
    left_join(cluster_class, by = "Cluster")

  return (data)
}

.classify_nimc <- function(data) {

  cell_class <- setNames(data$CellClass, data$CellID)

  nrm_in_mal <- data %>%
    filter(CellClass == "Nonmalignant" & ClusterClass == "Malignant")

  if (nrow(nrm_in_mal) > 0)
    cell_class[nrm_in_mal$CellID] <- "NiMC"

  data$CellClass <- cell_class

  return (data)
}

.classify_minc <- function(data) {

  cell_class <- setNames(data$CellClass, data$CellID)

  mal_in_nrm <- data %>%
    filter(CellClass == "Malignant" & ClusterClass == "Nonmalignant")

  if (nrow(mal_in_nrm) > 0)
    cell_class[mal_in_nrm$CellID] <- "MiNC"

  data$CellClass <- cell_class

  return (data)
}

.classify_intermediate_scores <- function(data, cnv_matrix) {

  cell_class <- setNames(data$CellClass, data$CellID)

  int_scores <- data %>%
    group_by(Cluster) %>%
    select(CellID, CNVDetected, Cluster, ClusterClass, CellClass) %>%
    filter(CNVDetected == "Low signal" | CNVDetected == "Low correlation", ClusterClass == "Malignant", CellClass == "Unresolved")

  nonmalignant_clusters <- unique(data %>% filter(ClusterClass == "Nonmalignant") %>% pull(Cluster))

  cnv_matrix <- t(cnv_matrix)

  int_scores$CNVCorOwn <- rep(0, nrow(int_scores))
  int_scores$CNVCorOth <- rep(0, nrow(int_scores))

  for (i in seq_len(nrow(int_scores))) {
    score_i <- int_scores[i, ]

    cnv_i <- cnv_matrix[, score_i$CellID]
    cnv_c <- cnv_matrix[, data %>% filter(Cluster == data$Cluster) %>% pull(CellID)]

    cnv_c <- rowMeans(cnv_c)

    int_scores$CNVCorOwn[i] <- cor(cnv_i, cnv_c)

    int_scores$CNVCorOth[i] <- max(sapply(nonmalignant_clusters, function(x) {
      cnv_x <- cnv_matrix[, data %>% filter(x == data$Cluster) %>% pull(CellID)]
      cnv_x <- rowMeans(cnv_x)
      cor(cnv_i, cnv_x)
    }))
  }

  cell_class[int_scores %>% filter(CNVCorOwn > 2*CNVCorOth) %>% pull(CellID)] <- "MbCC"
  cell_class[int_scores %>% filter(CNVCorOwn <= 2*CNVCorOth) %>% pull(CellID)] <- "Unresolved"

  data$CellClass <- cell_class

  return (data)
}

# Classify cells as Malignant\Nonmalignant\Unresolved
.classify_cells <- function(data) {

  malignant <- setNames(data$CellClass, data$CellID)

  malignant[malignant == "MbCC"] <- "Malignant"
  malignant[!(malignant %in% c("Malignant", "Nonmalignant"))] <- "Unresolved"

  data$Malignant <- malignant

  return (data)
}
