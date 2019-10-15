
### =========================================================================
### CNA matrix computation and plotting
### -------------------------------------------------------------------------
###

#'
#' @title CNA inference
#'
#' @description This function infers CNAs (chromosomal copy-number variations) from the
#' single-cell expression data. CNA inference is the main method of the scandal framework
#' for classifying malignant and non-malignant cells.
#'
#' @param object a \linkS4class{ScandalDataSet} object.
#' @param gene_positions_table a data frame containing all genes ordered by the position
#' of the gene on the chromosome and by the order of the chromosomes. The data frame should
#' contain the column names "Gene" and "CHR".
#' @param reference_cells a named vector of the cluster assignments of the reference cells.
#' The names should correspond to the cell IDs of the reference (non-malignant) cells. The
#' CNA matrix can be computed without a reference (with \code{reference=NULL}) but this is
#' not recommended as downstream comoutations using the inferred CNA matrix will be less
#' reliable.
#' @param max_genes maximal number of genes to use for computing the CNA matrix. Default
#' is 5000.
#' @param expression_limits a numeric vector with two elements representing the upper
#' and lower values with which to bound the centered expression matrix prior to
#' calculating the CNA matrix. This blunts the effect of noisy genes. Defaut is (-3, 3).
#' @param window number of genes to consider when calculating the running mean. Default
#' is a window of 100 genes.
#' @param scaling_factor a small constant by which to increase the calculated
#' (-BM, +BM) interval to compensate for possible noise. Default is 0.2.
#' @param initial_centering direction of centering the expression matrix (row-wise or
#' col-wise) prior to computing the CNA matrix. Accepts either strings "row" or "col",
#' default is "col".
#' @param base_metric a metric to use for calculating the (-BM, + BM) interval. Accepts
#' either strings "mean" or "median", default is "median".
#' @param verbose suppresses all messages from this function. Default is FALSE.
#'
#' @details The CNA algorithm is as follows:\cr
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
#' @return Returns the \linkS4class{ScandalDataSet} object with CNA matrix in the
#' "cna" element of the reducedDim slot (accessible by reducedDim(object, "cna")). Note
#' that the matrix is stored with cell IDs as row names and gene IDs as column names.
#'
#' @seealso The CNA inference method was defined and developed by **Dr. Itay Tirosh**
#' during his time at the *Broad Institute* and published in several high-impact papers
#' including the following paper from *Cell*:
#' https://www.cell.com/cell/fulltext/S0092-8674(17)31270-9.
#'
#' @author Avishay Spitzer
#'
#' @export
scandal_cna_infer <- function(object, gene_positions_table, reference_cells,
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

  cna_matrix <- scandal_cna_compute_matrix(x = x,
                                           gene_positions_table = gene_positions_table,
                                           reference_cells = reference_cells,
                                           cells_subset = NULL,
                                           max_genes = max_genes,
                                           expression_limits = expression_limits,
                                           window = window,
                                           scaling_factor = scaling_factor,
                                           initial_centering = initial_centering,
                                           base_metric = base_metric,
                                           verbose = verbose)

  if (is.null(cna_matrix))
    stop("An error has occured while computing the CNA matrix")

  reducedDim(object, "cna") <- t(cna_matrix)

  return (object)
}

#'
#' @title CNA matrix computation
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
scandal_cna_compute_matrix <- function(x, gene_positions_table, reference_cells, cells_subset = NULL,
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

  cna_matrix <- .scandal_compute_cna_matrix(x = x,
                                            gene_pos_tbl = gene_positions_table,
                                            reference_cells = reference_cells,
                                            cells_subset = cells_subset,
                                            max_genes = max_genes,
                                            expression_limits = expression_limits,
                                            window = window,
                                            scaling_factor = scaling_factor,
                                            initial_centering = initial_centering,
                                            base_metric = base_metric,
                                            verbose = verbose)

  if (is.null(cna_matrix))
    stop("An error has occured while computing the CNA matrix")

  return (cna_matrix)
}

#' @importFrom ComplexHeatmap Heatmap rowAnnotation HeatmapAnnotation draw decorate_heatmap_body
#' @importFrom circlize colorRamp2
#' @importFrom grDevices dev.off jpeg
#' @importFrom infercna splitGenes
#' @export
scandal_cna_plot <- function(object, cna_matrix = NULL, reference_cells = NULL, show_reference = FALSE,
                             row_annotation = NULL, row_annotation_cols = NULL, low = "dodgerblue", mid = "white", high = "red",
                             save_to_file = TRUE, show_plot = TRUE, filename = NULL, dirname = ".", width = 1600, height = 1000, quality = 100, verbose = FALSE,
                             cluster_rows = FALSE, ...) {

  if (!is.null(object))
    cna_matrix <- reducedDim(object, "cna")

  if (is.null(cna_matrix))
    stop("CNA matrix not found, did you forget to run scandal_infer_cna first?")

  if (isTRUE(verbose))
    message(sprintf("Plotting CNA: data range (%.2f, %.2f), limit set to (-1, 1)", min(cna_matrix), max(cna_matrix)))

  cell_ids <- NULL
  if (!is.null(reference_cells)) {

    malignant_cells <- rownames(cna_matrix)[!rownames(cna_matrix) %in% reference_cells]

    if (isTRUE(show_reference))
      cell_ids <- c(malignant_cells, reference_cells)
    else
      cell_ids <- malignant_cells
  } else
    cell_ids <- rownames(cna_matrix)

  plot_name <- paste0("cna", sample(1:10^9, 1))

  ra <- NULL
  if (!is.null(row_annotation)) {
    rann <- as.data.frame(row_annotation)[cell_ids, ]

    if (is.null(row_annotation_cols))
      ra <- rowAnnotation(df = rann, show_annotation_name = FALSE)
    else
      ra <- rowAnnotation(df = rann, col = row_annotation_cols, show_annotation_name = FALSE)
  }

  h <- Heatmap(cna_matrix[cell_ids, ], name = plot_name,
               col = colorRamp2(c(-1, 0, 1), c(low, mid, high)),
               show_row_names = FALSE, show_column_names = FALSE, cluster_columns = FALSE, cluster_rows = cluster_rows,
               left_annotation = ra,
               heatmap_legend_param = list(at = c(-1, -0.5, 0, 0.5, 1), title = "Inferred CNA\n(log2 ratio)"), ...)

  if (isTRUE(show_plot)) {
    draw(h)

    decorate_heatmap_body(plot_name, .draw_grid(cna_matrix))
  }

  if (isTRUE(save_to_file)) {

    if (is.null(filename))
      filename <- paste0(plot_name, ".jpeg")

    file <- paste0(dirname, "/", filename)

    jpeg(filename = file, width = width, height = height, quality = quality)
    draw(h)
    decorate_heatmap_body(plot_name, .draw_grid(cna_matrix))
    dev.off()
  }

  invisible(h)
}

### =========================================================================
### CNA scores computation, classification and plotting
### -------------------------------------------------------------------------
###

#'
#' @title Compute the cell CNA score
#'
#' @description Computes the CNA score for each cell which is composed of the CNA
#' signal (the sum of squares of computed CNAs) and CNA correlation (the pearson's
#' correlation between the cell CNA profile and the tumor's CNA profile).
#'
#' @param object
#' @param cna_matrix
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
scandal_cna_compute_scores <- function(object, cna_matrix = NULL, hotspot_threshold = .9, verbose = FALSE) {

  stopifnot(!(is.null(object) & is.null(cna_matrix)))
  stopifnot(is.null(object) | is_scandal_object(object))

  if (!is.null(object))
    cna_matrix <- reducedDim(object, "cna")

  if (is.null(cna_matrix))
    stop("CNA matrix not found, did you forget to run scandal_infer_cna first?")

  cna_matrix <- t(cna_matrix)

  # Square the CNA matrix
  tmp <- cna_matrix^2

  # Calculate the hotspots - areas with high signal
  pos_mean <- rowMeans(tmp)
  hotspots <- names(pos_mean[pos_mean >= quantile(pos_mean, hotspot_threshold)])
  tmp <- tmp[hotspots, ]

  # Compute the mean squared CNA signal per
  cna_signal <- colMeans(tmp)

  # Compute the tumor CNA profile
  tumor_cna0 <- rowMeans(cna_matrix)

  # Compute the correlation between each cell and the tumor CNA profile
  cna_correlation <- apply(cna_matrix, 2, cor, y = tumor_cna0, method = "pearson")

  # Return the result as a data frame with two variables
  res <- tibble(CellID = colnames(cna_matrix), Signal = cna_signal, Correlation = cna_correlation)

  return (res)
}

#'
#' @title Classify CNA scores
#'
#' @description This function classifies the CNA scores according to the supplied signal and correlation thresholds.
#'
#' @param cna_scores
#' @param signal_threshold
#' @param correlation_threshold
#' @param verbose
#'
#' @return
#'
#' @seealso \link{scandal_cna_compute_scores}
#'
#' @author Avishay Spitzer
#'
#' @export
scandal_cna_classify_scores <- function(cna_scores, signal_threshold = 0.05, correlation_threshold = 0.05, verbose = FALSE) {

  cna_detected <- rep("Non-classifiable", nrow(cna_scores))

  cna_detected[which(cna_scores$Signal >= signal_threshold & cna_scores$Correlation >= correlation_threshold)] <- "Detected"
  cna_detected[which(cna_scores$Signal <  signal_threshold & cna_scores$Correlation <  correlation_threshold)] <- "Not detected"
  cna_detected[which(cna_scores$Signal <  signal_threshold & cna_scores$Correlation >= correlation_threshold)] <- "Low signal"
  cna_detected[which(cna_scores$Signal >= signal_threshold & cna_scores$Correlation <  correlation_threshold)] <- "Low correlation"

  cna_scores$CNADetected <- cna_detected

  return (cna_scores)
}

#' @importFrom dplyr %>% group_by summarise n mutate filter select pull left_join
#' @importFrom tibble tibble
#' @importFrom stats setNames
#' @export
scandal_cna_classify_cells <- function(cna_scores, clusters, cna_matrix = NULL, cell_class_init = NULL, min_cluster_cc_freq = .5, return_all = FALSE, verbose = FALSE) {

  # Simple classification - leaving the low signal and low correlation cells as "Unresolved"
  if (is.null(clusters))
    data <- .classify_cells_simple(data = cna_scores)
  else # Complex classifiaction using cluster correlation
    data <- .classify_cells_cluster_correlation(data = cna_scores,
                                                clusters = clusters,
                                                cna_matrix = cna_matrix,
                                                cell_class_init = cell_class_init,
                                                min_cluster_cc_freq = min_cluster_cc_freq,
                                                verbose = verbose)

  if (isTRUE(return_all))
    return (data)
  else
    return (data$Malignant)
}

#' @export
scandal_cna_get_malignant_cells <- function(cna_scores) {

  res <- (cna_scores %>%
            select(CellID, Malignant) %>%
            filter(Malignant == "Malignant"))$CellID

  return (res)
}

#' @export
scandal_cna_scores_plot <- function(cna_scores, signal_threshold = 0.05, correlation_threshold = 0.5, title = NULL, labels = NULL, name = NULL, verbose = FALSE) {

  #cna_scores <- scandal_cna_classify_scores(cna_scores, signal_threshold = signal_threshold, correlation_threshold = correlation_threshold)

  #cna_detected <- cna_scores$CNADetected

  #non_classifiable <- (length(which(cna_detected == "Low signal" | cna_detected == "Low correlation")) / length(cna_detected)) * 100

  p <- scandal_scatter_plot(x = cna_scores$Correlation, y = cna_scores$Signal, labels = labels, color_legend_name = name,
                            title = title, xlab = "CNA Correlation", ylab = "CNA Signal", plot_ordered = FALSE) +
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
count_genes_per_chromosome <- function(cna_matrix, gene_pos_tbl) {

  chr_pos_tbl <- table(gene_pos_tbl[gene_pos_tbl$Gene %in% colnames(cna_matrix), "CHR"])[CHRs]

  return (chr_pos_tbl)
}

### =========================================================================
### Internal functions
### -------------------------------------------------------------------------
###

#CHRs <- c("chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8", "chr9", "chr10", "chr11", "chr12", "chr13", "chr14", "chr15", "chr16", "chr17", "chr18", "chr19", "chr20", "chr21", "chr22", "chrX", "chrY")
CHRs <- c("chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8", "chr9", "chr10", "chr11", "chr12", "chr13", "chr14", "chr15", "chr16", "chr17", "chr18", "chr19", "chr20", "chr21", "chr22", "chrX")

.scandal_compute_cna_matrix <- function(x, gene_pos_tbl, reference_cells, cells_subset, max_genes, expression_limits, window, scaling_factor, initial_centering, base_metric, verbose) {

  if (!is.null(cells_subset))
    x <- x[, unique(c(names(reference_cells), cells_subset))]

  cna_matrix <- .prepare_cna_matrix(x,
                                    gene_pos_tbl = gene_pos_tbl,
                                    max_genes = max_genes,
                                    expression_limits = expression_limits,
                                    initial_centering = initial_centering,
                                    verbose = verbose)

  gene_pos_tbl <- gene_pos_tbl[gene_pos_tbl$Gene %in% rownames(cna_matrix), ]

  if (!is.null(reference_cells)) {
    m_ref <- cna_matrix[, names(reference_cells)]
    m_mal <- cna_matrix[, !(colnames(cna_matrix) %in% names(reference_cells))]
  } else
    m_mal <- cna_matrix

  cna_ref <- NULL
  cna_mal <- NULL
  genes_per_chr <- replicate(length(CHRs), 0)
  names(genes_per_chr) <- CHRs

  for (c in CHRs) {

    if (isTRUE(verbose))
      message(paste0("Calculate CNA matrix for ", c))

    chr_pos_tbl <- gene_pos_tbl[as.character(gene_pos_tbl$CHR) == c, ]

    if (!is.null(reference_cells))
      m1 <- m_ref[rownames(m_ref) %in% chr_pos_tbl$Gene, ]

    m2 <- m_mal[rownames(m_mal) %in% chr_pos_tbl$Gene, ]

    if (!is.null(reference_cells))
      m1 <- .calc_cna_matrix(m1, window = window, verbose = verbose)

    m2 <- .calc_cna_matrix(m2, window = window, verbose = verbose)

    if (!is.null(reference_cells))
      cna_ref <- rbind(cna_ref, m1)

    cna_mal <- rbind(cna_mal, m2)

    genes_per_chr[c] <- nrow(m2)
  }

  if (!is.null(reference_cells))
    cna_matrix <- cbind(cna_ref, cna_mal)
  else
    cna_matrix <- cna_mal

  if (isTRUE(verbose))
    message(paste0("Median-centering CNA matrix cellwise"))

  cna_matrix <- center_matrix(cna_matrix, by = "col", method = "median", scale = FALSE, verbose = verbose)

  if (!is.null(reference_cells)) {

    if (isTRUE(verbose))
      message(paste0("Calculating CNA score"))

    cna_matrix <- .assign_cna_score(cna_matrix, ref_group = reference_cells, scaling_factor = scaling_factor, base_metric = base_metric, verbose = verbose)
  }

  cna_matrix <- cna_matrix[, colnames(x)]

  if (isTRUE(verbose))
    message(paste0("CNA matrix created!"))

  return (cna_matrix)
}

.prepare_cna_matrix <- function(cna_matrix, gene_pos_tbl, max_genes, expression_limits, initial_centering, verbose) {

  if (isTRUE(verbose))
    message("Preprocessing CNA matrix...")

  if (isTRUE(verbose))
    message(paste0("Expression matrix contains ", ncol(cna_matrix), " high quality cells"))

  mean_exp <- log2(rowMeans(cna_matrix) + 1)

  mean_exp <- mean_exp[order(mean_exp, decreasing = TRUE)]
  passQC <- head(mean_exp, n = min(max_genes, length(mean_exp)))

  cna_matrix <- cna_matrix[names(passQC), ]

  if (isTRUE(verbose))
    message(paste0("Keeping ", min(max_genes, length(mean_exp)), " highest expressed genes"))

  gene_pos_tbl <- gene_pos_tbl[gene_pos_tbl$Gene %in% rownames(cna_matrix), ]

  # Reorder the expression matrix according to the positions of the genes in each chromosome
  if (isTRUE(verbose))
    message(paste0("Reordering genes according to chromosomal location"))

  cna_matrix <- cna_matrix[gene_pos_tbl$Gene, ]

  # Log-transform the matrix
  if (isTRUE(verbose))
    message(paste0("Log2 expression matrix [TPM/10]"))

  cna_matrix <- log_transform(cna_matrix)

  # Center the matrix
  if (isTRUE(verbose))
    message(sprintf("Mean-centering expression matrix %s-wise", initial_centering))

  cna_matrix <- center_matrix(cna_matrix, by = initial_centering, method = "mean", scale = FALSE, verbose = verbose)

  # Bound the expression matrix to reduce the effect of specific highly expressed genes on the mean expression
  if (!is.null(expression_limits)) {

    if (isTRUE(verbose))
      message(paste0("Expresion matrix is in the range (", min(cna_matrix), ", ", max(cna_matrix), ")"))

    if (isTRUE(verbose))
      message(paste0("Limiting matrix to values betwwen (", expression_limits[1], ", ", expression_limits[2], ")"))

    cna_matrix <-.bound_matrix(cna_matrix, expression_limits = expression_limits, verbose = verbose)

    if (max(cna_matrix) > expression_limits[2] | min(cna_matrix) < expression_limits[1])
      stop("Matrix limiting failed")
  }

  if (isTRUE(verbose))
    message(paste0("Preprocessing done!"))

  return (cna_matrix)
}

.bound_matrix <- function(cna_matrix, expression_limits, verbose) {

  cna_matrix[cna_matrix > expression_limits[2]] <- expression_limits[2]
  cna_matrix[cna_matrix < expression_limits[1]] <- expression_limits[1]

  return (cna_matrix)
}

#' @importFrom caTools runmean
.calc_cna_matrix = function(m, window, verbose) {

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
.assign_cna_score <- function(m, ref_group, scaling_factor, base_metric, verbose) {

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
.draw_grid <- function(cna_matrix) {

  gpc <- infercna::splitGenes(colnames(cna_matrix)) %>% lengths

  #gpc <- count_genes_per_chromosome(cna_matrix, gene_positions_table)

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

  malignant[data$CNADetected == "Detected"] <- "Malignant"
  malignant[data$CNADetected == "Not detected"] <- "Nonmalignant"

  data$Malignant <- malignant

  return (data)
}

.classify_cells_cluster_correlation <- function(data, clusters, cna_matrix, cell_class_init, min_cluster_cc_freq, verbose) {

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
  if (!is.null(cna_matrix)) {

    if (isTRUE(verbose))
      message("Classifying intermediate scores (low signal/correlation)")

    data <- .classify_intermediate_scores(data, cna_matrix)
  }

  if (isTRUE(verbose))
    message("Classifying malignant cells")

  data <- .classify_cells(data)

  return (data)
}

.init_cell_class <- function(data) {

  cell_class <- setNames(rep("Unresolved", nrow(data)), data$CellID)

  malignant_cells <- data %>% filter(CNADetected == "Detected") %>% pull(CellID)

  nonmalignant_cells <- data %>% filter(CNADetected == "Not detected") %>% pull(CellID)

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

.classify_intermediate_scores <- function(data, cna_matrix) {

  cell_class <- setNames(data$CellClass, data$CellID)

  int_scores <- data %>%
    group_by(Cluster) %>%
    select(CellID, CNADetected, Cluster, ClusterClass, CellClass) %>%
    filter(CNADetected == "Low signal" | CNADetected == "Low correlation", CellClass == "Unresolved")

  nonmalignant_clusters <- unique(data %>% filter(ClusterClass == "Nonmalignant") %>% pull(Cluster))
  malignant_clusters <- unique(data %>% filter(ClusterClass == "Malignant") %>% pull(Cluster))

  cna_matrix <- t(cna_matrix)

  int_scores$CNACorOwn <- rep(0, nrow(int_scores))
  int_scores$CNACorOth <- rep(0, nrow(int_scores))

  for (i in seq_len(nrow(int_scores))) {
    score_i <- int_scores[i, ]

    cna_i <- cna_matrix[, score_i$CellID]
    cna_c <- cna_matrix[, data %>% filter(Cluster == score_i$Cluster) %>% pull(CellID)]

    cna_c <- rowMeans(cna_c)

    int_scores$CNACorOwn[i] <- cor(cna_i, cna_c)

    if (score_i$ClusterClass == "Malignant")
      cc <- nonmalignant_clusters
    else if (score_i$ClusterClass == "Nonmalignant")
      cc <- malignant_clusters
    else
      stop("Unclassified cluster")

    int_scores$CNACorOth[i] <- max(sapply(cc, function(x) {
      cna_x <- cna_matrix[, data %>% filter(x == data$Cluster) %>% pull(CellID)]
      cna_x <- rowMeans(cna_x)
      cor(cna_i, cna_x)
    }))
  }

  cell_class[int_scores %>% filter(CNACorOwn > 2*CNACorOth, ClusterClass == "Malignant") %>% pull(CellID)] <- "MbCC"
  cell_class[int_scores %>% filter(CNACorOwn > 2*CNACorOth, ClusterClass == "Nonmalignant") %>% pull(CellID)] <- "NbCC"
  cell_class[int_scores %>% filter(CNACorOwn <= 2*CNACorOth) %>% pull(CellID)] <- "Unresolved"

  data$CellClass <- cell_class

  return (data)
}

# Classify cells as Malignant\Nonmalignant\Unresolved
.classify_cells <- function(data) {

  malignant <- setNames(data$CellClass, data$CellID)

  malignant[malignant == "MbCC"] <- "Malignant"
  malignant[malignant == "NbCC"] <- "Nonmalignant"
  malignant[!(malignant %in% c("Malignant", "Nonmalignant"))] <- "Unresolved"

  data$Malignant <- malignant

  return (data)
}
