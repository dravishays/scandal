
### =========================================================================
### Exported functions
### -------------------------------------------------------------------------
###

#' @export
scandal_infer_cnv <- function(object, gene_positions_table, reference_cells,
                              max_genes = 5000, expression_limits = c(-3, 3), window = 100, scaling_factor = 0.2,
                              initial_centering = "col", base_metric = "median",
                              verbose = FALSE) {
  stopifnot(is_scandal_object(object))

  x <- as.matrix(unprocessedData(object))[, colnames(object)]

  cnv_matrix <- scandal_compute_cnv_matrix(x = x,
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

#' @export
scandal_compute_cnv_matrix <- function(x, gene_positions_table, reference_cells,
                                       max_genes = 5000, expression_limits = c(-3, 3), window = 100, scaling_factor = 0.2,
                                       initial_centering = "col", base_metric = "median",
                                       verbose = FALSE) {

  stopifnot(is_valid_assay(x))
  stopifnot(!is.null(gene_positions_table), is.data.frame(gene_positions_table), !is.null(colnames(gene_positions_table)), c("Gene", "CHR") %in% colnames(gene_positions_table))
  stopifnot(!is.null(reference_cells), is.vector(reference_cells), is.numeric(reference_cells) | is.character(reference_cells), !is.null(names(reference_cells)))
  stopifnot(is.numeric(max_genes), max_genes > 0)
  stopifnot(!is.null(expression_limits), is.vector(expression_limits), is.numeric(expression_limits), length(expression_limits) == 2)
  stopifnot(is.numeric(window), window > 0)
  stopifnot(is.numeric(scaling_factor), scaling_factor > 0)
  stopifnot(initial_centering %in% c("col", "row"))
  stopifnot(base_metric %in% c("mean", "median"))

  x <- as.matrix(x)

  cnv_matrix <- .infer_cnv(x = x,
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

#' @importFrom ComplexHeatmap Heatmap rowAnnotation HeatmapAnnotation
#' @importFrom circlize colorRamp2
#' @export
scandal_cnv_plot <- function(object, cnv_matrix = NULL, gene_positions_table, reference_cells = NULL, show_reference = FALSE,
                             row_annotation = NULL, row_annotation_cols = NULL,
                             low = "dodgerblue", mid = "white", high = "red", verbose = FALSE,
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

  h <- Heatmap(cnv_matrix[cell_ids, ], name = plot_name,
               col = colorRamp2(c(-1, 0, 1), c(low, mid, high)),
               show_row_names = FALSE, show_column_names = FALSE, cluster_columns = FALSE, cluster_rows = cluster_rows,
               heatmap_legend_param = list(at = c(-1, -0.5, 0, 0.5, 1), title = "Inferred CNV\n(log2 ratio)"), ...)

  if (!is.null(row_annotation)) {
    rann <- as.data.frame(row_annotation)[cell_ids, ]

    ra <- rowAnnotation(df = rann, col = row_annotation_cols)

    p <- ra + h
  } else
    p <- h

  draw(p)

  decorate_heatmap_body(plot_name, {

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
  })

  invisible(p)
}

### =========================================================================
### Internal functions
### -------------------------------------------------------------------------
###

#CHRs <- c("chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8", "chr9", "chr10", "chr11", "chr12", "chr13", "chr14", "chr15", "chr16", "chr17", "chr18", "chr19", "chr20", "chr21", "chr22", "chrX", "chrY")
CHRs <- c("chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8", "chr9", "chr10", "chr11", "chr12", "chr13", "chr14", "chr15", "chr16", "chr17", "chr18", "chr19", "chr20", "chr21", "chr22", "chrX")

.infer_cnv <- function(x, gene_pos_tbl, reference_cells, max_genes, expression_limits, window, scaling_factor, initial_centering, base_metric, verbose) {

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
  message(paste0("Reordering genes according to chromosomal location"))
  cnv_matrix <- cnv_matrix[gene_pos_tbl$Gene, ]

  # Log-transform the matrix
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
