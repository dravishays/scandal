
### =========================================================================
### Exported algorithm functions (start)
### -------------------------------------------------------------------------
###

#'
#' @title Programs of intra-sample heterogeneity
#'
#' @description This function extracts transpriptomic programs
#'
#' @param object
#' @param samples
#' @param clustering_data
#' @param algorithm
#' @param rank
#' @param ngenes1
#' @param ngenes2
#' @param sd_threshold
#' @param filter_method
#' @param bin_control
#' @param n_control_bins
#' @param n_bin_genes
#' @param return_all
#' @param verbose
#' @param ... further arguments passed to nmf function
#'
#' @details Generally, the algorithm extracts programs in a bottoms-up approach
#' starting from extracting programs within each individual sample. The algorithm
#' then detecs highly variable programs that best represent coherent program clusters
#' and aggregates these clusters into programs that generalize accross all samples.
#' \cr
#' The algorithm performs the following steps:
#' \cr
#' \enumerate{
#'   \item Call \link{prepare_samples} to prepare a \linkS4class{ScandalDataSet} for
#'   each individual sample. This step can be bypassed by supplying a valid \code{samples}
#'   argument.
#'   \item For each individual sample generate a predefined number of clusters
#'   (configurable by the \code{rank} parameter). At the moment the only supported
#'   clustering algorithm is NMF via \link{nmf_run}. This step can be bypassed by supplying
#'   a valid \code{clustering_data} argument.
#'   \item For each sample extract an initial list of within-sample programs by calling
#'   \link{nmf_extract_programs}. The number of genes in each program is configurable by
#'   the \code{ngenes1} parameter. The number of initial within-sample programs correpsonds
#'   to \code{rank}.
#'   \item For each sample score the cells for each of the within-sample programs by calling
#'   \link{score_within_samples}. See \link{scrabble::score} for more details about how to
#'   score cells while controling for differences in cell complexities.
#'   \item For each sample compute the standard deviation of scores for each within-sample
#'   program by calling \link{compute_programs_sd}.
#'   \item
#' }
#' \cr
#'
#' @return
#'
#' @seealso
#'
#' @author Avishay Spitzer
#'
#' @export
scandal_programs_of_intra_sample_heterogeneity <- function(object, samples = NULL, clustering_data = NULL,
                                                           algorithm = "nmf", rank = 10, ngenes1 = 50, ngenes2 = 30, sd_threshold = 0.8, filter_method = "relative",
                                                           bin_control = TRUE, n_control_bins = 25, n_bin_genes = 100, return_all = FALSE, verbose = FALSE, ...) {

  if (is.null(samples))
    samples <- prepare_samples(object = object, verbose = verbose)
  else {
    # Do sanity checks

    if (isTRUE(verbose))
      message("Sample data supplied, skipping prepare_samples")
  }

  if (algorithm == "nmf") {

    if (is.null(clustering_data))
      ws_clustering_data <- nmf_run(samples = samples, rank = rank, ..., verbose = verbose)
    else {
      # Do sanity checks

      if (isTRUE(verbose))
        message("NMF clustering data supplied, skipping nmf_run")

      ws_clustering_data <- clustering_data
    }

    ws_programs <- nmf_extract_programs(nmf_data = ws_clustering_data, n = ngenes1, verbose = verbose)
  }
  else
    stop("Unsupported algorithm ", algorithm)

  ws_scores <- score_within_samples(samples = samples, ws_programs = ws_programs, bin_control = bin_control, n_control_bins = n_control_bins, n_bin_genes = n_bin_genes, verbose = verbose)

  ws_score_sd <- compute_programs_sd(ws_scores = ws_scores, verbose = verbose)

  # Plot here SD distribution

  variable_programs <- within_sample_variable_programs(ws_score_sd = ws_score_sd, ws_programs = ws_programs, sd_threshold = sd_threshold, filter_method = filter_method, verbose = verbose)

  variable_programs <- unlist(variable_programs, recursive = FALSE)

  bs_scores <- score_between_samples(x = assay(object), programs = variable_programs,
                                     bin_control = bin_control, n_control_bins = n_control_bins, n_bin_genes = n_bin_genes, verbose = verbose)

  program_clusters <- cluster_variable_programs(bs_scores = bs_scores, rank = rank, verbose = verbose)

  # Plot here program correlation

  if (algorithm == "nmf") {
    gene_scores <- nmf_score_genes(nmf_data = ws_clustering_data, variable_programs = variable_programs, program_clusters = program_clusters, verbose = verbose)
  } else
    stop("Unsupported algorithm ", algorithm)

  metaprograms <- metaprograms(gene_scores = gene_scores, n_sig_genes = ngenes2, verbose = verbose)

  metaprogram_scores <- score_between_samples(x = assay(object), programs = metaprograms,
                                              bin_control = bin_control, n_control_bins = n_control_bins, n_bin_genes = n_bin_genes, verbose = verbose)

  assigned_metaprogram <- assign_metaprograms(metaprogram_scores = metaprogram_scores, verbose = verbose)

  scandal_results <- list()
  scandal_results$samples <- samples
  scandal_results$wsClusteringData <- ws_clustering_data
  scandal_results$wsPrograms <- ws_programs
  scandal_results$wsScores <- ws_scores
  scandal_results$wsScoreSDs <- ws_score_sd
  scandal_results$bsScores <- bs_scores
  scandal_results$variablePrograms <- variable_programs
  scandal_results$programClusters <- program_clusters
  scandal_results$geneScores <- gene_scores
  scandal_results$metaPrograms <- metaprograms
  scandal_results$mpScores <- metaprogram_scores
  scandal_results$assignedMetaprogram <- assigned_metaprogram

  return (scandal_results)
}

#' @author Avishay Spitzer
#'
#' @importFrom stats setNames
#'
#' @export
prepare_samples <- function(object, verbose = FALSE) {

  if (isTRUE(verbose))
    message("Generating individual samples")

  samples <- setNames(lapply(sampleIDs(object), function(sname) {
    sdata <- scandal_inspect_samples(object, sample_ids = sname, node_id = sname, verbose = verbose)
    sdata <- scandal_setup_analysis(sdata, verbose = verbose)
    return (sdata)
  }), sampleIDs(object))

  return (samples)
}

#' @author Avishay Spitzer
#'
#' @importFrom NMF nmf
#'
#' @export
nmf_run <- function(samples, rank = 10, verbose = FALSE, ...) {

  res <- list()

  run_time <- 0

  for (sname in names(samples)) {

    sdata <- samples[[sname]]

    m <- .nmf_matrix(sdata)

    if (isTRUE(verbose))
      message("Running NMF for ", sname)

    nmf_res <- nmf(x = m, rank = rank, ...)

    if (isTRUE(verbose))
      message("Done - ", nmf_res@runtime[3], "s")

    run_time <- run_time + nmf_res@runtime[3]

    res[[sname]] <- nmf_res
  }

  if (isTRUE(verbose))
    message("NMF done - overall runtime ", run_time, "s")

  return (res)
}

#' @author Avishay Spitzer
#'
#' @importFrom stats setNames
#'
#' @export
nmf_extract_programs <- function(nmf_data, n = 50, verbose = FALSE) {

  if (isTRUE(verbose))
    message("Extracting within-sample programs from NMF result")

  programs <- setNames(lapply(names(nmf_data),
                              function(sname) .extract_programs(.W(nmf_data[[sname]], sname), n = n)),
                         names(nmf_data))

  return (programs)
}

#' @author Avishay Spitzer
#'
#' @importFrom scrabble score
#'
#' @export
score_within_samples <- function(samples, ws_programs, bin_control = TRUE, n_control_bins = 25, n_bin_genes = 100, ..., verbose = FALSE) {

  if (isTRUE(verbose))
    message("Scoring programs within each sample")

  res <- list()

  for (sname in names(samples)) {

    sdata <- samples[[sname]]

    programs <- ws_programs[[sname]]

    m <- center_matrix(as.matrix(assay(sdata)), by = "row", method = "mean", scale = FALSE)

    scores <- score(mat = m, groups = programs, center = FALSE, bin.control = bin_control, nbin = n_control_bins, n = n_bin_genes, binmat = as.matrix(assay(sdata)), ...)

    res[[sname]] <- scores
  }

  return (res)
}

#' @author Avishay Spitzer
#'
#' @importFrom matrixStats colSds
#'
#' @export
compute_programs_sd <- function(ws_scores, verbose = FALSE) {

  if (isTRUE(verbose))
    message("Computing programr standard deviation (within each sample)")

  res <- list()

  for (sname in names(ws_scores)) {

    scores <- ws_scores[[sname]]

    res[[sname]] <- setNames(colSds(scores), colnames(scores))
  }

  return (res)
}

#' @author Avishay Spitzer
#'
#' @export
within_sample_variable_programs <- function(ws_score_sd, ws_programs, sd_threshold, filter_method = "relative", verbose = FALSE) {

  if (isTRUE(verbose))
    message("Filtering out invariant programs (within each sample")

  if (filter_method == "relative") {
    threshold <- quantile(unlist(ws_score_sd, recursive = FALSE), sd_threshold)
  } else
    threshold <- sd_threshold

  ws_score_sd <- lapply(ws_score_sd, function(x) x[x >= threshold])

  ws_programs <- setNames(lapply(names(ws_programs), function(sname) ws_programs[[sname]][names(ws_score_sd[[sname]])]), names(ws_programs))

  return (ws_programs)
}

#' @author Avishay Spitzer
#'
#' @importFrom scrabble score
#'
#' @export
score_between_samples <- function(x, programs, bin_control = TRUE, n_control_bins = 25, n_bin_genes = 100, ..., verbose = FALSE) {

  if (isTRUE(verbose))
    message("Scoring programs between samples")

  x <- as.matrix(x)

  centered_data <- center_matrix(x, by = "row", method = "mean", scale = FALSE)

  scores <- score(mat = centered_data, groups = programs, center = FALSE, bin.control = bin_control, nbin = n_control_bins, n = n_bin_genes, binmat = x)

  return (scores)
}

#' @author Avishay Spitzer
#'
#' @importFrom ConsensusClusterPlus ConsensusClusterPlus calcICL
#'
#' @export
cluster_variable_programs <- function(bs_scores, rank = 10, verbose = FALSE) {

  if (isTRUE(verbose))
    message("Computing program clusters")

  ccp <- ConsensusClusterPlus(bs_scores, maxK = rank, pItem = 1, pFeature = 1, clusterAlg = "hc", distance = "pearson")
  icl <- calcICL(ccp)

  sum_ic <- setNames(rep(0, length(2:rank)), as.character(2:rank))
  for (k in 2:rank) {
    k_cc <- icl$itemConsensus[icl$itemConsensus[, "k"] == k, ]
    sum_ic[as.character(k)] <- sum(k_cc[, "itemConsensus"], na.rm = TRUE)
  }

  max_sum <- max(sum_ic)

  best_cc <- as.numeric(names(sum_ic[sum_ic == max_sum]))
  best_cc <- best_cc[length(best_cc)]

  if (isTRUE(verbose))
    message(sprintf("Detected %d consensus clusters of programs", best_cc))

  return (ccp[[best_cc]])
}

#' @author Avishay Spitzer
#'
#' @export
nmf_score_genes <- function(nmf_data, variable_programs, program_clusters, verbose = FALSE) {

  if (isTRUE(verbose))
    message("Scoring genes within each program cluster")

  cc <- program_clusters$consensusClass

  genes <- list()
  for (i in unique(cc)) {
    pnames <- names(cc[which(cc == i)])
    genes[[paste0("P", i)]] <- unique(unlist(variable_programs[pnames], recursive = FALSE))
  }

  vp_nmf_data <- matrix(0, nrow = length(unique(unlist(genes, recursive = FALSE))), ncol = length(variable_programs),
                        dimnames = list(unique(unlist(genes, recursive = FALSE)), names(variable_programs)))

  for (cname in colnames(vp_nmf_data)) {
    sample <- gsub("\\..*", "", cname)
    pname <- gsub(".*\\.", "", cname)

    p <- variable_programs[[cname]]

    W <- .W(nmf_data[[sample]])[p, pname]

    vp_nmf_data[names(W), cname] <- W
  }

  vp_nmf_data <- log2(vp_nmf_data + 1)

  nmf_scores <- matrix(0, nrow = nrow(vp_nmf_data), ncol = length(unique(cc)),
                       dimnames = list(rownames(vp_nmf_data), paste0("MP", unique(cc))))

  for (c in unique(cc)) {
    pnames <- names(which(cc == c))
    mpname <- paste0("MP", c)

    if (length(pnames) > 1)
      nmf_scores[, mpname] <- rowMeans(vp_nmf_data[, pnames])
    else
      nmf_scores[, mpname] <- vp_nmf_data[, pnames]
  }

  return (nmf_scores)
}

#' @author Avishay Spitzer
#'
#' @export
metaprograms <- function(gene_scores, n_sig_genes = 30, verbose = FALSE) {

  if (isTRUE(verbose))
    message("Computing meta-programs")

  metaprograms <- list()
  for (mpname in colnames(gene_scores)) {
    metaprograms[[mpname]] <- head(names(sort(gene_scores[, mpname], decreasing = TRUE)), n = n_sig_genes)
  }

  return (as.data.frame(metaprograms, stringsAsFactors = FALSE))
}

#' @author Avishay Spitzer
#'
#' @export
assign_metaprograms <- function(metaprogram_scores, verbose = FALSE) {

  max_scores <- apply(metaprogram_scores, 1, which.max)
  max_scores <- paste0("MP", max_scores)
  names(max_scores) <- rownames(metaprogram_scores)

  return(max_scores)
}

### -------------------------------------------------------------------------
### Exported algorithm functions (end)
### =========================================================================
###

### =========================================================================
### Exported plotting functions (start)
### -------------------------------------------------------------------------
###

#' @author Avishay Spitzer
#'
#' @importFrom reshape2 melt
#'
#' @export
scandal_programs_sd_plot <- function(ws_score_sd, sd_threshold = 0.5, title = NULL) {

  df <- as.data.frame(ws_score_sd)

  p <- ggplot(melt(df), aes(x = variable, y = value, fill = variable)) +
        geom_boxplot() +
        geom_point(size = 2) +
        geom_hline(yintercept = sd_threshold, size = 1, color = "black", linetype = "dashed") +
        scale_fill_hue(name = "Sample") +
        labs(x = "Sample", y = "SD", title = title) +
        theme_classic() +
        theme(plot.title = element_text(hjust = 0.5, size = 20), axis.text.x = element_blank(), axis.ticks.x = element_blank())

  return (p)
}

#' @author Avishay Spitzer
#'
#' @importFrom stats cor
#' @importFrom scales hue_pal
#'
#' @export
scandal_program_clusters_plot <- function(bs_scores, program_clusters) {

  prog_cor <- cor(bs_scores, use = "all.obs", method = "pearson")

  cc <- program_clusters$consensusClass
  co <- program_clusters$consensusTree$order

  p <- Heatmap(prog_cor[co, co], col = colorRamp2(c(-1, 0, 1), c("dodgerblue", "white", "red")), cluster_columns = FALSE, cluster_rows = FALSE,
               show_column_names = TRUE, show_row_names = TRUE,
               heatmap_legend_param = list(title = "Score", at = c(-1, -0.5, 0, 0.5, 1)),
               top_annotation = HeatmapAnnotation("PC" = as.character(cc)[co], show_annotation_name = FALSE,
                                                  col = list("PC" = setNames(hue_pal()(n = length(unique(cc))), unique((cc))))))

  return (p)
}

### -------------------------------------------------------------------------
### Exported plotting functions (end)
### =========================================================================
###

### =========================================================================
### Internal functions (start)
### -------------------------------------------------------------------------
###

# Convenience method for transforming an expression matrix in logTPM units to a mtrxi compatible with
# running the NMF algorithm
.nmf_matrix <- function(object) {
  m <- as.matrix(assay(object))
  m <- center_matrix(m, by = "row", method = "mean", scale = FALSE)
  m[m < 0] <- 0
  return (m)
}

#' @importFrom NMF basis
.W <- function(nmf_res, name) {
  w <- basis(nmf_res)
  colnames(w) <- paste0("P", 1:ncol(w))
  return (w)
}

.extract_programs <- function(w, n = 50) {
  programs <- list()
  for (i in 1:ncol(w))
    programs[[colnames(w)[i]]] <- head(rownames(w)[order(w[, i], decreasing = TRUE)], n = n)
  return (programs)
}

### -------------------------------------------------------------------------
### Internal functions (end)
### =========================================================================
###
