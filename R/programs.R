
### =========================================================================
### Exported algorithm functions (start)
### -------------------------------------------------------------------------
###

#'
#' @author Avishay Spitzer
#'
#' @export
scandal_metaprograms <- function(object, nmf_data, samples = NULL, n_wsp_genes = 50, n_features = 500,
                                 cc_res = NULL, maxK = 10, reps = 1000, distance = "euclidean", override_best_k = NA,
                                 n_mp_genes = 50, mpss = "ws", score_threshold = .5, mp_map = NULL, verbose = FALSE, ...) {

  stopifnot(is_scandal_object(object))
  stopifnot(!is.null(nmf_data))
  stopifnot(mpss %in% c("ws", "bs"))

  .msg("Computing meta-programs", verbose = verbose)

  if (is.null(samples))
    samples <- prepare_samples(object = object, verbose = verbose)

  ws_programs <- nmf_extract_programs(nmf_data = nmf_data, n = n_wsp_genes)

  ws_scores <- score_within_samples(samples = samples, ws_programs = ws_programs)

  # Extract the top features from all samples
  features <- .features(object = object, nmf_data = nmf_data, n_features = n_features, verbose = verbose)

  # Extract the within-sample cluster (maximal coefficient/factor) of each cell
  mcs <- .mcs(nmf_data = nmf_data, verbose = verbose)

  # Compute the log2-ratio vector for each factor within each sample
  l2r <- .l2r(object = object, samples = samples, features = features, mcs = mcs, verbose = verbose)

  # Compute the pairwise correlations between L2R vectors
  cor_l2r <- cor(x = l2r, method = "pearson")

  # Run consensus clustering on the L2R correlation matrix and find the best K (number of program clusters)
  cc <- .consensus_cluster(cor_l2r = cor_l2r, maxK = maxK, reps = reps, distance = distance, cc_res = cc_res, override_best_k = override_best_k, verbose = verbose, ...)

  # Plot the L2R correlation matrix
  .plot_l2r_cor(cor_l2r = cor_l2r, ccp = cc$ccp[[cc$best_k]], clusters = cc$clusters)

  # Compute the metaprograms by DE using the L2R vectors
  mp_l2r <- .mp_l2r(l2r = l2r, clusters = cc$clusters, verbose = verbose)

  # Take the top N genes in each metaprogram
  mps <- .metaprograms(mp_l2r = mp_l2r, n_mp_genes = n_mp_genes, verbose = verbose)

  # score all metaprograms
  mp_scores <- .score_mps(object = object, samples = samples, metaprograms = mps, mpss = mpss, verbose = verbose)

  # Assign a metaprogram to each cell
  amp <- .assign_mps(mp_scores = mp_scores, score_threshold = score_threshold, verbose = verbose)

  #as.data.frame(mps)

  # Package all the products as a ScandalMetaprograms object
  res <- ScandalMetaprograms(wsPrograms = ws_programs, wsScores = ws_scores, l2R = l2r, corL2R = cor_l2r, consensusClusters = cc,
                             mpL2R = mp_l2r, metaPrograms = mps, mpScores = mp_scores,
                             mpAssigned = amp, mpMap = mp_map, scoringStrategy = mpss, scoreThreshold = score_threshold,
                             nodeID = nodeID(object), projectID = projectID(object))

  .msg("Done", verbose = verbose)

  return (res)
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
#' @importFrom stats setNames
#'
#' @export
prepare_nmf_matrix <- function(samples) {

  mats <- setNames(lapply(samples), function(s) as.matrix(assay(s)),
                   nm = names(samples))

  nmf_mats <- lapply(names(mats), function(n) {

    m <- as.matrix(mats[[n]])
    m <- center_matrix(m, by = "row", method = "mean", scale = FALSE)
    m[m < 0] <- 0
    m
  })
  names(nmf_mats) <- names(mats)

  return (nmf_mats)
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

### -------------------------------------------------------------------------
### Exported algorithm functions (end)
### =========================================================================
###

### =========================================================================
### Internal functions (start)
### -------------------------------------------------------------------------
###

.msg <- function(str, verbose = verbose) {
  if (isTRUE(verbose))
    message(str)
}

.features <- function(object, nmf_data, n_features, verbose) {

  .msg(paste0("Pulling ", n_features, " top features from each NMF factor"), verbose = verbose)

  # Pull genes from top N features for each factor in each sample
  genes <- lapply(nmf_data, function(d) {
    fs <- NMF::extractFeatures(d, method = n_features)
    unique(unlist(lapply(fs, function(f) rownames(NMF::basis(d))[f]), recursive = FALSE))
  })

  # Pull all genes to a single vector
  all_genes <- unique(unlist(lapply(genes, function(x) x), recursive = FALSE))

  .msg(sprintf("Pulled %d features (out of %d existing features) from %d NMF objects",
               length(all_genes),
               nrow(object),
               length(nmf_data)), verbose = verbose)

  return (all_genes)
}

.mcs <- function(nmf_data, verbose) {

  .msg("Pulling maximal coefficient for each cell", verbose = verbose)

  # Pull the maximal coefficients (e.g. the selected factor) for each cell
  mcs <- lapply(nmf_data, function(d) {
    setNames(paste0("P", apply(NMF::coefficients(d), 2, which.max)), colnames(NMF::coefficients(d)))
  })

  return (mcs)
}

.l2r <- function(object, samples, features, mcs, verbose) {

  .msg("Computing L2R vector for each NMF factor (within sample)", verbose = verbose)

  # Compute the L2R (log2 ratio) vector for each factor within each sample
  l2rs <- lapply(samples, function(s) {

    m <- center_matrix(as.matrix(assay(object)[features, colnames(s)]))
    mc <- mcs[[nodeID(s)]]

    l2r <- setNames(lapply(unique(mc), function(p) {
      g1 <- names(mc[mc %in% p])
      g2 <- names(mc[!(mc %in% p)])
      rowMeans(m[, g1]) - rowMeans(m[, g2])
    }), paste0(nodeID(s), ".", unique(mc)))
    l2r <- as.matrix(as.data.frame(l2r))
    l2r
  })

  # Bind all L2R vectors as a single matrix
  l2r <- do.call(cbind, l2rs)

  return (l2r)
}

.consensus_cluster <- function(cor_l2r, maxK, reps, distance, cc_res, override_best_k, verbose, ...) {

  .msg("Clustering the L2R pairwise correlation matrix", verbose = verbose)

  if (is.null(cc_res)) {

    # Compute the consensus clustering with 1-cor as the distance object
    ccp <- ConsensusClusterPlus::ConsensusClusterPlus(d = 1 - cor_l2r, maxK = maxK, reps = reps, distance = distance, ...)
    icl <- ConsensusClusterPlus::calcICL(ccp)

    # Compute the cluster consensus and item consensus scores for each cluster
    ccs <- left_join(as_tibble(icl$clusterConsensus) %>% group_by(k) %>% summarise(CC = mean(clusterConsensus, na.rm = TRUE)),
                     as_tibble(icl$itemConsensus) %>% group_by(item, k) %>% summarise(MIC = max(itemConsensus)) %>% group_by(k) %>% summarise(IC = mean(MIC, na.rm = TRUE)),
                     by = "k")
  } else {

    .msg("Got Consensus Clustering result, skipping step", verbose = verbose)

    ccp <- cc_res$ccp
    icl <- cc_res$icl
    ccs <- cc_res$ccs
  }

  if (is.na(override_best_k)) {
    # If CC and IC scores agree on the same cluster then this is chosen as the best K (number of program clusters)
    if (which.max(ccs$CC) == which.max(ccs$IC))
      best_k <- ccs$k[which.max(ccs$CC)]
    else {
      best_k <- ccs$k[which.max(ccs$CC)]

      .msg("No agreement between Item and Cluster Consensus scores - using Cluster Consensus score to set best K. It is adviseabele to review the clustering results", verbose = verbose)
    }

    .msg(paste0("Best K set to ", best_k), verbose = verbose)

  } else {

    .msg(paste0("Overriding best K set to ", override_best_k), verbose = verbose)

    best_k <- override_best_k
  }

  clusters <- ccp[[best_k]]$consensusClass

  res <- list(ccp = ccp, icl = icl, ccs = ccs, best_k = best_k, clusters = clusters)

  return (res)
}

.plot_l2r_cor <- function(cor_l2r, ccp, clusters) {

  h <- Heatmap(cor_l2r,
               col = circlize::colorRamp2(c(-0.5, -0.25, 0, 0.25, 0.5), c("dodgerblue", "lightblue", "white", "red" ,"darkred")),
               heatmap_legend_param = list(title = "Pearson\ncorrelation", at = c(-0.5, -0.25, 0, 0.25, 0.5), labels = c("-0.5", "", "0", "", "0.5")),
               cluster_rows = ccp$consensusTree, cluster_columns = ccp$consensusTree,
               column_title = "Program pairwise correlation", column_title_side = "top", column_title_gp = gpar(fontsize = 20, fontface = "bold"),
               show_row_names = FALSE, show_column_names = FALSE, row_names_gp = gpar(fontsize = 7, fontface = "bold"),
               top_annotation = HeatmapAnnotation(Cluster = as.character(clusters),
                                                  col = list(Cluster = setNames(scales::hue_pal()(length(unique(clusters))), unique(clusters))),
                                                  show_annotation_name = FALSE))

  draw(h)

  invisible(NULL)
}

.mp_l2r <- function(l2r, clusters, verbose) {

  .msg("Computing the metaprogram L2R vectors", verbose = verbose)

  # Compute the metaprograms by DE using the L2R vectors
  mp_l2r <- setNames(lapply(unique(clusters), function(i) {

    prgs <- names(clusters[clusters == i])

    sample_l2r <- l2r[, prgs]
    ref_l2r <- l2r[, !colnames(l2r) %in% prgs]

    if (!is.null(dim(sample_l2r)))
      res <- rowMeans(sample_l2r) - rowMeans(ref_l2r)
    else
      res <- sample_l2r - rowMeans(ref_l2r)

    #res <- sort(res, decreasing = TRUE)

    return (res)
  }), paste0("MP", unique(clusters)))

  return (mp_l2r)
}

.metaprograms <- function(mp_l2r, n_mp_genes, verbose) {

  .msg(paste0("Extracting the metaprogram genes (the top ", n_mp_genes, " genes)"), verbose = verbose)
  # Take the top N genes in each metaprogram
  mps <- setNames(lapply(mp_l2r, function(x) names(head(sort(x, decreasing = TRUE), n = n_mp_genes))), names(mp_l2r))

  return(mps)
}

.score_mps <- function(object, samples, metaprograms, mpss, verbose) {

  if (mpss == "ws") {

    .msg("Scoring cells for the metaprograms within samples", verbose = verbose)

    # score all metaprograms within sample
    mp_scores <- lapply(samples, function(s) {
      m <- as.matrix(assay(object))[, colnames(s)]
      scalop::score(mat = center_matrix(m), groups = metaprograms, binmat = m, bin.control = TRUE, center = FALSE)
    })
    mp_scores <- do.call(rbind, mp_scores)
    mp_scores <- mp_scores[colnames(object), ]
  } else {

    .msg("Scoring cells for the metaprograms between samples", verbose = verbose)

    # score all metaprograms between samples
    m <- as.matrix(assay(object))

    mp_scores <- scalop::score(mat = center_matrix(m), groups = metaprograms, binmat = m, bin.control = TRUE, center = FALSE)
  }

  return (mp_scores)
}

.assign_mps <- function(mp_scores, score_threshold, verbose) {

  .msg(paste0("Assigning metaprograms with score threshold ", score_threshold), verbose = verbose)

  # Assign a metaprogram to each cell
  amp <- apply(mp_scores, 1, function(x) {
    ms <- which.max(x)
    ifelse(x[ms] >= score_threshold, paste0("MP", ms), "NA")
  })

  .msg(sprintf("Assigned MP to %d/%d of cells, NA rate %.2f",
               length(which(amp != "NA")),
                      length(amp),
                      length(which(amp == "NA")) / length(amp)), verbose = verbose)

  return (amp)
}

### -------------------------------------------------------------------------
### Internal functions (end)
### =========================================================================
###

### =========================================================================
### Old version of the algorithm, will not be exported anymore (start)
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
#'   \link{score_within_samples}. See \link{scalop::score} for more details about how to
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
scandal_programs_of_intra_sample_heterogeneity <- function(object, samples = NULL, clustering_data = NULL,
                                                           algorithm = "nmf", rank = 10, ngenes1 = 50, ngenes2 = 30, sd_threshold = 0.8, filter_method = "relative",
                                                           n_clusters_min = 2, n_clusters_max = 10,
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

  program_clusters <- cluster_variable_programs(bs_scores = bs_scores, n_clusters_min = n_clusters_min, n_clusters_max = n_clusters_max, verbose = verbose)

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
#' @importFrom scalop score
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
#' @importFrom scalop score
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
cluster_variable_programs <- function(bs_scores, n_clusters_min = 2, n_clusters_max = 10, verbose = FALSE) {

  if (isTRUE(verbose))
    message("Computing program clusters")

  ccp <- ConsensusClusterPlus(bs_scores, maxK = n_clusters_max, pItem = 1, pFeature = 1, clusterAlg = "hc", distance = "pearson")
  icl <- calcICL(ccp)

  sum_ic <- setNames(rep(0, length(n_clusters_min:n_clusters_max)), as.character(n_clusters_min:n_clusters_max))
  for (k in n_clusters_min:n_clusters_max) {
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
assign_metaprograms <- function(metaprogram_scores, verbose = FALSE) {

  max_scores <- apply(metaprogram_scores, 1, which.max)
  max_scores <- paste0("MP", max_scores)
  names(max_scores) <- rownames(metaprogram_scores)

  return(max_scores)
}

#' @author Avishay Spitzer
#'
#' @importFrom reshape2 melt
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
### Old version of the algorithm, will not be exported anymore (end)
### =========================================================================
###
